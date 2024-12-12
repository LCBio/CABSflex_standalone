"""
Module for running CABS jobs.
"""
import operator
import os
import re
import tarfile
import glob
import numpy as np
from tempfile import mkstemp
from time import strftime
import random
from functools import reduce

from abc import ABCMeta, abstractmethod
from CABS import logger, pdblib, cabs, utils
from CABS.align import save_csv, AlignError, align_to
from CABS.cluster import Clustering
from CABS.cmap import ContactMapFactory
from CABS.filter import Filter
from CABS.plots import graph_RMSF, plot_E_RMSD, plot_RMSD_N
from CABS.protein import ProteinComplex
from CABS.restraints import Restraints
from CABS.trajectory import Trajectory
from CABS.pdblib import Pdb
import CABS.optparser as opt_parser
from CABS.cmap import ContactMap
from CABS.utils import convert_cg_to_all, sync_residues

_name = 'JOB'
_CABS_files = ["TRAF", "SEQ", "INP", "OUT", "FCHAINS", "PAIRMOD"]
DEFAULT_COLORS = ['#ffffff', '#f2d600', '#4b8f24', '#666666', '#e80915', '#000000']


class CABSTask(object):
    """Abstract CABS job instance."""

    __metaclass__ = ABCMeta

    def __init__(self, **kwargs):
        self.aa_method = kwargs.get('aa_method')
        self.aa_rebuild = kwargs.get('aa_rebuild')
        self.add_peptide = kwargs.get('add_peptide')
        self.align = kwargs.get('align')
        self.align_options = dict(kwargs.get('align_options', []))
        self.align_peptide_options = dict(kwargs.get('align_peptide_options', []))
        self.binding_interactions = kwargs.get('binding_interactions')
        self.ca_rest_add = kwargs.get('ca_rest_add')
        self.ca_rest_file = kwargs.get('ca_rest_file')
        self.ca_rest_weight = kwargs.get('ca_rest_weight')
        self.category_output = kwargs.get('category_output')
        self.clustering_iterations = kwargs.get('clustering_iterations')
        self.clustering_medoids = kwargs.get('clustering_medoids')
        self.contact_map_colors = kwargs.get('contact_map_colors')
        self.contact_maps = kwargs.get('contact_maps')
        self.contact_output = kwargs.get('contact_output')
        self.contact_threshold = kwargs.get('contact_threshold')
        self.contact_threshold_aa = kwargs.get('contact_threshold_aa')
        self.cyclization = kwargs.get('backbone_cyclization')
        self.disulfide_bonds = kwargs.get('disulfide_bonds')
        self.dssp_command = kwargs.get('dssp_command')
        self.dssp_output = kwargs.get('dssp_output')
        self.exclude = kwargs.get('exclude')
        self.excluding_distance = kwargs.get('excluding_distance')
        self.filtering_count = kwargs.get('filtering_count')
        self.filtering_mode = kwargs.get('filtering_mode')
        self.fortran_command = kwargs.get('fortran_command')
        self.gauss_iterations = kwargs.get('gauss_iterations')
        self.image_file_format = kwargs.get('image_file_format')
        self.input_protein = kwargs.get('input_protein')
        self.insertion_attempts = kwargs.get('insertion_attempts')
        self.insertion_clash = kwargs.get('insertion_clash')
        self.load_cabs_files = kwargs.get('load_cabs_files')
        self.mc_annealing = kwargs.get('mc_annealing')
        self.mc_cycles = kwargs.get('mc_cycles')
        self.mc_steps = kwargs.get('mc_steps')
        self.modeller_iterations = kwargs.get('modeller_iterations')
        self.pairmod = kwargs.get('pairmod')
        self.pdb_cache = kwargs.get('pdb_cache_dir')
        self.pdb_output = kwargs.get('pdb_output')
        self.peptide = kwargs.get('peptide')
        self.plddt_output = kwargs.get('plddt_output')
        self.protein_category = kwargs.get('protein_category')
        self.protein_flexibility = kwargs.get('protein_flexibility')
        self.protein_plddt = kwargs.get('protein_plddt')
        self.protein_restraints = kwargs.get('protein_restraints')
        self.protein_restraints_reduce = kwargs.get('protein_restraints_reduce')
        self.no_protein_restraints = kwargs.get('no_protein_restraints')
        self.random_seed = kwargs.get('random_seed')
        self.receptor_ss = kwargs.get('receptor_ss')
        self.reference_pdb = kwargs.get('reference_pdb')
        self.remote = kwargs.get('log')
        self.renumber = kwargs.get('renumber_residues_to_original')
        self.replicas = kwargs.get('replicas')
        self.replicas_dtemp = kwargs.get('replicas_dtemp')
        self.restraints_output = kwargs.get('restraints_output')
        self.save_cabs_files = kwargs.get('save_cabs_files')
        self.save_config = kwargs.get('save_config')
        self.sc_rest_add = kwargs.get('sc_rest_add')
        self.sc_rest_file = kwargs.get('sc_rest_file')
        self.sc_rest_weight = kwargs.get('sc_rest_weight')
        self.separation = kwargs.get('separation')
        self.temperature = kwargs.get('temperature')
        self.verbose = kwargs.get('verbose')
        self.work_dir = kwargs.get('work_dir')
        self.weighted_fit = kwargs.get('weighted_fit')

        # Job attributes collected.
        self.config = kwargs
        self.initial_complex = None
        self.restraints = None
        self.cabsrun = None
        self.trajectory = None
        self.filtered_trajectory = None
        self.filtered_ndx = None
        self.medoids = None
        self.clusters_dict = None
        self.clusters = None
        self.rmslst = {}
        self.results = None
        self.reference = None

        # seeding RNG
        random.seed(self.random_seed)

        # Work_dir processing: making sure work_dir is abspath
        self.work_dir = os.path.abspath(self.work_dir)

        try:
            logger.setup(log_level=self.verbose, remote=self.remote, work_dir=self.work_dir,
                         save_dssp=self.dssp_output, save_restraints=self.restraints_output,
                         save_plddt=self.plddt_output, save_category=self.category_output,
                         save_contact=self.contact_output)
            os.makedirs(self.work_dir)
        except OSError:
            if os.path.isdir(self.work_dir):
                logger.warning(_name, '{} already exists. Output data will be overwritten.'.format(self.work_dir))
            else:
                logger.exit_program(
                    _name, '{} already exists and is not a directory. Choose different name.'.format(self.work_dir)
                )

        if self.dssp_command:
            pdblib.Pdb.DSSP_COMMAND = self.dssp_command

        if self.fortran_command:
            cabs.CabsRun.FORTRAN_COMMAND = self.fortran_command

        self.file_TRAF = self.file_SEQ = None
        if self.load_cabs_files:
            try:
                self.load_cabs_results()
                self.file_TRAF = os.path.join(self.work_dir, "TRAF")
                self.file_SEQ = os.path.join(self.work_dir, "SEQ")
            except (ValueError, TypeError, IOError) as e:
                logger.exit_program(
                    module_name=_name,
                    msg="Could not load CABS files from %s. An error occurred: %s" % (self.load_cabs_files, e),
                    exc=e
                )

        # self.peptide + self.add_peptide -> self.ligand
        self.peptides = []
        if self.peptide:
            self.peptides.extend([[p, 'random', 'random'] for p in self.peptide])
        if self.add_peptide:
            self.peptides.extend([p for p in self.add_peptide if p])

        valid_letters = set('RFCMSAN')

        try:
            if not all(letter in valid_letters for letter in self.pdb_output):
                raise ValueError("Contains letters outside of 'RFCMSAN'.")

            # Process 'A' or 'N' in pdb_output
            if 'A' in self.pdb_output:
                self.pdb_output = 'RFCMS'
            elif 'N' in self.pdb_output:
                self.pdb_output = ''

        except ValueError as e:
            logger.exit_program(
                module_name=_name,
                msg="Invalid pdb_output. An error occurred: %s" % e,
                exc=e
            )

        if self.contact_map_colors:
            self.colors = self.contact_map_colors
        else:
            self.colors = DEFAULT_COLORS

        # Flag to check if dynamic weights should be used
        self.gauss = self.weighted_fit == 'gauss'
        if self.gauss and self.gauss_iterations:

            utils.GAUSS_MAX_ITER = self.gauss_iterations

        # Check whether to use restraints based on pLDDT
        if not self.no_protein_restraints:
            mode, gap, min_d, max_d = self.protein_restraints
            if (
            mode in ['min', 'max', 'mean', 'plddt1', 'plddt2'] and not self.protein_plddt
            ) or (
            mode == 'category' and not (self.protein_plddt or self.protein_category)
            ):
                logger.warning(
                    _name, 'No information about pLDDT or flexibility categories provided. '
                           'Changing protein restraints  mode to \'all\'. '
                           'If you want to use restraints based on pLDDT or flexibility categories, '
                           'please provide the necessary data.')
                self.protein_restraints = ('all', gap, min_d, max_d)

        # pairwise potential modification
        if self.pairmod:
            pairmod = {}
            try:
                with open(self.pairmod, 'r') as f:
                    for line in f:
                        r, w, s = line.replace('PEP ', 'PEP1 ').split()
                        pairmod[r] = (float(w), float(s))
                    logger.info(
                        module_name=_name,
                        msg='Using pairwise potential modification from ' + self.pairmod
                    )
                    self.pairmod = pairmod
            except (IOError, ValueError):
                logger.warning(_name, 'Error while reading file: ' + self.pairmod)
                self.pairmod = None

    def run(self):
        file_traf = self.file_TRAF
        file_seq = self.file_SEQ
        self.setup_job()
        with_cabs = None in (file_traf, file_seq)

        if with_cabs:
            self.setup_cabs_run()
            self.execute_cabs_run()
        if self.save_cabs_files:
            self.save_cabs_res()
        self.load_output(file_traf, file_seq)
        if self.reference_pdb:
            self.parse_reference(self.reference_pdb, self.pdb_cache)
        self.score_results(
            n_filtered=self.filtering_count,
            number_of_medoids=self.clustering_medoids,
            number_of_iterations=self.clustering_iterations
        )

        self.save_models()
        if self.reference:
            try:
                self.calculate_rmsd()
            except (ValueError, AlignError) as e:
                logger.critical(module_name=_name, msg=e.message)
        self.save_config_file()
        self.load_all_atom_tops()
        self.draw_plots(colors=self.colors)
        if self.load_cabs_files:
            for _file in _CABS_files:
                try:
                    os.remove(os.path.join(self.work_dir, _file))
                except OSError:
                    pass
        logger.info(module_name=_name, msg='Simulation completed successfully')

    def save_cabs_res(self):
        tar_dir = mkstemp(
            prefix=strftime(cabs.CabsRun.CABS_DIR_FMT),
            dir=self.work_dir,
            suffix='.cbs'
        )[1]
        with tarfile.open(tar_dir, "w:gz") as tar:
            logger.log_file(
                _name, "Saving CABS simulation files to: %s" % tar_dir)
            for file_name in _CABS_files:
                try:
                    tar.add(os.path.join(
                        self.cabsrun.cfg['cwd'], file_name), arcname=file_name)
                except OSError:
                    pass
                    # print("Could not add %s to the tar file" % file_name)

    def load_cabs_results(self):
        if not os.path.exists(self.load_cabs_files):
            logger.exit_program(
                module_name=_name,
                msg="Provided CABS files path does not exist (%s)" % self.load_cabs_files,
                traceback=False
            )
        try:
            files = glob.glob(os.path.join(self.load_cabs_files, "*.cbs"))
            if len(files) > 1:
                logger.critical(module_name=_name,
                                msg="More than one .cbs file in provided directory %s "
                                    % " \n".join(files))
                logger.exit_program(module_name=_name,
                                    msg="Please re-run with --load-cabs-files <filename> "
                                        "or remove the files you do not need. Quiting.",
                                    traceback=False)
            elif len(files) == 1:
                logger.info(module_name=_name,
                            msg="Loading CABS files from %s" % files[0])
                with tarfile.open(files[0], "r:gz") as f:
                    f.extractall(os.path.join(self.work_dir))

            else:
                raise IOError

        except IOError:
            files_loc = self.load_cabs_files
            try:
                with tarfile.open(files_loc, "r:gz") as f:
                    logger.info(module_name=_name,
                                msg="Loading CABS files from %s" % files_loc)
                    f.extractall(os.path.join(self.work_dir))
            except IOError:
                raise
        return

    @abstractmethod
    def setup_job(self):
        pass

    @abstractmethod
    def calculate_rmsd(self):
        pass

    @abstractmethod
    def parse_reference(self, ref, pdb_cache):
        mtx_q, mtx_p, dummy_aln = align_to(
            self.reference[0],
            self.reference[1],
            self.initial_complex.protein,
            self.initial_complex.protein.list_chains().keys(),
            self.align,
            self.align_options
        )
        mtx_p = mtx_p.to_numpy()
        mtx_q = mtx_q.to_numpy()
        dummy_rmsd, rot, t_com, q_com = utils.dynamic_kabsch(mtx_p, mtx_q)
        self.reference[0].from_np(np.dot(self.reference[0].to_numpy() - q_com, rot) + t_com)

    def load_all_atom_tops(self):
        if self.aa_rebuild:
            pth = os.path.join(self.work_dir, "output_pdbs", "model_%i.pdb")
            med_traj = np.array([Pdb(pth % i).atoms.to_numpy() for i in range(10)])
            med_traj = np.expand_dims(med_traj, axis=1)
            mod0 = Pdb(os.path.join(self.work_dir, "output_pdbs", "model_0.pdb"))
            self.medoids = Trajectory(mod0, med_traj, self.medoids.headers)
        else:
            self.medoids.coordinates = np.swapaxes(self.medoids.coordinates, 0, 1)

    @abstractmethod
    def draw_plots(self, plots_dir=None, colors=DEFAULT_COLORS):
        pass

    @abstractmethod
    def mk_cmaps(self, ca_traj, meds, clusts, top1k_inds, thr, thra, plots_dir):
        """
        Arguments:
        ca_traj -- np.array of coordinates in shape (replicas, steps, atoms, 3).
        meds -- np.array as ca_traj representing medoids in shape (n_of_meds, 1, atoms, 3).
        clusts -- np.array of clusters in shape (n_of_clusts, ?, atoms, 3).
        top1k_inds -- indices of 1k top scored models in ca_traj.
        thr -- atom distance cutoff for side chain distance contact criterion.
        thra -- atom distance cutoff for heavy atoms distance contact criterion.
        plots_dir -- path to directory storing plots.
        """
        scmodeler = utils.SCModeler(ca_traj.template)
        sc_traj_full = scmodeler.calculate_sc_traj(ca_traj.coordinates)

        if self.aa_rebuild:
            sc_med = meds.coordinates
        else:
            sc_med = scmodeler.calculate_sc_traj(meds.coordinates)

        cmapdir = os.path.join(self.work_dir, 'contact_maps')
        try:
            os.mkdir(cmapdir)
        except OSError:
            pass

        return sc_traj_full, sc_med, cmapdir

    def prepare_restraints(self, output=''):

        # generate protein restraints
        if self.no_protein_restraints:
            protein_restraints = Restraints(None)
        else:
            protein_restraints = Restraints(
                self.initial_complex.protein.generate_restraints(
                    *self.protein_restraints)
            )

        # reduce number of restraints
        if self.protein_restraints_reduce:
            protein_restraints.reduce_by(self.protein_restraints_reduce)

        # additional restraints
        add_restraints = Restraints('')

        if self.ca_rest_add:
            add_restraints += Restraints.from_parser(self.ca_rest_add)

        if self.sc_rest_add:
            add_restraints += Restraints.from_parser(self.sc_rest_add, sg=True)

        if self.ca_rest_file:
            for filename in self.ca_rest_file:
                add_restraints += Restraints.from_file(filename)

        if self.sc_rest_file:
            for filename in self.sc_rest_file:
                add_restraints += Restraints.from_file(filename, sg=True)

        if self.cyclization:
            add_restraints += Restraints(
                self.initial_complex.protein.generate_backbone_restraints(
                    self.cyclization)
            )

        if self.disulfide_bonds:
            add_restraints += Restraints(
                self.initial_complex.protein.generate_disulfide_restraints(
                    self.disulfide_bonds), sg=True
            )

        protein_restraints += add_restraints.update_id(
            self.initial_complex.new_ids)

        if output and logger.output_restraints():
            output_restraints = os.path.join(output, 'output_data', 'restraints.txt')
            d = os.path.dirname(output_restraints)
            if not os.path.isdir(d):
                os.makedirs(d)
            logger.to_file(
                filename=output_restraints,
                content=str(protein_restraints),
                msg='Saving restraints output to %s' % output_restraints
            )

        return protein_restraints

    def save_config_file(self):
        if self.save_config:
            with open(os.path.join(self.work_dir, 'config.ini'), 'w') as configfile:
                configfile.write(utils.CONFIG_HEADER)
                for k in sorted(self.config):
                    value = self.config[k]
                    name = re.sub("_", "-", str(k))
                    option = opt_parser.option_formatter(name, value)
                    try:
                        configfile.write(option)
                    except Exception as e:
                        logger.warning(
                            _name, "Failed to save %s option to config file. Reason: %s." % (
                                name, e)
                        )

    def setup_cabs_run(self):
        logger.info(module_name="CABS", msg='Setting up CABS simulation.')
        # Initializing CabsRun instance
        self.cabsrun = cabs.CabsRun(
            protein_complex=self.initial_complex,
            restraints=self.prepare_restraints(output=self.work_dir),
            work_dir=self.work_dir,
            excluding_distance=self.excluding_distance,
            replicas=self.replicas,
            replicas_dtemp=self.replicas_dtemp,
            binding_interactions=self.binding_interactions,
            temperature=self.temperature,
            ca_rest_weight=self.ca_rest_weight,
            sc_rest_weight=self.sc_rest_weight,
            mc_annealing=self.mc_annealing,
            mc_cycles=self.mc_cycles,
            mc_steps=self.mc_steps,
            pairmod=self.pairmod
        )
        return self.cabsrun

    def execute_cabs_run(self):
        self.cabsrun.run()

    @abstractmethod
    def load_output(self, ftraf=None, fseq=None):
        """
        Method for loading previously done simulation results. Stores the results to self.trajectory.
        :param ftraf: path to TRAF file
        :param fseq: path to SEQ file
        :return: returns trajectory.Trajectory instance
        """
        if ftraf is not None and fseq is not None:
            logger.debug(
                module_name=_name, msg="Loading trajectories from: %s, %s" % (ftraf, fseq))
            self.trajectory = Trajectory.read_trajectory(ftraf, fseq)
        else:
            logger.debug(module_name=_name, msg="Loading trajectories from the CABS run")
            self.trajectory = self.cabsrun.get_trajectory()
        self.trajectory.weights = self.initial_complex.protein.weights
        self.trajectory.template.update_ids(self.initial_complex.protein.old_ids, pedantic=False)
        self.initial_complex.protein.update_ids(self.initial_complex.protein.old_ids)
        chs = ''.join(self.initial_complex.protein_chains)
        tchs = ''.join(sorted(set(chs).intersection(self.trajectory.template.list_chains())))
        self.trajectory.tmp_target_chs = tchs
        ic_stc, tt_stc, dummy_aln = self.trajectory.align_to(
            self.initial_complex.protein, chs, tchs, align_mth='trivial'
        )
        self.trajectory.superimpose_to(ic_stc, tt_stc)
        logger.info(module_name=_name, msg="Trajectories loaded successfully")
        return self.trajectory

    @abstractmethod
    def score_results(self, n_filtered, number_of_medoids,
                      number_of_iterations):
        pass

    def save_models(self):
        # output folder
        output_folder = os.path.join(self.work_dir, 'output_pdbs')
        logger.log_file(module_name=_name,
                        msg="Saving pdb files to " + str(output_folder))
        try:
            os.mkdir(output_folder)
        except OSError:
            logger.warning(
                _name, "Possibly overwriting previous pdb files. Use --work-dir <DIR> to avoid that.")
        # Saving the trajectory to PDBs:
        if 'R' in self.pdb_output:
            logger.log_file(module_name=_name, msg='Saving replicas...')
            self.trajectory.to_pdb(mode='replicas', to_dir=output_folder)
        # Saving top1000 models to PDB:
        if 'F' in self.pdb_output:
            logger.log_file(module_name=_name, msg='Saving filtered models...')
            self.filtered_trajectory.to_pdb(
                mode='replicas', to_dir=output_folder, name='top1000')
        # Saving clusters in CA representation
        if 'C' in self.pdb_output:
            logger.log_file(module_name=_name, msg='Saving clusters...')
            for i, cluster in enumerate(self.clusters):
                cluster.to_pdb(mode='replicas', to_dir=output_folder,
                               name='cluster_{0}'.format(i))
        if 'S' in self.pdb_output:
            logger.log_file(module_name=_name,
                            msg='Saving starting structure...')
            self.initial_complex.save_to_pdb(
                os.path.join(output_folder, 'start.pdb'))

        # Saving final models:
        if 'M' in self.pdb_output:
            if self.aa_rebuild:
                if self.aa_method == 'modeller':
                    logger.log_file(
                        module_name=_name, msg='Saving final models (in AA representation)')
                    pdb_medoids = self.medoids.to_pdb()
                    from CABS.ca2all import ca2all
                    from CABS.pdblib import Pdb
                    logger.log_file(
                        module_name=_name, msg='Running Modeller to rebuild models')
                    for i, fname in enumerate(pdb_medoids):
                        ca2all(
                            fname,
                            output=os.path.join(
                                output_folder, 'model_{0}.pdb'.format(i)),
                            iterations=self.modeller_iterations,
                            out_mdl=os.path.join(
                                self.work_dir, 'output_data', 'modeller_output_{0}.txt'.format(i)),
                            work_dir=self.work_dir
                        )
                        pth_tmp = os.path.join(
                            self.work_dir, 'output_pdbs', 'model_{0}.pdb'.format(i))
                        mod = Pdb(pth_tmp)
                        ssh = mod.mk_ss_header()
                        mod.atoms.save_to_pdb(pth_tmp, header=ssh)
                elif self.aa_method == 'cg2all':
                    logger.log_file(
                        module_name=_name, msg='Saving final models (in AA representation)')
                    pdb_medoids = self.medoids.to_pdb()
                    logger.log_file(
                        module_name=_name, msg='Running cg2all to rebuild models')
                    for i, fname in enumerate(pdb_medoids):
                        convert_cg_to_all(fname, work_dir=self.work_dir,
                                          iter=i,
                                          reference_pdb=self.input_protein,
                                          renumber_flag=self.renumber)

                else:
                    logger.warning(
                        module_name=_name, msg='Unknown AA method: %s' % self.aa_method)
                    logger.log_file(
                        module_name=_name, msg='Saving final models (in CA representation)')
                    self.medoids.to_pdb(
                        mode='models', to_dir=output_folder, name='model')
            else:
                logger.log_file(
                    module_name=_name, msg='Saving final models (in CA representation)')
                self.medoids.to_pdb(
                    mode='models', to_dir=output_folder, name='model')

            if self.renumber:
                logger.log_file(
                    module_name=_name, msg='Renumbering residues to original numbering')
                pass


class DockTask(CABSTask):
    """Class representing single CABS job."""

    def setup_job(self):
        if not self.peptides:
            raise ValueError('No peptide given')
        self.initial_complex = ProteinComplex(
            protein=self.input_protein,
            flexibility=self.protein_flexibility,
            exclude=self.exclude,
            weights=self.weighted_fit,
            plddt=self.protein_plddt,
            category=self.protein_category,
            peptides=self.peptides,
            replicas=self.replicas,
            separation=self.separation,
            insertion_attempts=self.insertion_attempts,
            insertion_clash=self.insertion_clash,
            work_dir=self.work_dir,
            receptor_ss=self.receptor_ss,
            pdb_cache=self.pdb_cache,
        )

    def load_output(self, ftraf=None, fseq=None):
        """
        Method for loading previously done simulation results. Stores the results to self.trajectory.
        :param ftraf: path to TRAF file
        :param fseq: path to SEQ file
        :return: returns trajectory.Trajectory instance
        """
        ret = super(DockTask, self).load_output(ftraf, fseq)
        ret.number_of_peptides = len(self.peptides)
        return ret

    def calculate_rmsd(self, save=True):
        logger.debug(module_name=_name, msg="RMSD calculations starting...")
        sfname = None
        if save:
            odir = os.path.join(self.work_dir, 'output_data')
            try:
                os.mkdir(odir)
            except OSError:
                pass
        all_results = {}
        ref_trg_stc, self_trg_stc, trg_aln = self.trajectory.align_to(
            self.reference[0], self.reference[1], self.trajectory.tmp_target_chs,
            align_mth=self.align, kwargs=self.align_options
        )
        if save:
            sfname = os.path.join(
                self.work_dir, 'output_data', 'reference_alignment')
            paln_trg = sfname + '_target.csv'
            save_csv(paln_trg, ('reference', 'template'), trg_aln)
        for pept_chain, ref_pept_chain in zip(self.initial_complex.peptide_chains, self.reference[2]):
            ref_pep_stc, self_pep_stc, pep_aln = self.trajectory.align_to(
                self.reference[0], ref_pept_chain, pept_chain, align_mth=self.align,
                kwargs=self.align_peptide_options
            )
            if save:
                paln_pep = sfname + '_%s.csv' % pept_chain
                save_csv(paln_pep, ('reference', 'template'), pep_aln)
            self.rmslst[pept_chain] = self.trajectory.rmsd_to_reference(
                ref_pep_stc, self_pep_stc)
            rmsds = [header.rmsd for header in self.medoids.headers]
            results = {
                'rmsds_all': [header.rmsd for header in self.trajectory.headers],
                'rmsds_filtered': [header.rmsd for header in self.filtered_trajectory.headers],
                'rmsds_medoids': rmsds
            }
            results['lowest_all'] = sorted(results['rmsds_all'])[0]
            results['lowest_filtered'] = sorted(results['rmsds_filtered'])[0]
            results['lowest_medoids'] = sorted(results['rmsds_medoids'])[0]
            # Saving rmsd results
            if save:
                with open(os.path.join(odir, 'lowest_rmsds_%s.txt' % pept_chain), 'w') as outfile:
                    outfile.write(
                        'lowest_all; lowest_filtered; lowest_medoids\n {0};{1};{2}'.format(
                            results['lowest_all'], results['lowest_filtered'], results['lowest_medoids']
                        )
                    )
                for _type in ['all', 'filtered', 'medoids']:
                    with open(os.path.join(odir, '{0}_rmsds_{1}.txt'.format(_type, pept_chain)), 'w') as outfile:
                        for rmsd in results['rmsds_' + _type]:
                            outfile.write(str(rmsd) + ';\n')
            all_results[pept_chain] = results
        logger.info(module_name=_name, msg="RMSD successfully saved")
        return all_results

    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        logger.debug(module_name=_name, msg="Scoring results")
        # Filtering the trajectory
        self.filtered_trajectory, self.filtered_ndx = Filter(
            self.trajectory, n_filtered).cabs_filter()
        # Clustering the trajectory
        self.medoids, self.clusters_dict, self.clusters = Clustering(
            self.filtered_trajectory,
            'chain ' + ','.join(
                self.initial_complex.peptide_chains,
            )
        ).cabs_clustering(number_of_medoids=number_of_medoids, number_of_iterations=number_of_iterations)
        logger.info(module_name=_name, msg="Scoring results successful")

    def draw_plots(self, plots_dir=None, colors=None):
        logger.debug(module_name=_name, msg="Drawing plots")
        super(DockTask, self).draw_plots()
        # set the plots dir
        if plots_dir is None:
            pltdir = os.path.join(self.work_dir, 'plots')
            try:
                os.mkdir(pltdir)
            except OSError:
                pass
        else:
            pltdir = plots_dir
        logger.log_file(module_name=_name, msg="Saving plots to %s" % pltdir)

        graph_RMSF(self.trajectory, self.initial_complex.protein_chains,
                   os.path.join(pltdir, 'RMSF'))

        # RMSD-based graphs
        if self.reference_pdb:
            logger.log_file(module_name=_name, msg="Saving RMSD plots")
            for k, rmslst in self.rmslst.items():
                plot_E_RMSD(
                    [self.trajectory, self.filtered_trajectory],
                    [rmslst, rmslst[self.filtered_ndx, ]],
                    ['all models', 'top 1000 models'],
                    os.path.join(pltdir, 'E_RMSD_%s' % k)
                )
                plot_RMSD_N(
                    rmslst.reshape(self.replicas, -1),
                    os.path.join(pltdir, 'RMSD_frame_%s' % k)
                )

        # Contact maps
        if self.contact_maps:
            logger.log_file(module_name=_name, msg="Saving contact maps")
            self.mk_cmaps(
                self.trajectory,
                self.medoids,
                self.clusters_dict,
                self.filtered_ndx,
                self.contact_threshold,
                self.contact_threshold_aa,
                pltdir,
                colors=colors
            )
        logger.info(module_name=_name, msg="Plots successfully saved")

    @staticmethod
    def _add_cmaps(mk_cmap_output):
        # breakpoint()
        map_1, map_2 = mk_cmap_output
        return ContactMap(map_1.cmtx + map_2.cmtx, map_1.s1, map_2.s2, map_1.n + map_2.n)

    def mk_cmaps(self, ca_traj, meds, clusts, top1k_inds, thr, thra, plots_dir, colors=DEFAULT_COLORS):
        sc_traj_full, sc_med, cmapdir = super(DockTask, self).mk_cmaps(
            ca_traj, meds, clusts, top1k_inds, thr, thra, plots_dir
        )

        thrt = thra if self.aa_rebuild else thr

        sc_traj_1k = sc_traj_full.reshape(
            1, -1, len(ca_traj.template), 3)[:, top1k_inds, :, :]

        rchs = self.initial_complex.protein_chains
        lchs = self.initial_complex.peptide_chains

        targ_cmf = ContactMapFactory(rchs, rchs, ca_traj.template)

        cmfs = {lig: ContactMapFactory(
            rchs, lig, ca_traj.template) for lig in lchs}
        # cmap10ktarg = self._add_cmaps(targ_cmf.mk_cmap(sc_traj_full, thr))
        cmap10ktarg = reduce(operator.add, targ_cmf.mk_cmap(sc_traj_full, thr))
        cmap10ktarg.zero_diagonal()
        cmap10ktarg.save_all(cmapdir + '/target_all',
                             break_long_x=0, norm_n=True, colors=colors)

        for lig, cmf in cmfs.items():
            cmaps = cmf.mk_cmap(sc_traj_full, thr)
            for n, cmap in enumerate(cmaps):
                cmap.save_all(cmapdir + '/replica_%i_ch_%s' %
                              (n + 1, lig), norm_n=True, colors=colors)
            cmap10k = reduce(operator.add, cmaps)
            cmap10k.save_all(cmapdir + '/all_ch_%s' %
                             lig, norm_n=True, colors=colors)
            cmap10k.save_histo(plots_dir + '/all_contacts_histo_%s' % lig)
            cmap1k = cmf.mk_cmap(sc_traj_1k, thr)[0]
            cmap1k.save_all(cmapdir + '/top1000_ch_%s' %
                            lig, norm_n=True, colors=colors)
            if self.aa_rebuild:
                cmft = ContactMapFactory(rchs, rchs, self.medoids.template)
                #  ref to self.medoids is confusing -- remove in future
            else:
                cmft = cmf
            cmaps_top = cmft.mk_cmap(sc_med, thrt)
            for n, cmap in enumerate(cmaps_top):
                cmap.save_all(cmapdir + '/top_%i_ch_%s' %
                              (n + 1, lig), norm_n=True, colors=colors)
            for cn, clust in clusts.items():
                ccmap = cmf.mk_cmap(sc_traj_1k, thr, frames=clust)[0]
                ccmap.save_all(cmapdir + '/cluster_%i_ch_%s' %
                               (cn, lig), norm_n=True, colors=colors)

    def parse_reference(self, ref, pdb_cache):
        try:
            source, rec, pep = ref.split(':')
            self.reference = (pdblib.Pdb(
                ref, pdb_cache=pdb_cache, selection='name CA and (chain %s)' % ','.join(rec + pep),
                no_exit=True, verify=True
            ).atoms, rec, pep)
            super(DockTask, self).parse_reference(ref, pdb_cache)
            if len(self.initial_complex.peptide_chains) != len(self.reference[2]):
                raise ValueError
            logger.info(_name, 'Reference {} loaded.'.format(ref))
        except (ValueError, pdblib.Pdb.InvalidPdbInput):
            logger.warning(_name, 'Invalid reference {}'.format(ref))
            self.reference = None


class FlexTask(CABSTask):
    """Class of CABSFlex jobs."""

    def setup_job(self):
        self.initial_complex = ProteinComplex(
            protein=self.input_protein,
            flexibility=self.protein_flexibility,
            exclude=self.exclude,
            weights=self.weighted_fit,
            plddt=self.protein_plddt,
            category=self.protein_category,
            peptides=self.peptides,
            replicas=self.replicas,
            separation=self.separation,
            insertion_attempts=self.insertion_attempts,
            insertion_clash=self.insertion_clash,
            work_dir=self.work_dir,
            receptor_ss=self.receptor_ss,
            pdb_cache=self.pdb_cache,
        )

        if self.reference_pdb is None:
            self.reference_pdb = True

        # remove filtered trajectory from pdb saving
        self.pdb_output = self.pdb_output.replace('F', '')

    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        # Clustering the trajectory
        clst = Clustering(self.trajectory, 'chain ' +
                          ','.join(self.initial_complex.protein_chains))
        self.medoids, self.clusters_dict, self.clusters = clst.cabs_clustering(
            number_of_medoids=number_of_medoids,
            number_of_iterations=number_of_iterations
        )
        self.rmslst = {
            self.initial_complex.protein_chains: clst.distance_matrix[0]}

    def load_output(self, *args, **kwargs):
        ret = super(FlexTask, self).load_output(*args, **kwargs)
        ret.number_of_peptides = 0
        return ret

    def calculate_rmsd(self, reference_pdb=None, save=True):
        logger.debug(module_name=_name, msg="RMSD calculations starting...")
        odir = None
        if save:
            odir = os.path.join(self.work_dir, 'output_data')
            try:
                os.mkdir(odir)
            except OSError:
                pass

        chs_ids = self.trajectory.tmp_target_chs
        ref_trg_stc, self_trg_stc, trg_aln = self.trajectory.align_to(
            self.reference[0], self.reference[1], chs_ids,
            align_mth=self.align, kwargs=self.align_options
        )
        if save:
            sfname = os.path.join(
                self.work_dir, 'output_data', 'reference_alignment')
            paln_trg = sfname + '_target.csv'
            save_csv(paln_trg, ('reference', 'template'), trg_aln)
        self.rmslst[chs_ids] = self.trajectory.rmsd_to_reference(
            ref_trg_stc, self_trg_stc)
        rmsds = [header.rmsd for header in self.medoids.headers]
        results = {
            'rmsds_all': [header.rmsd for header in self.trajectory.headers],
            'rmsds_medoids': rmsds
        }
        results['lowest_all'] = sorted(results['rmsds_all'])[0]
        results['lowest_medoids'] = sorted(results['rmsds_medoids'])[0]
        # Saving rmsd results
        if save:
            with open(os.path.join(odir, 'lowest_rmsds_%s.txt' % chs_ids), 'w') as outfile:
                outfile.write(
                    'lowest_all; lowest_medoids\n {0};{1}'.format(
                        results['lowest_all'], results['lowest_medoids']
                    )
                )
            for _type in ['all', 'medoids']:
                with open(os.path.join(odir, '{0}_rmsds_{1}.txt'.format(_type, chs_ids)), 'w') as outfile:
                    for rmsd in results['rmsds_' + _type]:
                        outfile.write(str(rmsd) + ';\n')
        logger.info(module_name=_name, msg="RMSD successfully saved")
        return {chs_ids: results}

    def mk_cmaps(self, ca_traj, meds, clusts, top1k_inds, thr, thra, plots_dir, colors=DEFAULT_COLORS):
        sc_traj_full, sc_med, cmapdir = super(FlexTask, self).mk_cmaps(
            ca_traj, meds, clusts, top1k_inds, thr, thra, plots_dir
        )
        # Here add saving the contact maps

        thrt = thra if self.aa_rebuild else thr

        rchs = self.initial_complex.protein_chains

        cmf = ContactMapFactory(rchs, rchs, ca_traj.template)
        if self.aa_rebuild:
            cmft = ContactMapFactory(rchs, rchs, self.medoids.template)
            #  ref to self.medoids is confusing -- remove in future
        else:
            cmft = cmf
        cmap_all = reduce(operator.add, cmf.mk_cmap(sc_traj_full, thr))

        topscms = cmft.mk_cmap(sc_med, thrt)
        cmaptop = reduce(operator.add, topscms)

        for cmap, fname in zip(
                (cmap_all, cmaptop) + tuple(topscms),
                ('all', 'top10') + tuple(["top_%i" % (i + 1) for i, dummy in enumerate(topscms)])):
            cmap.zero_diagonal()
            cmap.save_all(cmapdir + '/' + fname, break_long_x=0,
                          norm_n=True, colors=colors)

    def parse_reference(self, ref, pdb_cache):
        try:
            try:
                dummy, trg_chids = ref.split(":")
                self.reference = (
                    pdblib.Pdb(ref, pdb_cache=pdb_cache, selection='name CA',
                               no_exit=True, verify=True).atoms, trg_chids
                )
                super(FlexTask, self).parse_reference(ref, pdb_cache)
                logger.info(_name, 'Reference {} loaded.'.format(ref))
            except AttributeError:  # if ref is None it has no split mth
                ref_stc = self.initial_complex.select(
                    'name CA and (chain %s)' % ','.join(self.initial_complex.protein_chains)
                )
                ref_stc.update_ids(self.initial_complex.protein.old_ids)
                self.reference = (ref_stc, "".join(sorted(set([i.chid for i in ref_stc]))))
                logger.info(_name, 'Input loaded as reference.')
        except (pdblib.Pdb.InvalidPdbInput, ValueError):
            logger.warning(_name, 'Invalid reference {}'.format(ref))

    def draw_plots(self, plots_dir=None, colors=DEFAULT_COLORS):
        super(FlexTask, self).draw_plots()
        # set the plots dir
        if plots_dir is None:
            pltdir = os.path.join(self.work_dir, 'plots')
            try:
                os.mkdir(pltdir)
            except OSError:
                pass
        else:
            pltdir = plots_dir

        graph_RMSF(self.trajectory,
                   self.initial_complex.protein_chains,
                   os.path.join(pltdir, 'RMSF'),
                   fmt=self.image_file_format)

        # RMSD-based graphs
        if self.reference_pdb:
            for k, rmslst in self.rmslst.items():
                plot_E_RMSD(
                    [self.trajectory],
                    [rmslst],
                    ['all models'],
                    os.path.join(pltdir, 'E_RMSD_%s' % k),
                    self.image_file_format,
                    interaction=False
                )
                plot_RMSD_N(
                    rmslst.reshape(self.replicas, -1),
                    os.path.join(pltdir, 'RMSD_frame_%s' % k),
                    self.image_file_format
                )

        # Contact maps
        if self.contact_maps:
            self.mk_cmaps(
                self.trajectory,
                self.medoids,
                self.clusters_dict,
                self.filtered_ndx,
                self.contact_threshold,
                self.contact_threshold_aa,
                pltdir,
                colors=colors
            )
