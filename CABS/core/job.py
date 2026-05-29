"""
Module for running CABS jobs with comprehensive type annotations.
"""

from abc import ABCMeta, abstractmethod
import csv
from copy import deepcopy
from functools import reduce
import glob
import operator
import os
import random
import re
import tarfile
from tempfile import NamedTemporaryFile
from time import strftime
from typing import Any, Dict, List, Literal, Optional, Union

import numpy as np
import numpy.typing as npt

from CABS.analysis.cluster import Clustering
from CABS.analysis.cmap import ContactMap, ContactMapFactory
from CABS.analysis.plots import drop_csv_file, graph_RMSF, plot_E_RMSD, plot_RMSD_N
from CABS.analysis.restraints import Restraints
from CABS.constants import (
    ALLOWED_AA_METHODS,
    CABS_FILES,
    CONFIG_HEADER,
    DEFAULT_COLORS,
    ColorHex,
)
from CABS.core import cabs
from CABS.core.trajectory import Trajectory
from CABS.io import logger
import CABS.io.optparser as opt_parser
from CABS.structures import pdblib, protein
from CABS.structures.atom import Atoms
from CABS.structures.pdblib import Pdb
from CABS.structures.protein import ProteinComplex
from CABS.utils import utils
from CABS.utils.align import AlignError, align_to, save_csv
from CABS.utils.filter import Filter
from CABS.reconstruction.cg2all import convert_cg_to_all
from CABS.reconstruction.cg2all_trajectory import reconstruct_job_outputs
from CABS.utils import chimerax as chimerax_utils
from CABS.utils import notebook as notebook_utils
from CABS.utils import pymol as pymol_utils
from CABS.config_loader import get_config_section, get_cg2all_env_prefix

_name = "JOB"


class CABSTask(metaclass=ABCMeta):
    """Abstract CABS job instance with comprehensive type annotations."""

    def __init__(self, **kwargs: Any) -> None:
        # Type annotations for all attributes
        self.aa_method: Optional[Literal["modeller", "cg2all"]] = kwargs.get(
            "aa_method"
        )
        self.aa_rebuild: str = str(kwargs.get("aa_rebuild", "M"))
        self.aa_rebuild_workers: Optional[int] = kwargs.get("aa_rebuild_workers")
        self.aa_minimize: Optional[bool] = kwargs.get("aa_minimize")
        self.add_peptide: Optional[str] = kwargs.get("add_peptide")
        self.align: Optional[bool] = kwargs.get("align")
        self.align_options: Dict[str, Any] = dict(kwargs.get("align_options", []))
        self.align_peptide_options: Dict[str, Any] = dict(
            kwargs.get("align_peptide_options", [])
        )
        self.binding_interactions: Optional[bool] = kwargs.get("binding_interactions")
        self.ca_rest_add: Optional[str] = kwargs.get("ca_rest_add")
        self.ca_rest_file: Optional[str] = kwargs.get("ca_rest_file")
        self.ca_rest_weight: Optional[float] = kwargs.get("ca_rest_weight")
        self.cg2all_env_prefix: Optional[str] = kwargs.get("cg2all_env_prefix")
        self.cg2all_representation: str = kwargs.get(
            "cg2all_representation", "calpha-sc"
        )
        self.clustering_iterations: Optional[int] = kwargs.get("clustering_iterations")
        self.clustering_medoids: Optional[int] = kwargs.get("clustering_medoids")
        self.contact_map_colors: Optional[List[ColorHex]] = kwargs.get(
            "contact_map_colors"
        )
        self.contact_maps: Optional[bool] = kwargs.get("contact_maps")
        self.contact_threshold: Optional[float] = kwargs.get("contact_threshold")
        self.contact_threshold_aa: Optional[float] = kwargs.get("contact_threshold_aa")
        self.csv_output: Optional[bool] = kwargs.get("csv_output")
        self.cyclization: Optional[bool] = kwargs.get("backbone_cyclization")
        self.disable_centro: Optional[bool] = kwargs.get("disable_centro")
        self.disulfide_bonds: Optional[bool] = kwargs.get("disulfide_bonds")
        self.dssp_output: Optional[bool] = kwargs.get("dssp_output")
        self.exclude: Optional[Union[str, List[str]]] = kwargs.get("exclude")
        self.excluding_distance: Optional[float] = kwargs.get("excluding_distance")
        self.filtering_count: Optional[int] = kwargs.get("filtering_count")
        self.filtering_mode: Optional[str] = kwargs.get("filtering_mode")
        self.fortran_command: Optional[str] = kwargs.get("fortran_command")
        self.gauss_iterations: Optional[int] = kwargs.get("gauss_iterations")
        self.generate_pymol_visualizations: Optional[bool] = kwargs.get(
            "generate_pymol_visualizations"
        )
        self.generate_chimera_visualizations: Optional[bool] = kwargs.get(
            "generate_chimera_visualizations"
        )
        self.generate_notebook: Optional[bool] = kwargs.get("generate_notebook")
        self.image_file_format: Optional[str] = kwargs.get("image_file_format")
        self.input_protein: Optional[str] = kwargs.get("input_protein")
        self.insertion_attempts: Optional[int] = kwargs.get("insertion_attempts")
        self.insertion_clash: Optional[float] = kwargs.get("insertion_clash")
        self.json_output: Optional[bool] = kwargs.get("json_output")
        self.load_cabs_files: Optional[str] = kwargs.get("load_cabs_files")
        self.mc_annealing: Optional[bool] = kwargs.get("mc_annealing")
        self.mc_cycles: Optional[int] = kwargs.get("mc_cycles")
        self.mc_steps: Optional[int] = kwargs.get("mc_steps")
        self.modeller_iterations: Optional[int] = kwargs.get("modeller_iterations")
        self.nsp3_model_path: Optional[str] = kwargs.get("nsp3_model_path")
        self.pairmod: Optional[str] = kwargs.get("pairmod")
        self.pdb_cache: Optional[str] = kwargs.get("pdb_cache_dir")
        self.pdb_bfac_output: Optional[str] = kwargs.get("pdb_bfac_output")
        self.pdb_output: Optional[str] = kwargs.get("pdb_output")
        self.peptide: Optional[str] = kwargs.get("peptide")
        self.peptide_structure_prediction: Optional[str] = kwargs.get(
            "peptide_structure_prediction"
        )
        self.peptide_rmsd_mode: str = kwargs.get("peptide_rmsd_mode", "auto")
        self.peptide_min_aligned_length: int = kwargs.get(
            "peptide_min_aligned_length", 5
        )
        self.peptide_min_aligned_fraction: float = kwargs.get(
            "peptide_min_aligned_fraction", 0.5
        )
        self.peptide_min_identity: float = kwargs.get("peptide_min_identity", 0.5)
        self.peptide_max_gap_fraction: float = kwargs.get(
            "peptide_max_gap_fraction", 0.4
        )
        self.peptide_min_contiguous_block: int = kwargs.get(
            "peptide_min_contiguous_block", 3
        )
        self.protein_category: Optional[str] = kwargs.get("protein_category")
        self.protein_flexibility: Optional[str] = kwargs.get("protein_flexibility")
        self.protein_plddt: Optional[str] = kwargs.get("protein_plddt")
        self.protein_restraints: Optional[str] = kwargs.get("protein_restraints")
        self.protein_restraints_retain: Optional[bool] = kwargs.get(
            "protein_restraints_retain"
        )
        self.no_protein_restraints: Optional[bool] = kwargs.get("no_protein_restraints")
        self.no_progress_bar: Optional[bool] = kwargs.get("no_progress_bar")
        self.random_seed: Optional[int] = kwargs.get("random_seed")
        self.receptor_ss: Optional[str] = kwargs.get("receptor_ss")
        self.reference_pdb: Optional[str] = kwargs.get("reference_pdb")
        self.remote: Optional[bool] = kwargs.get("log")
        self.renumber: Optional[bool] = kwargs.get("renumber_residues_to_original")
        self.replicas: Optional[int] = kwargs.get("replicas")
        self.replicas_dtemp: Optional[float] = kwargs.get("replicas_dtemp")
        self.restraints_output: Optional[bool] = kwargs.get("restraints_output")
        self.save_cabs_files: Optional[str] = kwargs.get("save_cabs_files")
        self.save_config: Optional[bool] = kwargs.get("save_config")
        self.sc_rest_add: Optional[str] = kwargs.get("sc_rest_add")
        self.sc_rest_file: Optional[str] = kwargs.get("sc_rest_file")
        self.sc_rest_weight: Optional[float] = kwargs.get("sc_rest_weight")
        self.weighted_fit: Optional[str] = kwargs.get("weighted_fit")
        self.disable_side_chain_centers = kwargs.get("disable_side_chain_centers")
        self.separation: Optional[str] = kwargs.get("separation")
        self.ss_output: Optional[bool] = kwargs.get("ss_output")
        self.temperature: Optional[float] = kwargs.get("temperature")
        self.verbose: Optional[int] = kwargs.get("verbose")
        self.work_dir: Optional[str] = kwargs.get("work_dir")
        self.write_sc_start_pdbs = kwargs.get("write_sc_start_pdbs")

        if self.cg2all_representation == "calpha" and self.disable_side_chain_centers is None:
            self.disable_side_chain_centers = True
            logger.info(_name, "Side-chain centers disabled by default for 'calpha' representation.")
        elif self.disable_side_chain_centers is None:
            self.disable_side_chain_centers = False

        if self.write_sc_start_pdbs is None:
            self.write_sc_start_pdbs = False

        # Job attributes collected.
        self.config: Dict[str, Any] = kwargs
        self.initial_complex: Optional[Any] = None
        self.restraints: Optional[Any] = None
        self.cabsrun: Optional[Any] = None
        self.trajectory: Optional[Trajectory] = None
        self.filtered_trajectory: Optional[Trajectory] = None
        self.filtered_ndx: Optional[npt.NDArray[np.int_]] = None
        self.medoids: Optional[Trajectory] = None
        self.clusters_dict: Optional[Dict[str, Any]] = None
        self.clusters: Optional[npt.NDArray[np.float64]] = None
        self.rmslst: Dict[str, Any] = {}
        self.results: Optional[Dict[str, Any]] = None
        self.reference: Optional[Any] = None
        self.file_TRAF: Optional[str] = None
        self.file_SEQ: Optional[str] = None

        # seeding RNG
        random.seed(self.random_seed)

        # Work_dir processing: making sure work_dir is abspath
        self.work_dir = os.path.abspath(self.work_dir)

        self.generate_visualizations = any(
            (
                self.generate_pymol_visualizations,
                self.generate_chimera_visualizations,
                self.generate_notebook,
            )
        )
        if self.generate_visualizations:
            if not self.aa_rebuild:
                self.aa_rebuild = "M"
            elif "M" not in self.aa_rebuild and "A" not in self.aa_rebuild and "N" not in self.aa_rebuild:
                self.aa_rebuild += "M"
            self.pdb_output = "A"
            self.pdb_bfac_output = "A"
            self.restraints_output = True
            self.contact_maps = True
            self.ss_output = True

        if "N" in self.aa_rebuild:
            self.aa_rebuild = ""

        if "T" in self.aa_rebuild or "A" in self.aa_rebuild:
            if not self.pdb_output or "N" in self.pdb_output:
                self.pdb_output = "R"
            elif "A" not in self.pdb_output and "R" not in self.pdb_output:
                self.pdb_output += "R"

        try:
            logger.setup(
                log_level=self.verbose,
                remote=self.remote,
                work_dir=self.work_dir,
                save_dssp=self.dssp_output,
                save_ss=self.ss_output,
                save_restraints=self.restraints_output,
                progress_bar=not self.no_progress_bar,
            )
            os.makedirs(self.work_dir)
        except OSError:
            if os.path.isdir(self.work_dir):
                logger.warning(
                    _name,
                    f"{self.work_dir} already exists. Output data will be overwritten.",
                )
            else:
                logger.exit_program(
                    _name,
                    f"{self.work_dir} already exists and is not a directory. Choose different name.",
                )

        if not self.cg2all_env_prefix:
            try:
                # Use the new dedicated utility function
                discovered_path = get_cg2all_env_prefix()
                if discovered_path:
                    self.cg2all_env_prefix = discovered_path
                    logger.info(_name, f"Discovered cg2all environment path: {self.cg2all_env_prefix}")
            except Exception:
                logger.debug(_name, "Skipping automatic cg2all environment discovery.")
                pass
        if self.fortran_command:
            cabs.CabsRun.FORTRAN_COMMAND = self.fortran_command

        if self.disable_centro:
            cabs.CabsRun.FORCE_FIELD = tuple(
                0.0 if i == 3 else cabs.CabsRun.FORCE_FIELD[i]
                for i, _ in enumerate(cabs.CabsRun.FORCE_FIELD)
            )

        if self.nsp3_model_path:
            protein.Protein.NSP3_MODEL_PATH = self.nsp3_model_path

        if "T" in self.aa_rebuild or "A" in self.aa_rebuild:
            logger.warning(
                _name,
                "All-atom reconstruction of trajectory replicas may take significant time "
                "(hours to days), depending on the system size and trajectory length.",
            )

        self.file_TRAF = self.file_SEQ = None
        if self.load_cabs_files:
            try:
                self.load_cabs_results()
                self.file_TRAF = os.path.join(self.work_dir, "TRAF")
                self.file_SEQ = os.path.join(self.work_dir, "SEQ")
            except (OSError, ValueError, TypeError) as e:
                logger.exit_program(
                    module_name=_name,
                    msg=f"Could not load CABS files from {self.load_cabs_files}. An error occurred: {e}",
                    exc=e,
                )

        # self.peptide + self.add_peptide -> self.ligand
        self.peptides = []
        if self.peptide:
            self.peptides.extend([[p, "random", "random"] for p in self.peptide])
        if self.add_peptide:
            self.peptides.extend([p for p in self.add_peptide if p])

        valid_letters = set("RFCMSAN")

        try:
            if self.pdb_output:
                if not all(letter in valid_letters for letter in self.pdb_output):
                    raise ValueError("Contains letters outside of 'RFCMSAN'.")

                # Process 'A' or 'N' in pdb_output
                if "A" in self.pdb_output:
                    self.pdb_output = "RFCMS"
                elif "N" in self.pdb_output:
                    self.pdb_output = ""

        except ValueError as e:
            logger.exit_program(
                module_name=_name,
                msg="Invalid pdb_output. An error occurred: %s" % e,
                exc=e,
            )

        valid_bfac_letters = set("ABCPRSN")

        try:
            if self.pdb_bfac_output:
                if not all(letter in valid_bfac_letters for letter in self.pdb_bfac_output):
                    raise ValueError("Contains letters outside of 'ABCPRSN'.")

                if "A" in self.pdb_bfac_output:
                    self.pdb_bfac_output = "BCPRS"
                elif "N" in self.pdb_bfac_output:
                    self.pdb_bfac_output = ""

            if self.pdb_bfac_output or "S" in self.pdb_output:
                self.save_initial_pdb = True
            else:
                self.save_initial_pdb = False

        except ValueError as e:
            logger.exit_program(
                module_name=_name,
                msg="Invalid pdb_bfac_output. An error occurred: %s" % e,
                exc=e,
            )

        valid_csv_letters = set("ABCPSN")

        try:
            if not all(letter in valid_csv_letters for letter in self.csv_output):
                raise ValueError("Contains letters outside of 'ABCPSN'.")

            if "A" in self.csv_output:
                self.csv_output = "BCPS"
            elif "N" in self.csv_output:
                self.csv_output = ""

        except ValueError as e:
            logger.exit_program(
                module_name=_name,
                msg="Invalid csv_output. An error occurred: %s" % e,
                exc=e,
            )

        if self.contact_map_colors:
            self.colors = self.contact_map_colors
        else:
            self.colors = DEFAULT_COLORS

        # Flag to check if dynamic weights should be used
        self.gauss = self.weighted_fit == "gauss"
        if self.gauss and self.gauss_iterations:
            utils.GAUSS_MAX_ITER = self.gauss_iterations

        allowed_modes = ["rigid", "plddt", "manual", "flexible", "ss1", "none", "unleashed", "no-protein-restraints"]
        if not self.no_protein_restraints:
            args = self.protein_restraints
            if isinstance(args, str):
                args = [args]
                
            mode = args[0]
            mode_lower = mode.lower()
            if mode_lower == "all":
                logger.warning(
                    _name,
                    "Protein restraints mode 'all' is legacy. Mapping to the new equivalent 'rigid'.",
                )
                mode = "rigid"
            elif mode_lower == "ss2":
                logger.warning(
                    _name,
                    "Protein restraints mode 'ss2' is legacy. Mapping to the new equivalent 'flexible'.",
                )
                mode = "flexible"

            if len(args) == 4:
                _, gap, min_d, max_d = args
                gap = int(gap)
                min_d = float(min_d)
                max_d = float(max_d)
            elif len(args) == 1:
                # Provide application specific defaults based on child class name
                if self.__class__.__name__ == "DockTask":
                    gap, min_d, max_d = 5, 5.0, 15.0
                else:  # FlexTask
                    if mode.lower() in ["flexible", "ss1"]:
                        gap, min_d, max_d = 3, 3.8, 11.5
                    elif mode.lower() == "rigid":
                        gap, min_d, max_d = 5, 5.0, 15.0
                    else: 
                        gap, min_d, max_d = 5, 5.0, 15.0
            else:
                logger.exit_program(_name, "Invalid number of arguments for --protein-restraints. Expected 1 or 4 arguments.")

            self.protein_restraints = (mode, gap, min_d, max_d)
            if mode in ["manual", "plddt"] and not (
                self.protein_plddt or self.protein_category
            ):
                logger.warning(
                    _name,
                    "No information about pLDDT or flexibility categories provided. "
                    "Changing protein restraints  mode to 'rigid'. "
                    "If you want to use restraints based on pLDDT or flexibility categories, "
                    "please provide the necessary data.",
                )
                self.protein_restraints = ("rigid", gap, min_d, max_d)
            elif (
                mode.lower() == "none"
                or mode.lower() == "unleashed"
                or mode.lower() == "no-protein-restraints"
            ):
                self.no_protein_restraints = True
            elif mode not in allowed_modes:
                # Use the predefined default for this task type
                parser = opt_parser.dock_parser if self.__class__.__name__ == "DockTask" else opt_parser.flex_parser
                default_mode = parser.get_default("protein_restraints")[0]
                
                logger.warning(_name, "Unknown mode: '%s'. Changing to '%s'." % (mode, default_mode))
                self.protein_restraints = (default_mode, gap, min_d, max_d)

        if self.no_protein_restraints:
            self.category_mode = "unleashed"
        else:
            if self.protein_restraints[0].lower() in ["flexible", "ss1"]:
                self.category_mode = "flexible"
            else:
                self.category_mode = "rigid"

        # pairwise potential modification
        if self.pairmod:
            pairmod = {}
            try:
                with open(self.pairmod) as f:
                    for line in f:
                        r, w, s = line.replace("PEP ", "PEP1 ").split()
                        pairmod[r] = (float(w), float(s))
                    logger.info(
                        module_name=_name,
                        msg="Using pairwise potential modification from "
                        + self.pairmod,
                    )
                    self.pairmod = pairmod
            except (OSError, ValueError):
                logger.warning(_name, "Error while reading file: " + self.pairmod)
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
        num_available = self.trajectory.coordinates.shape[0] * self.trajectory.coordinates.shape[1]
        if self.clustering_medoids > num_available:
            logger.warning(
                _name,
                f"The number of medoids {self.clustering_medoids} exceeds the number of structures to be clustered {num_available}. Reducing medoids count to {num_available}.",
            )
            self.clustering_medoids = num_available

        self.score_results(
            n_filtered=self.filtering_count,
            number_of_medoids=self.clustering_medoids,
            number_of_iterations=self.clustering_iterations,
        )
        if self.pdb_output:
            output_folder = self.save_models()
            reconstruct_job_outputs(self, output_folder)
        if self.reference:
            try:
                self.calculate_rmsd()
            except (ValueError, AlignError) as e:
                logger.critical(module_name=_name, msg=str(e))
        self.save_config_file()
        self.load_all_atom_tops()
        self.draw_plots(colors=self.colors)
        if self.pdb_bfac_output:
            self.save_bfac_models()
        if self.csv_output:
            self.save_csv_files()
            
        if self.generate_visualizations:
            self.generate_analysis_visualizations()

        if self.load_cabs_files:
            for _file in CABS_FILES:
                try:
                    os.remove(os.path.join(self.work_dir, _file))
                except OSError:
                    pass
        logger.info(module_name=_name, msg="Simulation completed successfully")

    def generate_analysis_visualizations(self):
        start_pdb_path = os.path.join(self.work_dir, "output_pdbs", "start.pdb")
        model_pdbs = [
            os.path.join(self.work_dir, "output_pdbs", f"model_{i}.pdb")
            for i in range(self.clustering_medoids)
        ]
        restraints_file = os.path.join(self.work_dir, "output_data", "restraints.txt")

        if self.generate_pymol_visualizations:
            logger.info(module_name=_name, msg="Generating PyMOL visualization scripts...")
            try:
                pymol_utils.generate_pymol_scripts(
                    work_dir=self.work_dir,
                    start_pdb_path=start_pdb_path,
                    models_pdbs=model_pdbs,
                    restraints_file=restraints_file
                )
            except Exception as e:
                logger.warning(_name, f"Failed to generate PyMOL scripts: {e}")

        if self.generate_chimera_visualizations:
            logger.info(module_name=_name, msg="Generating ChimeraX visualization scripts...")
            try:
                chimerax_utils.generate_chimerax_scripts(
                    work_dir=self.work_dir,
                    start_pdb_path=start_pdb_path,
                    models_pdbs=model_pdbs,
                    restraints_file=restraints_file,
                )
            except Exception as e:
                logger.warning(_name, f"Failed to generate ChimeraX scripts: {e}")

        if self.generate_notebook:
            logger.info(module_name=_name, msg="Generating analysis notebook...")
            try:
                notebook_utils.generate_notebook(self.work_dir, export_html=True)
            except Exception as e:
                logger.warning(_name, f"Failed to generate analysis notebook: {e}")

    def save_cabs_res(self):
        with NamedTemporaryFile(
            prefix=strftime(cabs.CabsRun.CABS_DIR_FMT),
            dir=self.work_dir,
            suffix=".cbs",
            delete=False,
        ) as temp_file:
            tar_dir = temp_file.name

        with tarfile.open(tar_dir, "w:gz") as tar:
            logger.log_file(_name, "Saving CABS simulation files to: %s" % tar_dir)
            for file_name in CABS_FILES:
                try:
                    tar.add(
                        os.path.join(self.cabsrun.cfg["cwd"], file_name),
                        arcname=file_name,
                    )
                except OSError:
                    pass
                    # print("Could not add %s to the tar file" % file_name)

    def load_cabs_results(self):
        if not os.path.exists(self.load_cabs_files):
            logger.exit_program(
                module_name=_name,
                msg="Provided CABS files path does not exist (%s)"
                % self.load_cabs_files,
                traceback=False,
            )
        try:
            files = glob.glob(os.path.join(self.load_cabs_files, "*.cbs"))
            if len(files) > 1:
                logger.critical(
                    module_name=_name,
                    msg="More than one .cbs file in provided directory %s "
                    % " \n".join(files),
                )
                logger.exit_program(
                    module_name=_name,
                    msg="Please re-run with --load-cabs-files <filename> "
                    "or remove the files you do not need. Quiting.",
                    traceback=False,
                )
            elif len(files) == 1:
                logger.info(
                    module_name=_name, msg="Loading CABS files from %s" % files[0]
                )
                with tarfile.open(files[0], "r:gz") as f:
                    f.extractall(os.path.join(self.work_dir))

            else:
                raise OSError

        except OSError:
            files_loc = self.load_cabs_files
            try:
                with tarfile.open(files_loc, "r:gz") as f:
                    logger.info(
                        module_name=_name, msg="Loading CABS files from %s" % files_loc
                    )
                    f.extractall(os.path.join(self.work_dir))
            except OSError:
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
            self.align_options,
        )
        mtx_p = mtx_p.to_numpy()
        mtx_q = mtx_q.to_numpy()
        dummy_rmsd, rot, t_com, q_com = utils.dynamic_kabsch(mtx_p, mtx_q)
        self.reference[0].from_numpy(
            np.dot(self.reference[0].to_numpy() - q_com, rot) + t_com
        )

    def load_all_atom_tops(self):
        if self.aa_rebuild:
            pth = os.path.join(self.work_dir, "output_pdbs", "model_%i.pdb")
            pdb_files = [pth % i for i in range(self.clustering_medoids)]
            coords_list = [Pdb(f, create_from_aa=True).atoms.to_numpy() for f in pdb_files]
            try:
                med_traj = np.array(coords_list)
            except ValueError:
                logger.warning(_name, "Inconsistent atom counts in all-atom models. Clustering/RMSD might be affected.")
                min_atoms = min(c.shape[0] for c in coords_list)
                med_traj = np.array([c[:min_atoms] for c in coords_list])
            
            med_traj = np.expand_dims(med_traj, axis=1)
            mod0 = Pdb(
                os.path.join(self.work_dir, "output_pdbs", "model_0.pdb"),
                create_from_aa=True,
            )
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

        cmapdir = os.path.join(self.work_dir, "contact_maps")
        try:
            os.mkdir(cmapdir)
        except OSError:
            pass

        return sc_traj_full, sc_med, cmapdir

    def prepare_restraints(self, output=""):
        # generate protein restraints
        if self.no_protein_restraints:
            protein_restraints = Restraints(None)
        else:
            protein_restraints = Restraints(
                self.initial_complex.protein.generate_restraints(
                    *self.protein_restraints
                )
            )

        # reduce number of restraints
        if self.protein_restraints_retain:
            if self.protein_restraints_retain < 100:
                protein_restraints.retain_percentage(self.protein_restraints_retain)
                logger.debug(
                    module_name=_name,
                    msg=f"Retaining {self.protein_restraints_retain}% of restraints.",
                )

        # additional restraints
        add_restraints = Restraints("")

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
            mapped_cyclization = [self.initial_complex.map_user_chain(c) for c in self.cyclization]
            protein_restraints += Restraints(
                self.initial_complex.generate_backbone_restraints(
                    mapped_cyclization
                )
            )

        if self.disulfide_bonds:
            mapped_disulfides = []
            for b1, b2 in self.disulfide_bonds:
                mapped_disulfides.append(
                    (
                        self.initial_complex.map_user_residue(b1),
                        self.initial_complex.map_user_residue(b2),
                    )
                )
            protein_restraints += Restraints(
                self.initial_complex.generate_disulfide_restraints(
                    mapped_disulfides
                ),
                sg=True,
            )

        protein_restraints += add_restraints.update_id(self.initial_complex.new_ids)

        if output and logger.output_restraints():
            restraints_for_output = deepcopy(protein_restraints).update_id(
                self.initial_complex.old_ids
            )
            output_restraints = os.path.join(output, "output_data", "restraints.txt")
            odir = os.path.dirname(output_restraints)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            logger.to_file(
                filename=output_restraints,
                content=str(restraints_for_output),
                msg="Saving restraints output to %s" % output_restraints,
            )

        return protein_restraints

    def save_config_file(self):
        if self.save_config:
            with open(os.path.join(self.work_dir, "config.ini"), "w") as configfile:
                configfile.write(CONFIG_HEADER)   
                for k in sorted(self.config):
                    value = self.config[k]
                    name = re.sub("_", "-", str(k))
                    option = opt_parser.option_formatter(name, value)
                    try:
                        configfile.write(option)
                    except Exception as e:
                        logger.warning(
                            _name,
                            "Failed to save %s option to config file. Reason: %s."
                            % (name, e),
                        )

    def setup_cabs_run(self):
        logger.info(module_name="CABS", msg="Setting up CABS simulation.")

        # --- Memory Warning Logic ---
        n_mols = len(self.initial_complex.chain_list)
        if n_mols > 10:
            logger.warning(
                "CABS",
                f"Large system detected ({n_mols} chains). This simulation may require "
                "significant RAM per replica. Ensure enough memory is available."
            )

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
            pairmod=self.pairmod,
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
                module_name=_name, msg=f"Loading trajectories from: {ftraf}, {fseq}"
            )
            self.trajectory = Trajectory.read_trajectory(ftraf, fseq)
        else:
            logger.debug(
                module_name=_name, msg="Loading trajectories from the CABS run"
            )
            self.trajectory = self.cabsrun.get_trajectory()
        self.trajectory.weights = self.initial_complex.weights
        chs = list(self.initial_complex.list_chains().keys())
        tchs = [c for c in chs if c in self.trajectory.template.list_chains()]
        self.trajectory.tmp_target_chs = tchs
        ref_stc = self.initial_complex.protein
        if not ref_stc:
            ref_stc = self.initial_complex.select("model 1")

        ic_stc, tt_stc, dummy_aln = self.trajectory.align_to(
            ref_stc, chs, tchs, align_mth="trivial"
        )
        self.trajectory.superimpose_to(ic_stc, tt_stc)

        # Update IDs after alignment to avoid issues with multi-character chain IDs (like PEP1) in selection strings
        self.trajectory.template.update_ids(
            self.initial_complex.old_ids, pedantic=False
        )
        self.initial_complex.update_ids(self.initial_complex.old_ids, pedantic=False)
        logger.info(module_name=_name, msg="Trajectories loaded successfully")
        return self.trajectory

    @abstractmethod
    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        pass

    def save_models(self):
        output_folder = os.path.join(self.work_dir, "output_pdbs")
        
        # Safe mapping of peptide aliases (PEP1, PEP2...) to actual chain IDs
        def map_resid(resid):
            if self.initial_complex and resid in self.initial_complex.new_ids:
                return self.initial_complex.new_ids[resid]
            return resid

        mapped_cyclization = [map_resid(f"1:{c}") for c in self.cyclization] if self.cyclization else []
        mapped_cyclization = [c.split(":")[1] if ":" in c else c for c in mapped_cyclization]

        mapped_disulfides = [
            (map_resid(b[0]), map_resid(b[1]))
            for b in self.disulfide_bonds
        ] if self.disulfide_bonds else []

        logger.log_file(
            module_name=_name, msg="Saving pdb files to " + str(output_folder)
        )
        try:
            os.mkdir(output_folder)
        except OSError:
            logger.warning(
                _name,
                "Possibly overwriting previous pdb files. Use --work-dir <DIR> to avoid that.",
            )
        # Saving the trajectory to PDBs:
        sc_flag = not self.disable_side_chain_centers
        if "R" in self.pdb_output:
            logger.log_file(module_name=_name, msg="Saving replicas...")
            self.trajectory.to_pdb(mode="replicas", to_dir=output_folder, sc=sc_flag)
            if ("T" in self.aa_rebuild or "A" in self.aa_rebuild) and self.aa_method == "modeller":
                self.save_all_atom_replicas(output_folder)
        # Saving top1000 models to PDB:
        if "F" in self.pdb_output:
            logger.log_file(module_name=_name, msg="Saving filtered models...")
            self.filtered_trajectory.to_pdb(
                mode="replicas", to_dir=output_folder, name="top1000", sc=sc_flag
            )
        # Saving clusters in CA representation
        if "C" in self.pdb_output:
            logger.log_file(module_name=_name, msg="Saving clusters...")
            for i, cluster in enumerate(self.clusters):
                cluster.to_pdb(
                    mode="replicas", to_dir=output_folder, name=f"cluster_{i}", sc=sc_flag
                )
        if "S" in self.pdb_output:
            logger.log_file(module_name=_name, msg="Saving starting structure...")
            if self.write_sc_start_pdbs:
                temp_complex = deepcopy(self.initial_complex)
                temp_complex.add_side_chain_centers()
                temp_complex.save_to_pdb(os.path.join(output_folder, "start.pdb"))
            else:
                self.initial_complex.save_to_pdb(os.path.join(output_folder, "start.pdb"))
        
        # Saving final models:
        if "M" in self.pdb_output:
            save_to_ca = True
            odir = os.path.join(self.work_dir, "output_data")
            if not os.path.isdir(odir):
                os.makedirs(odir)
            if "M" in self.aa_rebuild or "A" in self.aa_rebuild:
                if self.aa_method in ALLOWED_AA_METHODS:
                    logger.log_file(
                        module_name=_name,
                        msg=f"Saving final models (in AA representation) using {self.aa_method}.",
                    )
                    if self.aa_method == "modeller":
                        try:
                            from CABS.reconstruction.ca2all import ca2all
                        except ImportError:
                            logger.warning(
                                _name, msg="Modeller not found. Skipping AA rebuild."
                            )
                        else:
                            logger.log_file(
                                module_name=_name,
                                msg=f"Running Modeller to rebuild models (iterations: {self.modeller_iterations}).",
                            )
                            pdb_medoids = self.medoids.to_pdb(sc=False)
                            for i, fname in enumerate(pdb_medoids):
                                ca2all(
                                    fname,
                                    output=os.path.join(
                                        output_folder, f"model_{i}.pdb"
                                    ),
                                    iterations=self.modeller_iterations,
                                    out_mdl=os.path.join(
                                        self.work_dir,
                                        "output_data",
                                        f"modeller_output_{i}.txt",
                                    ),
                                    work_dir=self.work_dir,
                                    cyclization=mapped_cyclization,
                                    disulfide_bonds=mapped_disulfides,
                                )
                            save_to_ca = False
                    elif self.aa_method == "cg2all":
                        logger.log_file(
                            module_name=_name,
                            msg=f"Running cg2all to rebuild models (env_prefix: {self.cg2all_env_prefix}).",
                        )
                        attempt_cyclization = False
                        if self.cyclization or self.disulfide_bonds:
                            try:
                                from CABS.reconstruction.ca2all import ca2all
                            except ImportError:
                                logger.warning(
                                    _name,
                                    msg="Modeller not found. Skipping backbone and/or disulfide cyclization.",
                                )
                            else:
                                attempt_cyclization = True
                        pdb_medoids = self.medoids.to_pdb(sc=False)
                        original_chains = "".join(self.medoids.template.list_chains())
                        for i, fname in enumerate(pdb_medoids):
                            logger.info(
                                module_name="CG2ALL",
                                msg=f"Reconstructing model_{i}...",
                            )
                            convert_cg_to_all(
                                fname,
                                work_dir=self.work_dir,
                                iter=i,
                                reference_pdb=self.input_protein,
                                renumber_flag=self.renumber,
                                env_prefix=self.cg2all_env_prefix,
                                cg2all_representation=self.cg2all_representation,
                                minimize_flag=(self.aa_minimize is not False),
                                n_proc=self.aa_rebuild_workers,
                            )
                            if attempt_cyclization:
                                pth_tmp = os.path.join(
                                    self.work_dir, "output_pdbs", f"model_{i}.pdb"
                                )
                                with open(pth_tmp) as f:
                                    ca2all(
                                        f,
                                        output=os.path.join(
                                            output_folder, f"model_{i}.pdb"
                                        ),
                                        iterations=1,
                                        out_mdl=os.path.join(
                                            self.work_dir,
                                            "output_data",
                                            f"modeller_output_{i}.txt",
                                        ),
                                        work_dir=self.work_dir,
                                        cyclization=mapped_cyclization,
                                        disulfide_bonds=mapped_disulfides,
                                        only_cyclization=True,
                                    )
                        save_to_ca = False
                else:
                    logger.warning(
                        module_name=_name,
                        msg="Unknown AA method: %s. Skipping AA rebuild."
                        % self.aa_method,
                    )

            if save_to_ca:
                logger.log_file(
                    module_name=_name, msg="Saving final models (in CA representation)."
                )
                self.medoids.to_pdb(mode="models", to_dir=output_folder, name="model", sc=sc_flag)
                if self.renumber and self.input_protein:
                    start_all_path = os.path.join(self.work_dir, "output_pdbs", "start_all.pdb")
                    if not os.path.exists(start_all_path):
                        start_all_path = os.path.join(self.work_dir, "output_pdbs", "start.pdb")
                    if os.path.exists(start_all_path):
                        from CABS.reconstruction.cg2all import sync_residues_with_template
                        from pathlib import Path
                        for i in range(self.clustering_medoids):
                            model_path = os.path.join(output_folder, f"model_{i}.pdb")
                            if os.path.exists(model_path):
                                try:
                                    sync_residues_with_template(
                                        input_pdb_path=Path(self.input_protein.split(":")[0]),
                                        topology_pdb_path=Path(start_all_path),
                                        output_pdb_path=Path(model_path)
                                    )
                                except Exception as e:
                                    logger.warning(
                                        _name,
                                        msg=f"Could not renumber CA model {i}: {e}"
                                    )

            if self.json_output:
                json_file = os.path.join(self.work_dir, "output_data", "medoid.json")
                odir = os.path.dirname(json_file)
                if not os.path.isdir(odir):
                    os.makedirs(odir)
                logger.log_file(module_name=_name, msg="Saving JSON output.")
                medoids_ca_atoms_list = self.medoids.to_atoms_list(sc=sc_flag)
                medoids_ca_atoms_list[0].save_to_json(json_file)
        return output_folder

        def map_resid(resid):
            if self.initial_complex and resid in self.initial_complex.new_ids:
                return self.initial_complex.new_ids[resid]
            return resid

        mapped_cyclization = [map_resid(f"1:{c}") for c in self.cyclization] if self.cyclization else []
        mapped_cyclization = [c.split(":")[1] if ":" in c else c for c in mapped_cyclization]

        mapped_disulfides = [
            (map_resid(b[0]), map_resid(b[1]))
            for b in self.disulfide_bonds
        ] if self.disulfide_bonds else []

        if self.aa_method not in ALLOWED_AA_METHODS:
            logger.warning(
                _name,
                f"Unknown AA method: {self.aa_method}. Skipping replica reconstruction.",
            )
            return

        ca2all = None
        if self.aa_method == "modeller":
            try:
                from CABS.reconstruction.ca2all import ca2all
            except ImportError:
                logger.warning(
                    _name,
                    "Modeller not found. Skipping all-atom replica reconstruction.",
                )
                return

        logger.warning(
            _name,
            "Reconstructing all trajectory replica frames to all-atom representation. "
            "This may take hours to days for large systems.",
        )

        shape = self.trajectory.coordinates.shape
        for replica_idx in range(shape[0]):
            for model_idx in range(shape[1]):
                ca_model = self.trajectory.get_model(replica_idx * shape[1] + model_idx)
                ca_path = os.path.join(
                    output_folder, f"replica_{replica_idx}_model_{model_idx}_ca.pdb"
                )
                aa_name = f"replica_{replica_idx}_model_{model_idx}.pdb"
                aa_path = os.path.join(output_folder, aa_name)
                ca_model.save_to_pdb(ca_path)
                try:
                    if self.aa_method == "modeller":
                        ca2all(
                            ca_path,
                            output=aa_path,
                            iterations=self.modeller_iterations,
                            out_mdl=os.path.join(
                                self.work_dir,
                                "output_data",
                                f"modeller_replica_{replica_idx}_model_{model_idx}.txt",
                            ),
                            work_dir=self.work_dir,
                            cyclization=mapped_cyclization,
                            disulfide_bonds=mapped_disulfides,
                        )
                    elif self.aa_method == "cg2all":
                        attempt_cyclization = False
                        if self.cyclization or self.disulfide_bonds:
                            try:
                                from CABS.reconstruction.ca2all import ca2all
                                attempt_cyclization = True
                            except ImportError:
                                logger.warning(
                                    _name,
                                    msg="Modeller not found. Skipping backbone and/or disulfide cyclization for replica reconstruction.",
                                )
                        
                        logger.info(
                            module_name="CG2ALL",
                            msg=f"Reconstructing replica frame: {aa_name}...",
                        )
                        convert_cg_to_all(
                            ca_path,
                            work_dir=self.work_dir,
                            iter=model_idx,
                            reference_pdb=self.input_protein,
                            renumber_flag=self.renumber,
                            env_prefix=self.cg2all_env_prefix,
                            output_filename=aa_name,
                            cg2all_representation=self.cg2all_representation,
                            minimize_flag=(self.aa_minimize is True),
                            n_proc=self.aa_rebuild_workers,
                        )

                        if attempt_cyclization:
                            with open(aa_path) as f:
                                ca2all(
                                    f,
                                    output=aa_path,
                                    iterations=1,
                                    out_mdl=os.path.join(
                                        self.work_dir,
                                        "output_data",
                                        f"modeller_replica_{replica_idx}_model_{model_idx}.txt",
                                    ),
                                    work_dir=self.work_dir,
                                    cyclization=mapped_cyclization,
                                    disulfide_bonds=mapped_disulfides,
                                    only_cyclization=True,
                                )
                except Exception as e:
                    logger.warning(
                        _name,
                        f"Failed to reconstruct replica {replica_idx}, model {model_idx}: {e}",
                    )

    def save_bfac_models(self):
        pdb_output = os.path.join(self.work_dir, "output_pdbs")
        if not os.path.isdir(pdb_output):
            os.makedirs(pdb_output)
        logger.log_file(
            module_name=_name,
            msg="Saving starting structures with different beta factors to "
            + str(pdb_output),
        )

        try:
            initial_pdb = os.path.join(pdb_output, "start_all.pdb")
            if not os.path.exists(initial_pdb):
                raise FileNotFoundError
            initial_pdb_file = Pdb(initial_pdb, create_from_aa=True).atoms
            initial_pdb_file.update_occ(self.initial_complex.get_occ())
        except FileNotFoundError:
            logger.warning(
                _name, "No start_all.pdb file found. Skipping beta-factor calculations."
            )
            return

        if "B" in self.pdb_bfac_output:
            logger.log_file(
                module_name=_name, msg="Saving starting structure with beta-factors..."
            )
            bfac_update_dict = self.initial_complex.get_bfac()
            initial_pdb_file.update_bfac(bfac_update_dict)
            initial_pdb_file.save_to_pdb(os.path.join(pdb_output, "start_bfac.pdb"))

        if "C" in self.pdb_bfac_output:
            logger.log(
                module_name=_name,
                msg="Saving starting structure with flexibility categories...",
            )
            category_update_dict = self.initial_complex.get_category()
            initial_pdb_file.update_bfac(category_update_dict)
            initial_pdb_file.save_to_pdb(os.path.join(pdb_output, "start_category.pdb"))

        if "P" in self.pdb_bfac_output:
            logger.log(
                module_name=_name, msg="Saving starting structure with pLDDT values..."
            )
            plddt_update_dict = self.initial_complex.get_plddt()
            for key in plddt_update_dict:
                plddt_update_dict[key] = plddt_update_dict[key] * 100
            initial_pdb_file.update_bfac(plddt_update_dict)
            initial_pdb_file.save_to_pdb(os.path.join(pdb_output, "start_plddt.pdb"))

        if "R" in self.pdb_bfac_output:
            logger.log(
                module_name=_name, msg="Saving starting structure with RMSF values..."
            )
            chs = self.initial_complex.protein_chains or self.initial_complex.peptide_chains
            rmsfs = self.trajectory.rmsf(chs)
            rmsf_update_dict = {}
            atom_index = 0
            for atom in self.trajectory.template.atoms:
                if atom.chid in chs:
                    rmsf_update_dict[atom.resid_id()] = rmsfs[atom_index]
                    atom_index += 1
            initial_pdb_file.update_bfac(rmsf_update_dict)
            initial_pdb_file.save_to_pdb(os.path.join(pdb_output, "start_rmsf.pdb"))

        if "S" in self.pdb_bfac_output:
            logger.log(
                module_name=_name,
                msg="Saving starting structure with secondary structure...",
            )
            ss_update_dict = self.initial_complex.get_occ()
            initial_pdb_file.update_bfac(ss_update_dict)
            initial_pdb_file.save_to_pdb(os.path.join(pdb_output, "start_secstr.pdb"))

    def save_csv_files(self):
        csv_output = os.path.join(self.work_dir, "output_data")
        if not os.path.isdir(csv_output):
            os.makedirs(csv_output)

        logger.log_file(module_name=_name, msg="Saving csv files to " + str(csv_output))

        if "B" in self.csv_output:
            logger.log_file(
                module_name=_name, msg="Saving csv file with beta-factors..."
            )
            bfac_dict = self.initial_complex.get_bfac()
            drop_csv_file(
                os.path.join(csv_output, "bfactor"),
                [list(bfac_dict.keys()), list(bfac_dict.values())],
                fmts=["%s", "%s"],
            )

        if "C" in self.csv_output:
            logger.log(
                module_name=_name, msg="Saving csv file with flexibility categories..."
            )
            category_dict = self.initial_complex.get_category()
            drop_csv_file(
                os.path.join(csv_output, "category"),
                [list(category_dict.keys()), list(category_dict.values())],
                fmts=["%s", "%s"],
            )

        if "P" in self.csv_output:
            logger.log(module_name=_name, msg="Saving csv file with pLDDT values...")
            plddt_dict = self.initial_complex.get_plddt()
            for key in plddt_dict:
                plddt_dict[key] = plddt_dict[key] * 100
            drop_csv_file(
                os.path.join(csv_output, "plddt"),
                [list(plddt_dict.keys()), list(plddt_dict.values())],
                fmts=["%s", "%s"],
            )

        if "S" in self.csv_output:
            logger.log(
                module_name=_name, msg="Saving csv file with secondary structure..."
            )
            ss_dict = self.initial_complex.get_occ()
            resname_dict = self.initial_complex.get_resname()
            drop_csv_file(
                os.path.join(csv_output, "secstr"),
                [
                    list(ss_dict.keys()),
                    list(resname_dict.values()),
                    list(ss_dict.values()),
                ],
                fmts=["%s", "%s", "%s"],
            )


class DockTask(CABSTask):
    """Class representing single CABS job."""

    def setup_job(self):
        if not self.peptides and not self.load_cabs_files:
            raise ValueError("No peptide given")
        self.initial_complex = ProteinComplex(
            protein=self.input_protein,
            flexibility=self.protein_flexibility,
            exclude=self.exclude,
            weights=self.weighted_fit,
            plddt=self.protein_plddt,
            category=self.protein_category,
            mode=self.category_mode,
            peptides=self.peptides,
            replicas=self.replicas,
            separation=self.separation,
            insertion_attempts=self.insertion_attempts,
            insertion_clash=self.insertion_clash,
            work_dir=self.work_dir,
            receptor_ss=self.receptor_ss,
            pdb_cache=self.pdb_cache,
            save_initial_pdb=self.save_initial_pdb,
            json_output=self.json_output,
            predict_peptide_structure=self.peptide_structure_prediction,
            cg2all_env_prefix=self.cg2all_env_prefix,
            sc=self.write_sc_start_pdbs,
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

    @staticmethod
    def _longest_common_subsequence_alignment(ref_atoms, model_atoms):
        ref_seq = [atom.resname for atom in ref_atoms]
        model_seq = [atom.resname for atom in model_atoms]
        dp = [[0] * (len(model_seq) + 1) for _ in range(len(ref_seq) + 1)]

        for i in range(len(ref_seq) - 1, -1, -1):
            for j in range(len(model_seq) - 1, -1, -1):
                if ref_seq[i] == model_seq[j]:
                    dp[i][j] = dp[i + 1][j + 1] + 1
                else:
                    dp[i][j] = max(dp[i + 1][j], dp[i][j + 1])

        i = j = 0
        aligned_pairs = []
        while i < len(ref_seq) and j < len(model_seq):
            if ref_seq[i] == model_seq[j]:
                aligned_pairs.append((ref_atoms[i], model_atoms[j]))
                i += 1
                j += 1
            elif dp[i + 1][j] >= dp[i][j + 1]:
                i += 1
            else:
                j += 1

        if not aligned_pairs:
            raise AlignError("No common peptide subsequence between model and reference.")

        return (
            Atoms([ref_atom for ref_atom, _ in aligned_pairs]),
            Atoms([model_atom for _, model_atom in aligned_pairs]),
            tuple(aligned_pairs),
        )

    @staticmethod
    def _smith_waterman_peptide_alignment(ref_atoms, model_atoms):
        match_score = 2
        mismatch_score = -1
        gap_score = -1
        rows = len(ref_atoms)
        cols = len(model_atoms)
        scores = [[0] * (cols + 1) for _ in range(rows + 1)]
        moves = [[None] * (cols + 1) for _ in range(rows + 1)]
        best_i = best_j = 0
        best_score = 0

        for i in range(1, rows + 1):
            for j in range(1, cols + 1):
                diag = scores[i - 1][j - 1] + (
                    match_score
                    if ref_atoms[i - 1].resname == model_atoms[j - 1].resname
                    else mismatch_score
                )
                up = scores[i - 1][j] + gap_score
                left = scores[i][j - 1] + gap_score
                score, move = max(
                    ((0, None), (diag, "diag"), (up, "up"), (left, "left")),
                    key=lambda item: item[0],
                )
                scores[i][j] = score
                moves[i][j] = move
                if score > best_score:
                    best_score = score
                    best_i = i
                    best_j = j

        if best_score <= 0:
            raise AlignError("No sequential similarity.")

        i, j = best_i, best_j
        aligned_pairs = []
        gap_count = 0
        while i > 0 and j > 0 and scores[i][j] > 0:
            move = moves[i][j]
            if move == "diag":
                aligned_pairs.append((ref_atoms[i - 1], model_atoms[j - 1]))
                i -= 1
                j -= 1
            elif move == "up":
                gap_count += 1
                i -= 1
            elif move == "left":
                gap_count += 1
                j -= 1
            else:
                break

        if not aligned_pairs:
            raise AlignError("No sequential similarity.")

        aligned_pairs.reverse()
        return (
            Atoms([ref_atom for ref_atom, _ in aligned_pairs]),
            Atoms([model_atom for _, model_atom in aligned_pairs]),
            tuple(aligned_pairs),
            gap_count,
        )

    @staticmethod
    def _peptide_alignment_metrics(ref_atoms, model_atoms, pep_aln, gap_count=0):
        aligned_ref_ids = [ref_atom.resid_id() for ref_atom, _ in pep_aln]
        aligned_model_ids = [model_atom.resid_id() for _, model_atom in pep_aln]
        aligned_length = len(pep_aln)
        identical = sum(
            1 for ref_atom, model_atom in pep_aln if ref_atom.resname == model_atom.resname
        )
        mismatch_count = aligned_length - identical
        denominator = aligned_length + gap_count
        identity = identical / float(aligned_length) if aligned_length else 0.0
        longest_block = 0
        current_block = 0
        prev_ref = prev_model = None
        for ref_atom, model_atom in pep_aln:
            if (
                prev_ref is not None
                and ref_atom.resnum == prev_ref.resnum + 1
                and model_atom.resnum == prev_model.resnum + 1
                and ref_atom.chid == prev_ref.chid
                and model_atom.chid == prev_model.chid
            ):
                current_block += 1
            else:
                current_block = 1
            longest_block = max(longest_block, current_block)
            prev_ref = ref_atom
            prev_model = model_atom

        return {
            "reference_length": len(ref_atoms),
            "model_length": len(model_atoms),
            "aligned_length": aligned_length,
            "aligned_fraction_shorter": (
                aligned_length / float(min(len(ref_atoms), len(model_atoms)))
                if min(len(ref_atoms), len(model_atoms))
                else 0.0
            ),
            "identity": identity,
            "gap_count": gap_count,
            "gap_fraction": (gap_count / float(denominator)) if denominator else 0.0,
            "mismatch_count": mismatch_count,
            "longest_contiguous_block": longest_block,
            "excluded_reference_residues": len(ref_atoms) - len(set(aligned_ref_ids)),
            "excluded_model_residues": len(model_atoms) - len(set(aligned_model_ids)),
        }

    def _mode_order(self):
        mode = self.peptide_rmsd_mode
        if mode == "strict":
            return ["strict"]
        if mode == "overlap":
            return ["strict", "overlap"]
        return ["strict", "overlap", "mutational"]

    def _alignment_passes_thresholds(self, mode, metrics):
        if mode == "strict":
            return (
                metrics["aligned_length"] == metrics["reference_length"] == metrics["model_length"]
                and metrics["mismatch_count"] == 0
                and metrics["gap_count"] == 0
            )
        if metrics["aligned_length"] < self.peptide_min_aligned_length:
            return False
        if metrics["aligned_fraction_shorter"] < self.peptide_min_aligned_fraction:
            return False
        if metrics["longest_contiguous_block"] < self.peptide_min_contiguous_block:
            return False
        if mode == "mutational":
            if metrics["identity"] < self.peptide_min_identity:
                return False
            if metrics["gap_fraction"] > self.peptide_max_gap_fraction:
                return False
        return True

    @staticmethod
    def _save_peptide_alignment_csv(filename, pep_aln):
        save_csv(filename, ("reference", "template"), pep_aln)

    @staticmethod
    def _save_peptide_alignment_summary(filename, summary_rows):
        fieldnames = [
            "model_chain",
            "reference_chain",
            "selected_mode",
            "strict_status",
            "strict_reason",
            "overlap_status",
            "overlap_reason",
            "mutational_status",
            "mutational_reason",
            "reference_length",
            "model_length",
            "aligned_length",
            "aligned_fraction_shorter",
            "identity",
            "gap_count",
            "gap_fraction",
            "mismatch_count",
            "longest_contiguous_block",
            "excluded_reference_residues",
            "excluded_model_residues",
        ]
        with open(filename, "w", newline="") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in summary_rows:
                writer.writerow(row)

    def _log_peptide_alignment_warning(
        self,
        ref_pept_chain,
        pept_chain,
        ref_atoms,
        model_atoms,
        pep_aln,
        method_label,
    ):
        aligned_ref_ids = {ref_atom.resid_id() for ref_atom, _ in pep_aln}
        aligned_model_ids = {model_atom.resid_id() for _, model_atom in pep_aln}
        excluded_ref = len(ref_atoms) - len(aligned_ref_ids)
        excluded_model = len(model_atoms) - len(aligned_model_ids)
        if excluded_ref or excluded_model:
            logger.warning(
                module_name=_name,
                msg=(
                    f"Using partial peptide alignment for model chain {pept_chain} "
                    f"against reference chain {ref_pept_chain} via {method_label}: "
                    f"aligned {len(pep_aln)} residue(s), excluded {excluded_ref} "
                    f"reference residue(s) and {excluded_model} model residue(s)."
                ),
            )

    def _select_peptide_alignment(self, ref_pept_chain, pept_chain):
        ref_chain_atoms = self.reference[0].select(
            "name CA and not HETERO and chain %s" % ref_pept_chain
        )
        model_chain_atoms = self.trajectory.template.select(
            "name CA and not HETERO and chain %s" % pept_chain
        )
        results = {
            "strict": {"status": "not-run", "reason": "", "metrics": None, "alignment": None},
            "overlap": {"status": "not-run", "reason": "", "metrics": None, "alignment": None},
            "mutational": {"status": "not-run", "reason": "", "metrics": None, "alignment": None},
        }

        for mode in self._mode_order():
            try:
                if mode == "strict":
                    ref_pep_stc, self_pep_stc, pep_aln = self.trajectory.align_to(
                        self.reference[0],
                        ref_pept_chain,
                        pept_chain,
                        align_mth=self.align,
                        kwargs=self.align_peptide_options,
                    )
                    gap_count = 0
                elif mode == "overlap":
                    ref_pep_stc, self_pep_stc, pep_aln = (
                        self._longest_common_subsequence_alignment(
                            ref_chain_atoms, model_chain_atoms
                        )
                    )
                    gap_count = 0
                else:
                    ref_pep_stc, self_pep_stc, pep_aln, gap_count = (
                        self._smith_waterman_peptide_alignment(
                            ref_chain_atoms, model_chain_atoms
                        )
                    )
            except (ValueError, AlignError) as e:
                results[mode].update(status="failed", reason=str(e))
                continue

            metrics = self._peptide_alignment_metrics(
                ref_chain_atoms, model_chain_atoms, pep_aln, gap_count=gap_count
            )
            results[mode].update(
                metrics=metrics,
                alignment=(ref_pep_stc, self_pep_stc, pep_aln),
            )
            if self._alignment_passes_thresholds(mode, metrics):
                results[mode]["status"] = "selected"
                return {
                    "mode": mode,
                    "metrics": metrics,
                    "alignment": (ref_pep_stc, self_pep_stc, pep_aln),
                    "results": results,
                    "ref_atoms": ref_chain_atoms,
                    "model_atoms": model_chain_atoms,
                }
            results[mode].update(status="rejected", reason="Thresholds not satisfied")

        return {
            "mode": None,
            "metrics": None,
            "alignment": None,
            "results": results,
            "ref_atoms": ref_chain_atoms,
            "model_atoms": model_chain_atoms,
        }

    def calculate_rmsd(self, save=True):
        logger.debug(module_name=_name, msg="RMSD calculations starting...")
        sfname: str = ""
        if save:
            odir = os.path.join(self.work_dir, "output_data")
            try:
                os.mkdir(odir)
            except OSError:
                pass
        all_results = {}
        peptide_alignment_summary = []
        ref_trg_stc, self_trg_stc, trg_aln = self.trajectory.align_to(
            self.reference[0],
            self.reference[1],
            self.trajectory.tmp_target_chs,
            align_mth=self.align,
            kwargs=self.align_options,
        )
        if save:
            sfname = os.path.join(self.work_dir, "output_data", "reference_alignment")
            paln_trg = sfname + "_target.csv"
            save_csv(paln_trg, ("reference", "template"), trg_aln)
        for pept_chain, ref_pept_chain in zip(
            self.initial_complex.peptide_chains, self.reference[2]
        ):
            selection = self._select_peptide_alignment(ref_pept_chain, pept_chain)
            metrics = selection["metrics"]
            ref_chain_atoms = selection["ref_atoms"]
            model_chain_atoms = selection["model_atoms"]
            selected_mode = selection["mode"]
            row = {
                "model_chain": pept_chain,
                "reference_chain": ref_pept_chain,
                "selected_mode": selected_mode or "skipped",
                "strict_status": selection["results"]["strict"]["status"],
                "strict_reason": selection["results"]["strict"]["reason"],
                "overlap_status": selection["results"]["overlap"]["status"],
                "overlap_reason": selection["results"]["overlap"]["reason"],
                "mutational_status": selection["results"]["mutational"]["status"],
                "mutational_reason": selection["results"]["mutational"]["reason"],
                "reference_length": len(ref_chain_atoms),
                "model_length": len(model_chain_atoms),
                "aligned_length": 0,
                "aligned_fraction_shorter": 0.0,
                "identity": 0.0,
                "gap_count": 0,
                "gap_fraction": 0.0,
                "mismatch_count": 0,
                "longest_contiguous_block": 0,
                "excluded_reference_residues": len(ref_chain_atoms),
                "excluded_model_residues": len(model_chain_atoms),
            }
            if metrics:
                row.update(metrics)
            peptide_alignment_summary.append(row)

            if save:
                for mode_name, result in selection["results"].items():
                    if result["alignment"] is not None:
                        _, _, alignment_pairs = result["alignment"]
                        self._save_peptide_alignment_csv(
                            sfname + f"_{pept_chain}_{mode_name}.csv", alignment_pairs
                        )

            if selected_mode is None:
                logger.warning(
                    module_name=_name,
                    msg=(
                        f"Skipping peptide RMSD for model chain {pept_chain} against "
                        f"reference chain {ref_pept_chain}. "
                        f"Reference length {len(ref_chain_atoms)}, model length {len(model_chain_atoms)}."
                    ),
                )
                continue

            ref_pep_stc, self_pep_stc, pep_aln = selection["alignment"]
            if save:
                paln_pep = sfname + "_%s.csv" % pept_chain
                save_csv(paln_pep, ("reference", "template"), pep_aln)
            if selected_mode == "strict":
                logger.info(
                    module_name=_name,
                    msg=(
                        f"Peptide RMSD for model chain {pept_chain} used strict mode: "
                        f"aligned {metrics['aligned_length']}/{metrics['reference_length']} residues."
                    ),
                )
            else:
                logger.warning(
                    module_name=_name,
                    msg=(
                        f"Peptide RMSD for model chain {pept_chain} fell back to "
                        f"{selected_mode} mode: aligned {metrics['aligned_length']} "
                        f"residue(s), identity {metrics['identity']:.3f}, "
                        f"mismatches {metrics['mismatch_count']}, gaps {metrics['gap_count']}."
                    ),
                )
                self._log_peptide_alignment_warning(
                    ref_pept_chain,
                    pept_chain,
                    ref_chain_atoms,
                    model_chain_atoms,
                    pep_aln,
                    method_label=selected_mode,
                )
            if save:
                with open(
                    os.path.join(odir, f"peptide_alignment_{pept_chain}.txt"), "w"
                ) as outfile:
                    outfile.write(
                        "\n".join(
                            [
                                f"selected_mode: {selected_mode}",
                                f"reference_length: {metrics['reference_length']}",
                                f"model_length: {metrics['model_length']}",
                                f"aligned_length: {metrics['aligned_length']}",
                                f"aligned_fraction_shorter: {metrics['aligned_fraction_shorter']:.3f}",
                                f"identity: {metrics['identity']:.3f}",
                                f"gap_count: {metrics['gap_count']}",
                                f"gap_fraction: {metrics['gap_fraction']:.3f}",
                                f"mismatch_count: {metrics['mismatch_count']}",
                                f"longest_contiguous_block: {metrics['longest_contiguous_block']}",
                                f"excluded_reference_residues: {metrics['excluded_reference_residues']}",
                                f"excluded_model_residues: {metrics['excluded_model_residues']}",
                            ]
                        )
                    )
            self.rmslst[pept_chain] = self.trajectory.rmsd_to_reference(
                ref_pep_stc, self_pep_stc
            )
            rmsds = [header.rmsd for header in self.medoids.headers]
            results = {
                "rmsds_all": [header.rmsd for header in self.trajectory.headers],
                "rmsds_filtered": [
                    header.rmsd for header in self.filtered_trajectory.headers
                ],
                "rmsds_medoids": rmsds,
            }
            results["lowest_all"] = sorted(results["rmsds_all"])[0]
            results["lowest_filtered"] = sorted(results["rmsds_filtered"])[0]
            results["lowest_medoids"] = sorted(results["rmsds_medoids"])[0]
            # Saving rmsd results
            if save:
                with open(
                    os.path.join(odir, "lowest_rmsds_%s.txt" % pept_chain), "w"
                ) as outfile:
                    outfile.write(
                        "lowest_all; lowest_filtered; lowest_medoids\n {0};{1};{2}".format(
                            results["lowest_all"],
                            results["lowest_filtered"],
                            results["lowest_medoids"],
                        )
                    )
                for _type in ["all", "filtered", "medoids"]:
                    with open(
                        os.path.join(odir, f"{_type}_rmsds_{pept_chain}.txt"), "w"
                    ) as outfile:
                        for rmsd in results["rmsds_" + _type]:
                            outfile.write(str(rmsd) + ";\n")
            all_results[pept_chain] = results
        if save and peptide_alignment_summary:
            self._save_peptide_alignment_summary(
                os.path.join(odir, "peptide_alignment_summary.csv"),
                peptide_alignment_summary,
            )
        if all_results:
            logger.info(module_name=_name, msg="RMSD successfully saved")
        else:
            logger.warning(
                module_name=_name,
                msg="Target RMSD alignment succeeded, but no peptide RMSD results were saved.",
            )
        return all_results

    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        logger.debug(module_name=_name, msg="Scoring results")
        # Filtering the trajectory
        self.filtered_trajectory, self.filtered_ndx = Filter(
            self.trajectory, n_filtered
        ).cabs_filter()
        # Clustering the trajectory
        self.medoids, self.clusters_dict, self.clusters = Clustering(
            self.filtered_trajectory,
            "chain "
            + ",".join(
                self.initial_complex.peptide_chains,
            ),
        ).cabs_clustering(
            number_of_medoids=number_of_medoids,
            number_of_iterations=number_of_iterations,
        )
        logger.info(module_name=_name, msg="Scoring results successful")

    def draw_plots(self, plots_dir=None, colors=None):
        logger.debug(module_name=_name, msg="Drawing plots")
        super(DockTask, self).draw_plots()
        # set the plots dir
        if plots_dir is None:
            pltdir = os.path.join(self.work_dir, "plots")
            try:
                os.mkdir(pltdir)
            except OSError:
                pass
        else:
            pltdir = plots_dir
        logger.log_file(module_name=_name, msg="Saving plots to %s" % pltdir)

        graph_RMSF(
            self.trajectory,
            self.initial_complex.protein_chains,
            os.path.join(pltdir, "RMSF"),
        )

        # RMSD-based graphs
        if self.reference:
            logger.log_file(module_name=_name, msg="Saving RMSD plots")
            for k, rmslst in self.rmslst.items():
                plot_E_RMSD(
                    [self.trajectory, self.filtered_trajectory],
                    [rmslst, rmslst[self.filtered_ndx,]],
                    ["all models", "top 1000 models"],
                    os.path.join(pltdir, "E_RMSD_%s" % k),
                )
                plot_RMSD_N(
                    rmslst.reshape(self.replicas, -1),
                    os.path.join(pltdir, "RMSD_frame_%s" % k),
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
                colors=colors,
            )
        logger.info(module_name=_name, msg="Plots successfully saved")

    @staticmethod
    def _add_cmaps(mk_cmap_output):
        # breakpoint()
        map_1, map_2 = mk_cmap_output
        return ContactMap(
            map_1.cmtx + map_2.cmtx, map_1.s1, map_2.s2, map_1.n + map_2.n
        )

    def mk_cmaps(
        self,
        ca_traj,
        meds,
        clusts,
        top1k_inds,
        thr,
        thra,
        plots_dir,
        colors=DEFAULT_COLORS,
    ):
        sc_traj_full, sc_med, cmapdir = super(DockTask, self).mk_cmaps(
            ca_traj, meds, clusts, top1k_inds, thr, thra, plots_dir
        )

        thrt = thra if self.aa_rebuild else thr

        sc_traj_1k = sc_traj_full.reshape(1, -1, len(ca_traj.template), 3)[
            :, top1k_inds, :, :
        ]

        rchs = self.initial_complex.protein_chains
        lchs = self.initial_complex.peptide_chains

        # Map to sequential chain IDs for the sequential trajectory template ca_traj.template
        all_orig = (self.initial_complex.protein_chains or []) + (self.initial_complex.peptide_chains or [])
        all_temp = sorted(list(ca_traj.template.list_chains().keys()))
        chain_map = {}
        for i, orig_ch in enumerate(all_orig):
            if i < len(all_temp):
                chain_map[orig_ch] = all_temp[i]

        rchs_seq = [chain_map.get(ch, ch) for ch in rchs]
        lchs_seq = [chain_map.get(ch, ch) for ch in lchs]

        targ_cmf = ContactMapFactory(rchs_seq, rchs_seq, ca_traj.template)

        cmfs = {lig: ContactMapFactory(rchs_seq, lig, ca_traj.template) for lig in lchs_seq}
        # cmap10ktarg = self._add_cmaps(targ_cmf.mk_cmap(sc_traj_full, thr))
        cmap10ktarg = reduce(operator.add, targ_cmf.mk_cmap(sc_traj_full, thr))
        cmap10ktarg.zero_diagonal()
        cmap10ktarg.save_all(
            cmapdir + "/target_all", break_long_x=0, norm_n=True, colors=colors
        )

        for lig, cmf in cmfs.items():
            cmaps = cmf.mk_cmap(sc_traj_full, thr)
            for n, cmap in enumerate(cmaps):
                cmap.save_all(
                    cmapdir + "/replica_%i_ch_%s" % (n + 1, lig),
                    norm_n=True,
                    colors=colors,
                )
            cmap10k = reduce(operator.add, cmaps)
            cmap10k.save_all(cmapdir + "/all_ch_%s" % lig, norm_n=True, colors=colors)
            cmap10k.save_histo(plots_dir + "/all_contacts_histo_%s" % lig)
            cmap1k = cmf.mk_cmap(sc_traj_1k, thr)[0]
            cmap1k.save_all(
                cmapdir + "/top1000_ch_%s" % lig, norm_n=True, colors=colors
            )
            if self.aa_rebuild:
                cmft = ContactMapFactory(rchs, rchs, self.medoids.template)
            else:
                cmft = cmf
            cmaps_top = cmft.mk_cmap(sc_med, thrt)
            for n, cmap in enumerate(cmaps_top):
                cmap.save_all(
                    cmapdir + "/top_%i_ch_%s" % (n + 1, lig), norm_n=True, colors=colors
                )
            for cn, clust in clusts.items():
                ccmap = cmf.mk_cmap(sc_traj_1k, thr, frames=clust)[0]
                ccmap.save_all(
                    cmapdir + "/cluster_%i_ch_%s" % (cn, lig),
                    norm_n=True,
                    colors=colors,
                )

    def parse_reference(self, ref, pdb_cache):
        try:
            source, rec, pep = ref.split(":")
            self.reference = (
                pdblib.Pdb(
                    source,
                    pdb_cache=pdb_cache,
                    selection="name CA and (chain %s)" % ",".join(rec + pep),
                    no_exit=True,
                    verify=True,
                ).atoms,
                rec,
                pep,
            )
            super(DockTask, self).parse_reference(ref, pdb_cache)
            if len(self.initial_complex.peptide_chains) != len(self.reference[2]):
                raise ValueError
            logger.info(_name, f"Reference {ref} loaded.")
        except (ValueError, pdblib.Pdb.InvalidPdbInput):
            logger.warning(_name, f"Invalid reference {ref}")
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
            mode=self.category_mode,
            peptides=self.peptides,
            replicas=self.replicas,
            separation=self.separation,
            insertion_attempts=self.insertion_attempts,
            insertion_clash=self.insertion_clash,
            work_dir=self.work_dir,
            receptor_ss=self.receptor_ss,
            pdb_cache=self.pdb_cache,
            save_initial_pdb=self.save_initial_pdb,
            json_output=self.json_output,
            predict_peptide_structure=self.peptide_structure_prediction,
            cg2all_env_prefix=self.cg2all_env_prefix,
            sc=self.write_sc_start_pdbs,
        )

        # score results logic
        self.score_results_dict = dict()

        if isinstance(self.pdb_output, str):
            self.pdb_output = self.pdb_output.replace("F", "")

    def score_results(self, n_filtered, number_of_medoids, number_of_iterations):
        # Clustering the trajectory
        chs = self.initial_complex.protein_chains or self.initial_complex.peptide_chains
        clst = Clustering(
            self.trajectory, "chain " + ",".join(chs)
        )
        self.medoids, self.clusters_dict, self.clusters = clst.cabs_clustering(
            number_of_medoids=number_of_medoids,
            number_of_iterations=number_of_iterations,
        )
        self.rmslst = {"".join(chs): clst.distance_matrix[0]}

    def load_output(self, *args, **kwargs):
        ret = super(FlexTask, self).load_output(*args, **kwargs)
        ret.number_of_peptides = 0
        return ret

    def calculate_rmsd(self, reference_pdb=None, save=True):
        logger.debug(module_name=_name, msg="RMSD calculations starting...")
        odir = None
        if save:
            odir = os.path.join(self.work_dir, "output_data")
            try:
                os.mkdir(odir)
            except OSError:
                pass

        chs_ids = self.trajectory.tmp_target_chs
        ref_trg_stc, self_trg_stc, trg_aln = self.trajectory.align_to(
            self.reference[0],
            self.reference[1],
            chs_ids,
            align_mth=self.align,
            kwargs=self.align_options,
        )
        if save:
            sfname = os.path.join(self.work_dir, "output_data", "reference_alignment")
            paln_trg = sfname + "_target.csv"
            save_csv(paln_trg, ("reference", "template"), trg_aln)
        self.rmslst[chs_ids] = self.trajectory.rmsd_to_reference(
            ref_trg_stc, self_trg_stc
        )
        rmsds = [header.rmsd for header in self.medoids.headers]
        results = {
            "rmsds_all": [header.rmsd for header in self.trajectory.headers],
            "rmsds_medoids": rmsds,
        }
        results["lowest_all"] = sorted(results["rmsds_all"])[0]
        results["lowest_medoids"] = sorted(results["rmsds_medoids"])[0]
        # Saving rmsd results
        if save:
            with open(
                os.path.join(odir, "lowest_rmsds_%s.txt" % chs_ids), "w"
            ) as outfile:
                outfile.write(
                    "lowest_all; lowest_medoids\n {0};{1}".format(
                        results["lowest_all"], results["lowest_medoids"]
                    )
                )
            for _type in ["all", "medoids"]:
                with open(
                    os.path.join(odir, f"{_type}_rmsds_{chs_ids}.txt"), "w"
                ) as outfile:
                    for rmsd in results["rmsds_" + _type]:
                        outfile.write(str(rmsd) + ";\n")
        logger.info(module_name=_name, msg="RMSD successfully saved")
        return {chs_ids: results}

    def mk_cmaps(
        self,
        ca_traj,
        meds,
        clusts,
        top1k_inds,
        thr,
        thra,
        plots_dir,
        colors=DEFAULT_COLORS,
    ):
        sc_traj_full, sc_med, cmapdir = super(FlexTask, self).mk_cmaps(
            ca_traj, meds, clusts, top1k_inds, thr, thra, plots_dir
        )

        thrt = thra if self.aa_rebuild else thr

        rchs = self.initial_complex.protein_chains or self.initial_complex.peptide_chains

        # Map to sequential chain IDs for the sequential trajectory template ca_traj.template
        all_orig = (self.initial_complex.protein_chains or []) + (self.initial_complex.peptide_chains or [])
        all_temp = sorted(list(ca_traj.template.list_chains().keys()))
        chain_map = {}
        for i, orig_ch in enumerate(all_orig):
            if i < len(all_temp):
                chain_map[orig_ch] = all_temp[i]

        rchs_seq = [chain_map.get(ch, ch) for ch in rchs]

        cmf = ContactMapFactory(rchs_seq, rchs_seq, ca_traj.template)
        if self.aa_rebuild:
            cmft = ContactMapFactory(rchs, rchs, self.medoids.template)
        else:
            cmft = cmf
        cmap_all = reduce(operator.add, cmf.mk_cmap(sc_traj_full, thr))

        topscms = cmft.mk_cmap(sc_med, thrt)
        cmaptop = reduce(operator.add, topscms)

        for cmap, fname in zip(
            (cmap_all, cmaptop) + tuple(topscms),
            ("all", "top10")
            + tuple(["top_%i" % (i + 1) for i, dummy in enumerate(topscms)]),
        ):
            cmap.zero_diagonal()
            cmap.save_all(
                cmapdir + "/" + fname, break_long_x=0, norm_n=True, colors=colors
            )

    def parse_reference(self, ref, pdb_cache):
        try:
            try:
                if ":" in str(ref):
                    source, trg_chids = ref.split(":")
                else:
                    source = ref
                    trg_chids = None
                ref_stc = pdblib.Pdb(
                    source,
                    pdb_cache=pdb_cache,
                    selection="name CA",
                    no_exit=True,
                    verify=True,
                ).atoms
                if trg_chids is None:
                    trg_chids = "".join(sorted(set([i.chid for i in ref_stc])))
                self.reference = (ref_stc, trg_chids)
                super(FlexTask, self).parse_reference(ref, pdb_cache)
                logger.info(_name, f"Reference {ref} loaded.")
            except AttributeError:  # if ref is None it has no split mth
                chs = self.initial_complex.protein_chains or self.initial_complex.peptide_chains
                ref_stc = self.initial_complex.select(
                    "name CA and (chain %s)" % ",".join(chs)
                )
                if self.initial_complex.protein:
                    ref_stc.update_ids(self.initial_complex.protein.old_ids)
                self.reference = (
                    ref_stc,
                    "".join(sorted(set([i.chid for i in ref_stc]))),
                )
                logger.info(_name, "Input loaded as reference.")
        except (pdblib.Pdb.InvalidPdbInput, ValueError):
            logger.warning(_name, f"Invalid reference {ref}")

    def draw_plots(self, plots_dir=None, colors=DEFAULT_COLORS):
        super(FlexTask, self).draw_plots()
        # set the plots dir
        if plots_dir is None:
            pltdir = os.path.join(self.work_dir, "plots")
            try:
                os.mkdir(pltdir)
            except OSError:
                pass
        else:
            pltdir = plots_dir

        graph_RMSF(
            self.trajectory,
            self.initial_complex.protein_chains or self.initial_complex.peptide_chains,
            os.path.join(pltdir, "RMSF"),
            fmt=self.image_file_format,
        )

        # RMSD-based graphs
        if self.reference:
            for k, rmslst in self.rmslst.items():
                plot_E_RMSD(
                    [self.trajectory],
                    [rmslst],
                    ["all models"],
                    os.path.join(pltdir, "E_RMSD_%s" % k),
                    self.image_file_format,
                    interaction=False,
                )
                plot_RMSD_N(
                    rmslst.reshape(self.replicas, -1),
                    os.path.join(pltdir, "RMSD_frame_%s" % k),
                    self.image_file_format,
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
                colors=colors,
            )
