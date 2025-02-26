"""
Classes Protein, Peptide, ProteinComplex - prepares initial complex.
"""

import re
import os
import json
from collections import OrderedDict
from copy import deepcopy
from string import ascii_uppercase
from math import exp

from CABS import utils
from CABS import logger
from CABS.pdblib import Pdb
from CABS.atom import Atoms
from CABS.vector3d import Vector3d
from CABS import randinit

_name = 'Protein'


class Protein(Atoms):
    """
    Class for the protein molecule.
    """
    NSP3_MODEL_PATH = ''

    def __init__(self, source, flexibility=None, exclude=None, weights=None, plddt=None, work_dir='.', receptor_ss=None,
                 pdb_cache=None, category=None, save_initial_pdb=False, predict_peptide_structure=False):

        Atoms.__init__(self)

        logger.info(module_name=_name, msg="Loading %s as input protein" % source)

        # Only happens if user explicitly wants to predict peptide structure
        if predict_peptide_structure:
            try:
                self.atoms = randinit.RandomInitialStructure(source).pdb
            except Exception as e:
                logger.exit_program(
                    module_name=_name,
                    msg=f'Invalid input {source} for peptide structure prediction',
                    exc=e
                )
            predictor = None
            if ':' not in source:
                if self.NSP3_MODEL_PATH:
                        try:
                            from CABS.secstrpredictor import SecStrPredictor
                            predictor = SecStrPredictor(self.NSP3_MODEL_PATH)
                        except ImportError as e:
                            logger.warning(
                                module_name=_name,
                                msg=f'NetSurfP-3.0 library or its dependencies are missing: {str(e)}'
                            )
                        except Exception as e:
                            logger.warning(
                                module_name=_name,
                                msg=f'Cannot load NetSurfP-3.0 model: {str(e)}'
                            )
                else:
                    logger.warning(
                        module_name=_name,
                        msg='NSP3 model path not provided.'
                    )

                if not predictor:
                    logger.warning(
                        module_name=_name,
                        msg='Secondary structure prediction will not be performed.'
                    )

            if predictor:
                logger.info(
                    module_name=_name,
                    msg='Running secondary structure prediction for the peptide using NetSurfP-3.0.'
                )
                try:
                    sec_str = predictor.predict_q3(sequence_to_predict=source)
                    ss = OrderedDict((a.resid_id(), sec_str[i]) for i, a in enumerate(self.atoms))
                    logger.info(
                        module_name=_name,
                        msg='Secondary structure prediction for the peptide successful.'
                    )
                except Exception as e:
                    logger.warning(
                        module_name=_name,
                        msg=f'Secondary structure prediction for the peptide failed: {str(e)}'
                    )
                    CABS_SS = 'CHTE'
                    ss = OrderedDict((a.resid_id(), CABS_SS[int(a.occ) - 1]) for a in self.atoms)
            else:
                CABS_SS = 'CHTE'
                ss = OrderedDict((a.resid_id(), CABS_SS[int(a.occ) - 1]) for a in self.atoms)

        # This is the default case, the same as before
        else:
            try:
                self.atoms = randinit.RandomInitialStructure(source).pdb
                CABS_SS = 'CHTE'
                ss = OrderedDict((a.resid_id(), CABS_SS[int(a.occ) - 1]) for a in self.atoms)

            except:
                pdb = Pdb(source=source, selection='name CA', pdb_cache=pdb_cache)
                self.atoms = pdb.atoms.models()[0]
                ss = pdb.dssp(work_dir=work_dir)
                if save_initial_pdb:
                    pdb.save_initial_pdb(work_dir=work_dir)

        if receptor_ss:
            logger.info('Running manual assignment of receptor\'s II structure.')
            try:
                ss = ReceptorSS(current_ss=ss, receptor_ss=receptor_ss).ss
            except InvalidReceptorSS:
                logger.warning(msg='Invalid data for --receptor-ss option')

        # setup plddt
        if plddt:
            if plddt.lower() == 'pdb' or plddt.lower() == 'bf':
                for a in self.atoms:
                    a.plddt = a.bfac / 100
            else:
                try:
                    d, de = self.read_plddt(plddt)
                    self.atoms.update_plddt(d, de)
                except IOError:
                    logger.warning(_name, 'Could not read pLDDT file: %s' % plddt)
                    logger.warning(_name, 'Using default plddt(1.0) for all residues.')
                    self.atoms.set_plddt(1.0)
                except Exception as e:
                    logger.warning(_name, '%s' % e)
                    logger.warning(_name, 'Using default plddt(1.0) for all residues.')
                    self.atoms.set_plddt(1.0)
        else:
            self.atoms.set_plddt(1.0)

        # setup flexibility
        if flexibility:
            try:
                flex = float(flexibility)
                self.atoms.set_flexibility(flex)
            except ValueError:
                if flexibility.lower() == 'bf':
                    for a in self.atoms:
                        a.flexibility = a.bfac
                elif flexibility.lower() == 'bfi':
                    for a in self.atoms:
                        if a.bfac > 1.:
                            a.flexibility = 0.
                        elif a.bfac < 0.:
                            a.flexibility = 1.
                        else:
                            a.flexibility = 1. - a.flexibility
                elif flexibility.lower() == 'bfg':
                    for a in self.atoms:
                        if a.bfac < 0.:
                            a.flexibility = 1.
                        else:
                            a.flexibility = exp(-0.5 * a.bfac ** 2)
                else:
                    try:
                        d, de = self.read_flexibility(flexibility)
                        self.atoms.update_flexibility(d, de)
                    except IOError:
                        logger.warning(_name, 'Could not read flexibility file: %s' % flexibility)
                        logger.warning(_name, 'Using default flexibility(1.0) for all residues.')
                        self.atoms.set_flexibility(1.0)
        else:
            self.atoms.set_flexibility(1.0)

        # setup excluding
        self.exclude = {}
        if exclude:
            for s in exclude:
                words = s.split('@')
                if len(words) == 1:
                    key = 'ALL'
                else:
                    key = utils.pep2pep1(words[-1])
                if key in self.exclude:
                    self.exclude[key] += '+' + words[0]
                else:
                    self.exclude[key] = words[0]

            for k, v in self.exclude.items():
                self.exclude[k] = []
                for word in v.split('+'):
                    if ':' in word:
                        if '-' in word:
                            beg, end = word.split('-')
                            self.exclude[k].extend(self.atoms.atom_range(beg, end))
                        else:
                            self.exclude[k].append(word)
                    else:
                        chains = re.sub(r'[^%s]*' % word, '', ascii_uppercase)
                        self.exclude[k].extend(a.resid_id() for a in self.atoms.select('chain %s' % chains))

        self.atoms.update_sec(ss)

        if work_dir and logger.output_ss():
            output_ss = os.path.join(work_dir, 'output_data', 'ss.txt')
            odir = os.path.dirname(output_ss)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            logger.to_file(
                filename=output_ss,
                content=str(''.join(ss.values())),
                msg='Saving secondary structure output to %s' % output_ss
            )

        # setup categories
        if category:
            try:
                d, de = self.read_category(category)
                self.atoms.update_category(d, de)
            except IOError:
                logger.warning(_name, 'Could not read category file: %s' % category)
                logger.warning(_name, 'Using default categories based on pLDDT and SS for all residues.')
                self.atoms.determine_category()
        else:
            self.atoms.determine_category()

        self.old_ids = self.atoms.fix_broken_chains()

        todrop = [c for c, l in self.list_chains().items() if l < 4]
        if todrop:
            self.atoms = self.atoms.drop('chain ' + ','.join(todrop))

        self.new_ids = {v: k for k, v in self.old_ids.items()}

        for key, val in self.exclude.items():
            self.exclude[key] = [self.new_ids[r] for r in val]

        # setup rmsd_weights
        self.weights = None
        if weights and weights.lower() == 'flex':
            self.weights = [a.flexibility for a in self.atoms]
        if weights and weights.lower() == 'ss':
            self.weights = [(a.occ + 1.) % 2 for a in self.atoms]
        if weights and (weights.lower() == 'off' or weights.lower() == 'gauss'):
            self.weights = [1.0] * len(self.atoms)
        else:
            try:
                default = 1.0
                self.weights = []
                weights_dict = {}
                with open(weights, 'r') as _file:
                    for line in _file:
                        k, v = line.split()[:2]
                        weights_dict[k] = v

                if 'default' in weights_dict:
                    default = float(weights_dict['default'])

                for a in self.atoms:
                    w = weights_dict.get(a.resid_id())
                    w = float(w) if w else default
                    self.weights.append(w)

            except (IOError, ValueError):
                logger.warning(_name, 'Could not read weights file: %s' % weights)
                logger.warning(_name, 'Using default weights(1.0) for all atoms.')
                self.weights = [1.0] * len(self.atoms)

        self.center = self.cent_of_mass()
        self.dimension = self.max_dimension()
        self.patches = {}

    def convert_patch(self, location):
        if location not in self.patches:
            chains = {}
            for res in [self.new_ids[r] for r in location.split('+')]:
                num, chid = res.split(':')
                if chid in chains:
                    chains[chid].append(num)
                else:
                    chains[chid] = [num]
            s = " or ".join(["(chain " + ch + " and resnum " + ",".join(chains[ch]) + ")" for ch in chains])
            patch = self.select(s)
            self.patches[location] = (patch.cent_of_mass() - self.center).norm()
        return self.patches[location]

    @staticmethod
    def read_flexibility(filename):

        key = r'[0-9A-Z]+:[A-Z]'
        val = r'[0-9.]+'

        patt_range = re.compile('(%s) *-* *(%s) +(%s)' % (key, key, val))
        patt_single = re.compile('(%s) +(%s)' % (key, val))

        with open(filename) as f:
            d = {}
            def_val = 1.0
            for line in f:
                if re.search('default', line):
                    def_val = float(line.split()[-1])
                else:
                    match = re.search(patt_range, line)
                    if match:
                        n1, c1 = match.group(1).split(':')
                        n2, c2 = match.group(2).split(':')
                        n1 = int(n1)
                        n2 = int(n2)
                        if c1 != c2 or n1 > n2:
                            raise Exception('Invalid range: \'%s\' in file: %s!!!' % (line, filename))
                        for i in range(n1, n2 + 1):
                            d[str(i) + ':' + c1] = float(match.group(3))
                    else:
                        match = re.search(patt_single, line)
                        if match:
                            d[match.group(1)] = float(match.group(2))
                        else:
                            raise Exception('Invalid syntax in flexibility file!!!')
            return d, def_val

    @staticmethod
    def validate_plddt_file(file_path):
        """
        Validates a pLDDT file, ensuring it's either a valid JSON or TSV format containing pLDDT values.

        Parameters:
            file_path (str): Path to the file to validate.

        Returns:
            str: "valid_json" or "valid_tsv" if the file is valid.
            Raises a ValueError if the file is invalid.
        """
        try:
            # Try to validate as JSON
            with open(file_path, 'r') as file:
                json_dict = json.load(file)
                if "plddt" not in json_dict:
                    raise ValueError("Validation failed: Missing 'plddt' field in JSON.")
                plddt_values = json_dict["plddt"]
                if not isinstance(plddt_values, list):
                    raise ValueError("Validation failed: 'plddt' field must be a list.")
                if not all(isinstance(value, (int, float)) and 0 <= value <= 100 for value in plddt_values):
                    raise ValueError("Validation failed: All pLDDT scores must be numbers between 0 and 100.")
            return "valid_json"
        except json.JSONDecodeError:
            # Not JSON, try TSV
            try:
                with open(file_path, 'r') as file:
                    # Skip the header line if present
                    next(file)
                    for line in file:
                        columns = line.strip().split('\t')
                        if len(columns) < 2:
                            raise ValueError("Validation failed: TSV file must have at least two columns.")
                        float(columns[1])  # Ensure the second column is a valid number
                return "valid_tsv"
            except Exception as e:
                raise ValueError("Validation failed: File is neither valid JSON nor TSV.") from e

    def read_plddt(self, filename):
        """
        Reads pLDDT values from a validated JSON or TSV file and maps them to residues.

        Parameters:
            filename (str): Path to the pLDDT file.

        Returns:
            dict: Mapping of residue IDs to pLDDT values.
            float: Default pLDDT value.
        """

        d = {}
        def_val = 1.0
        plddt_values = []

        file_type = self.validate_plddt_file(filename)

        try:
            if file_type == "valid_json":
                with open(filename, 'r') as file:
                    json_dict = json.load(file)
                    plddt_values = [float(entry) / 100 for entry in json_dict['plddt']]
            elif file_type == "valid_tsv":

                with open(filename, 'r') as file:
                    next(file)
                    for line in file:
                        columns = line.strip().split('\t')
                        plddt_values.append(float(columns[1]) / 100)
        except Exception as e:
            raise ValueError("Error while reading pLDDT values.") from e

        for i, atom in enumerate(self.atoms):
            if i < len(plddt_values):
                d[atom.resid_id()] = plddt_values[i]
            else:
                raise ValueError("Mismatch between the number of residues and pLDDT values.")

        return d, def_val

    @staticmethod
    def read_category(filename):

        key = r'[0-9A-Z]+:[A-Z]'
        val = r'[0-9.]+'

        patt_range = re.compile('(%s) *-* *(%s) +(%s)' % (key, key, val))
        patt_single = re.compile('(%s) +(%s)' % (key, val))

        with open(filename) as f:
            d = {}
            def_val = 0.0
            for line in f:
                if re.search('default', line):
                    def_val = float(line.split()[-1])
                else:
                    match = re.search(patt_range, line)
                    if match:
                        n1, c1 = match.group(1).split(':')
                        n2, c2 = match.group(2).split(':')
                        n1 = int(n1)
                        n2 = int(n2)
                        if c1 != c2 or n1 > n2:
                            raise Exception('Invalid range: \'%s\' in file: %s!!!' % (line, filename))
                        for i in range(n1, n2 + 1):
                            d[str(i) + ':' + c1] = float(match.group(3))
                    else:
                        match = re.search(patt_single, line)
                        if match:
                            d[match.group(1)] = float(match.group(2))
                        else:
                            raise Exception('Invalid syntax in flexibility file!!!')
            return d, def_val

    def generate_restraints(self, mode, gap, min_d, max_d):
        gap = int(gap)
        min_d = float(min_d)
        max_d = float(max_d)
        restr = []
        _len = len(self.atoms)

        if mode in ['category', 'plddt']:
            for i in range(_len):
                a1 = self.atoms[i]
                for j in range(i + gap, _len):
                    a2 = self.atoms[j]
                    d = (a1.coord - a2.coord).length()
                    if min_d < d < max_d:
                        sum_of_category = a1.category + a2.category
                        if sum_of_category < 4:
                            w = 0.0
                        elif sum_of_category == 4:
                            w = 0.5
                        else:
                            w = 1.0
                        if w:
                            restr.append('%s %s %f %f' % (a1.resid_id(), a2.resid_id(), d, w))

        elif mode in ['min', 'max', 'mean', 'plddt1', 'plddt2']:
            for i in range(_len):
                a1 = self.atoms[i]
                for j in range(i + gap, _len):
                    a2 = self.atoms[j]
                    d = (a1.coord - a2.coord).length()
                    if min_d < d < max_d:
                        if a1.plddt > a2.plddt:
                            w_bigger = a1.plddt
                            w_smaller = a2.plddt
                        else:
                            w_bigger = a2.plddt
                            w_smaller = a1.plddt
                        if mode == 'max':
                            w = w_bigger if w_bigger > 0.5 else 0
                        elif mode == 'mean':
                            w = (w_bigger + w_smaller) / 2 if (w_bigger + w_smaller) / 2 > 0.5 else 0
                        elif mode == 'plddt1':
                            w = 1.0 if (w_smaller > 0.5 or w_smaller > 0.5) else 0
                        elif mode == 'plddt2':
                            w = 1.0 if w_bigger > 0.5 and w_bigger > 0.5 else 0
                        else:
                            w = w_smaller if w_smaller > 0.5 else 0
                        if w:
                            restr.append('%s %s %f %f' % (a1.resid_id(), a2.resid_id(), d, w))

        else:
            for i in range(_len):
                a1 = self.atoms[i]
                ssi = int(a1.occ) % 2
                if mode == 'ss2' and ssi:
                    continue
                for j in range(i + gap, _len):
                    a2 = self.atoms[j]
                    ssj = int(a2.occ) % 2
                    if (mode == 'ss2' and ssj) or (mode == 'ss1' and ssi * ssj):
                        continue
                    d = (a1.coord - a2.coord).length()
                    if min_d < d < max_d:
                        if a1.flexibility < a2.flexibility:
                            w = a1.flexibility
                        else:
                            w = a2.flexibility
                        if w:
                            restr.append('%s %s %f %f' % (a1.resid_id(), a2.resid_id(), d, w))
        return restr

    def generate_backbone_restraints(self, cyclic_chains):
        restr = []
        for chain in cyclic_chains:
            first_res = None
            for atom in self.atoms:
                if atom.chid == chain:
                    first_res = atom.resid_id()
                    break
            last_res = None
            for atom in reversed(self.atoms):
                if atom.chid == chain:
                    last_res = atom.resid_id()
                    break
            if first_res and last_res:
                restr.append('%s %s 3.8 1.0' % (first_res, last_res))
            else:
                logger.warning(module_name=_name, msg='Cyclic backbone could not be created in chain %s' % chain)

        return restr

    def generate_disulfide_restraints(self, disulfide_bonds):
        restr = []
        for bond in disulfide_bonds:
            res1 = None
            res2 = None
            for atom in self.atoms:
                if atom.resid_id() == bond[0] and atom.resname == 'CYS':
                    res1 = atom.resid_id()
                elif atom.resid_id() == bond[1] and atom.resname == 'CYS':
                    res2 = atom.resid_id()
            if res1 and res2:
                restr.append('%s %s 2.0 1.0' % (res1, res2))
            else:
                logger.warning(module_name=_name, msg='Disulfide bond between residues %s %s could not be created' % (bond[0], bond[1]))
        return restr

    def calculate_distances(self):
        """
        Generate a matrix of distances between each C-alpha in the protein (server uses this)
        :return: NxN matrix of distances (dict of dicts {'1:A': {'1:A' :20, '2:A': 30}, ...} )
        """
        out = OrderedDict()
        _len = len(self.atoms)
        for i in range(_len):
            a1 = self.atoms[i]
            out[a1.resid_id()] = {}
            for j in range(_len):
                a2 = self.atoms[j]
                d = (a1.coord - a2.coord).length()
                out[a1.resid_id()][a2.resid_id()] = "%6.4f" % d
        return out


class Peptide(Atoms):
    """
    Class for the peptides.
    """

    def __init__(self, source, conformation, location, work_dir='.', pdb_cache=None):
        logger.info(
            module_name=_name,
            msg='Loading ligand: {}, conformation - {}, location - {}'.format(
                source, conformation, location
            )
        )
        try:
            pdb = Pdb(source=source, selection='name CA', pdb_cache=pdb_cache, no_exit=True)
            atoms = pdb.atoms.models()[0]
            atoms.update_sec(pdb.dssp(work_dir=work_dir))
        except Pdb.InvalidPdbInput:
            atoms = Atoms(source)
        atoms.set_bfac(0.0)
        self.conformation = conformation
        self.location = location
        Atoms.__init__(self, atoms)


class ProteinComplex(Atoms):
    """
    Class that assembles the initial complex.
    """

    def __init__(self, protein, flexibility, exclude, weights, plddt, category, peptides, replicas,
                 separation, insertion_attempts, insertion_clash, work_dir, receptor_ss, pdb_cache, save_initial_pdb,
                 json_output=False, predict_peptide_structure=False):
        logger.debug(module_name=_name, msg="Preparing the complex")
        Atoms.__init__(self)

        self.protein = Protein(
            protein,
            flexibility=flexibility,
            exclude=exclude,
            weights=weights,
            plddt=plddt,
            category=category,
            work_dir=work_dir,
            receptor_ss=receptor_ss,
            pdb_cache=pdb_cache,
            save_initial_pdb=save_initial_pdb,
            predict_peptide_structure=predict_peptide_structure
        )
        self.chain_list = self.protein.list_chains()
        self.protein_chains = ''.join(self.chain_list.keys())
        self.old_ids = deepcopy(self.protein.old_ids)

        self.peptides = []
        self.peptide_chains = ''
        if peptides:
            taken_chains = self.protein_chains + 'X'
            for num, p in enumerate(peptides):
                peptide = Peptide(*p, work_dir=work_dir, pdb_cache=pdb_cache)
                if peptide[0].chid in taken_chains:
                    peptide.change_chid(peptide[0].chid, utils.next_letter(taken_chains))
                taken_chains += peptide[0].chid
                self.peptide_chains += peptide[0].chid
                self.peptides.append(peptide)
                update_dict = {}
                i = 1
                for atom in peptide:
                    update_dict[atom.resid_id()] = '%i:PEP%i' % (i, num + 1)
                    i += 1
                self.old_ids.update(update_dict)
                self.chain_list.update(peptide.list_chains())
        self.new_ids = {v: k for k, v in self.old_ids.items()}

        exclude = []
        for key, value in self.protein.exclude.items():
            if key == 'ALL':
                kword = 'PEP'
            else:
                kword = key
            keys = [v for k, v in self.new_ids.items() if re.search(kword, k)]
            exclude.extend((r1, r2) for r1 in keys for r2 in value)
        self.protein.exclude = list(set(exclude))

        for i in range(replicas):
            model = deepcopy(self.protein)
            model.set_model_number(i + 1)
            for peptide in self.peptides:
                for attempt in range(insertion_attempts):
                    self.insert_peptide(self.protein, peptide, separation)
                    if model.min_distance(peptide) > insertion_clash:
                        peptide = deepcopy(peptide)
                        peptide.set_model_number(i + 1)
                        model.atoms.extend(peptide)
                        break
                else:
                    raise Exception('Maximum number of attempts to insert peptide %s reached!!!' % peptide) #used to be peptide.name
            self.atoms.extend(model)
        logger.debug(module_name=_name, msg="Complex successfully created")

        if json_output:
            json_file = os.path.join(work_dir, 'output_data', 'atoms.json')
            odir = os.path.dirname(json_file)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            self.save_to_json(json_file)

        logger.debug(module_name=_name, msg="Atoms saved to JSON file")

    @staticmethod
    def insert_peptide(protein, peptide, separation):

        radius = 0.5 * protein.dimension + separation

        if peptide.location == 'keep':
            location = peptide.cent_of_mass()
        elif peptide.location == 'random':
            peptide.rotate_in_place(utils.random_rotation_matrix())
            location = Vector3d().random() * radius + protein.center
        else:
            location = protein.convert_patch(peptide.location) * radius + protein.center

        if peptide.conformation == 'random':
            peptide.random_conformation()

        peptide.move_to(location)


class InvalidReceptorSS(Exception):
    pass


class ReceptorSS:

    def __init__(self, current_ss, receptor_ss):

        self.ss = OrderedDict()
        chains = {}

        try:
            with open(receptor_ss) as f:
                for line in f:
                    chid, ss = line.replace(' ', '').split(':')
                    chains[chid] = ss
        except IOError:
            chains = dict(ch.split(':') for ch in receptor_ss.split('+'))

        if not current_ss:
            tmp = []
            for c in chains:
                s = list(chains[c])
                for l in range(1, len(s)+1):
                    st = str(l)+':'+str(c)
                    tmp.append((st, s[l-1]))
            current_ss = OrderedDict([(a[0], a[1]) for a in tmp])

        if not chains:
            raise InvalidReceptorSS

        current_chains = OrderedDict()
        for k, v in current_ss.items():
            chid = k.split(':')[-1]
            if chid not in current_chains:
                current_chains[chid] = OrderedDict()
            current_chains[chid][k] = v

        for chid in current_chains:
            if len(current_chains[chid]) == len(chains[chid]):
                self.ss.update(OrderedDict(zip(current_chains[chid].keys(), chains[chid])))
            else:
                raise InvalidReceptorSS


if __name__ == '__main__':
    ss = OrderedDict([
        ('1:A', 'C'),
        ('2:A', 'H'),
        ('3:A', 'H'),
        ('4:A', 'H'),
        ('5:A', 'C'),
        ('2:B', 'C'),
        ('3:B', 'E'),
        ('4:B', 'E'),
        ('5:B', 'C'),
        ('6:B', 'H'),
        ('7:B', 'H'),
        ('8:B', 'C'),
    ])

    r_ss = 'A:EEEEE+B:CHHCEEC'
    r = ReceptorSS(ss, r_ss)
    print(r.ss)
