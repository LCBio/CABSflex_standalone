"""
Classes Protein, Peptide, ProteinComplex - prepares initial complex.
"""

import re
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

    def __init__(self, source, flexibility=None, exclude=None, weights=None, plddt=None, work_dir='.', receptor_ss=None,
                 pdb_cache=None):

        Atoms.__init__(self)

        logger.info(module_name=_name, msg="Loading %s as input protein" % source)

        try:
            self.atoms = randinit.RandomInitialStructure(source).pdb
            CABS_SS = 'CHTE'
            ss = OrderedDict([(a.resid_id(), CABS_SS[int(a.occ) - 1]) for a in self.atoms])

        except:
            pdb = Pdb(source=source, selection='name CA', pdb_cache=pdb_cache)
            self.atoms = pdb.atoms.models()[0]
            ss = pdb.dssp(output=work_dir)

        if receptor_ss:
            logger.info('Running manual assignment of receptor\'s II structure.')
            try:
                ss = ReceptorSS(current_ss=ss, receptor_ss=receptor_ss).ss
            except InvalidReceptorSS:
                logger.warning(msg='Invalid data for --receptor-ss option')

        if plddt:
            if plddt.lower() == 'pdb':
                for a in self.atoms:
                    a.bfac = a.bfac / 100
            else:
                try:
                    d, de = self.read_plddt(plddt)
                    self.atoms.update_bfac(d, de)
                except IOError:
                    logger.warning(_name, 'Could not read pLDDT file: %s' % plddt)
                    logger.warning(_name, 'Using default flexibility(1.0) for all residues.')
                    self.atoms.set_bfac(1.0)
                except Exception as e:
                    logger.warning(_name, '%s' % e)
                    logger.warning(_name, 'Using default flexibility(1.0) for all residues.')
                    self.atoms.set_bfac(1.0)
        else:
            # setup flexibility
            if flexibility:
                try:
                    bfac = float(flexibility)
                    self.atoms.set_bfac(bfac)
                except ValueError:
                    if flexibility.lower() == 'bf':
                        pass
                    elif flexibility.lower() == 'bfi':
                        for a in self.atoms:
                            if a.bfac > 1.:
                                a.bfac = 0.
                            elif a.bfac < 0.:
                                a.bfac = 1.
                            else:
                                a.bfac = 1. - a.bfac
                    elif flexibility.lower() == 'bfg':
                        for a in self.atoms:
                            if a.bfac < 0.:
                                a.bfac = 1.
                            else:
                                a.bfac = exp(-0.5 * a.bfac ** 2)
                    else:
                        try:
                            d, de = self.read_flexibility(flexibility)
                            self.atoms.update_bfac(d, de)
                        except IOError:
                            logger.warning(_name, 'Could not read flexibility file: %s' % flexibility)
                            logger.warning(_name, 'Using default flexibility(1.0) for all residues.')
                            self.atoms.set_bfac(1.0)
            else:
                self.atoms.set_bfac(1.0)

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

        self.old_ids = self.atoms.update_sec(ss).fix_broken_chains()

        todrop = [c for c, l in self.list_chains().items() if l < 4]
        if todrop:
            self.atoms = self.atoms.drop('chain ' + ','.join(todrop))

        self.new_ids = {v: k for k, v in self.old_ids.items()}

        for key, val in self.exclude.items():
            self.exclude[key] = [self.new_ids[r] for r in val]

        # setup rmsd_weights
        self.weights = None
        if weights and weights.lower() == 'flex':
            self.weights = [a.bfac for a in self.atoms]
        if weights and weights.lower() == 'ss':
            self.weights = [(a.occ + 1.) % 2 for a in self.atoms]
        if weights and (weights.lower() == 'off' or weights.lower() == 'gauss'):
            self.weights = [1.0] * len(self.atoms)
        else:
            try:
                default = 1.0
                self.weights = []
                weights_dict = {}
                with open(weights, 'rb') as _file:
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

    def read_plddt(self, filename):
        d = {}
        def_val = 1.0
        plddt_values = []
        if filename.endswith('.json'):
            with open(filename, 'r') as file:
                json_dict = json.load(file)
                parsed = json_dict['plddt']
                for i, entry in enumerate(parsed):
                    plddt_values.append(float(entry) / 100)

        elif filename.endswith('.tsv'):
            with open(filename, 'r') as file:
                next(file)
                for i, line in enumerate(file):
                    columns = line.strip().split('\t')
                    plddt_values.append(float(columns[1]) / 100)

        else:
            raise Exception("pLLDT file format must end with '.json' or '.tsv'")

        if len(plddt_values) != len(self.atoms):
            raise Exception(f"Number of pLLDT values: {len(plddt_values)}, does not match"
                            f"the number of atoms: {len(self.atoms)} in the protein")
        else:
            for i, a in enumerate(self.atoms):
                d[a.resid_id()] = plddt_values[i]

        return d, def_val

    def generate_restraints(self, mode, gap, min_d, max_d):
        gap = int(gap)
        min_d = float(min_d)
        max_d = float(max_d)
        restr = []
        _len = len(self.atoms)

        if mode == 'category':
            for i in range(_len):
                a1 = self.atoms[i]
                for j in range(i + gap, _len):
                    a2 = self.atoms[j]
                    d = (a1.coord - a2.coord).length()
                    if min_d < d < max_d:
                        w = self.generate_weight(a1, a2)
                        if w:
                            restr.append('%s %s %f %f' % (a1.resid_id(), a2.resid_id(), d, w))

        elif mode in ['min', 'max', 'mean', 'plddt1', 'plddt2', 'plddt_bins']:
            for i in range(_len):
                a1 = self.atoms[i]
                for j in range(i + gap, _len):
                    a2 = self.atoms[j]
                    d = (a1.coord - a2.coord).length()
                    if min_d < d < max_d:
                        if a1.bfac > a2.bfac:
                            w_bigger = a1.bfac
                            w_smaller = a2.bfac
                        else:
                            w_bigger = a2.bfac
                            w_smaller = a1.bfac
                        if mode == 'min':
                            w = w_smaller if w_smaller > 0.5 else 0
                        elif mode == 'max':
                            w = w_bigger if w_bigger > 0.5 else 0
                        elif mode == 'mean':
                            w = (w_bigger + w_smaller) / 2 if (w_bigger + w_smaller) / 2 > 0.5 else 0
                        elif mode == 'plddt1':
                            w = 1.0 if (w_smaller > 0.5 or w_smaller > 0.5) else 0
                        elif mode == 'plddt2':
                            w = 1.0 if w_bigger > 0.5 and w_bigger > 0.5 else 0
                        elif mode == 'plddt_bins':
                            w = self.generate_bin_weight(a1, a2)
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
                        if a1.bfac < a2.bfac:
                            w = a1.bfac
                        else:
                            w = a2.bfac
                        if w:
                            restr.append('%s %s %f %f' % (a1.resid_id(), a2.resid_id(), d, w))
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

    @staticmethod
    def generate_weight(a1, a2):

        sum_of_restraint_strength = 0
        for atom in [a1, a2]:
            bfac = atom.bfac
            occ = atom.occ
            restraint_strength = 0
            if bfac < 0.5:
                restraint_strength -= 1
            elif bfac < 0.7:
                pass
            elif bfac < 0.9:
                restraint_strength += 1
            else:
                restraint_strength += 2

            if occ == 1:
                restraint_strength += 0
            elif occ == 3:
                restraint_strength += 1
            else:
                restraint_strength += 2

            if restraint_strength < 0:
                pass
            if restraint_strength > 3:
                sum_of_restraint_strength += 3
            else:
                sum_of_restraint_strength += restraint_strength

        if sum_of_restraint_strength < 4:
            return 0.0
        if sum_of_restraint_strength == 4:
            return 0.5
        else:
            return 1.0

    @staticmethod
    def generate_bin_weight(a1, a2):
        sum_of_restraint_strength = 0
        for atom in [a1, a2]:
            bfac = atom.bfac
            if bfac < 0.5:
                restraint_strength = 0
            elif bfac < 0.7:
                restraint_strength = 1
            elif bfac < 0.9:
                restraint_strength = 2
            else:
                restraint_strength = 3

            sum_of_restraint_strength += restraint_strength

        if sum_of_restraint_strength < 4:
            return 0.0
        if sum_of_restraint_strength == 4:
            return 0.5
        else:
            return 1.0


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
            atoms.update_sec(pdb.dssp(output=work_dir))
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

    def __init__(self, protein, flexibility, exclude, weights, plddt, peptides, replicas,
                 separation, insertion_attempts, insertion_clash, work_dir, receptor_ss, pdb_cache):
        logger.debug(module_name=_name, msg="Preparing the complex")
        Atoms.__init__(self)

        self.protein = Protein(
            protein,
            flexibility=flexibility,
            exclude=exclude,
            weights=weights,
            plddt=plddt,
            work_dir=work_dir,
            receptor_ss=receptor_ss,
            pdb_cache=pdb_cache,
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
                self.old_ids.update({atom.resid_id(): '%i:PEP%i' % (i + 1, num + 1) for i, atom in enumerate(peptide)})
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
                    raise Exception('Maximum number of attempts to insert peptide %s reached!!!' % peptide.name)
            self.atoms.extend(model)
        logger.debug(module_name=_name, msg="Complex successfully created")

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
