"""
Classes Protein, Peptide, ProteinComplex - prepares initial complex.
"""

from copy import deepcopy
import json
from math import exp
import os
import re
from string import ascii_uppercase
from typing import Dict, List, Literal, Optional

from CABS.io import logger
from CABS.prediction import randinit
from CABS.structures.atom import Atoms
from CABS.structures.pdblib import Pdb
from CABS.structures.vector3d import Vector3d
from CABS.utils import utils

_name = "Protein"


class Protein(Atoms):
    """
    Class for the protein molecule.
    """

    NSP3_MODEL_PATH: str = ""

    def __init__(
        self,
        source: str,
        flexibility: Optional[Dict[str, float]] = None,
        exclude: Optional[List[str]] = None,
        weights: Optional[Dict[str, float]] = None,
        plddt: Optional[Dict[str, float]] = None,
        mode: Literal[
            "rigid", "flexible", "ss1", "no-protein-restraints", "unleashed", "none"
        ] = "rigid",
        work_dir: str = ".",
        receptor_ss: Optional[Dict[str, str]] = None,
        pdb_cache: Optional[str] = None,
        category: Optional[Dict[str, int]] = None,
        save_initial_pdb: bool = False,
        predict_peptide_structure: bool = False,
        cg2all_env_prefix: Optional[str] = None,
        sc: bool = False,
    ) -> None:
        Atoms.__init__(self)

        if source is None:
            self.source = None
            self.old_ids = {}
            self.new_ids = {}
            self.exclude = {}
            self.weights = []
            self.cg2all_env_prefix = cg2all_env_prefix
            self.center = Vector3d(0, 0, 0)
            self.dimension = 0.0
            self.patches = {}
            return

        logger.info(module_name=_name, msg=f"Loading {source} as input protein")

        # Only happens if user explicitly wants to predict peptide structure
        if predict_peptide_structure:
            try:
                self.atoms = randinit.RandomInitialStructure(source).atoms
            except Exception as e:
                logger.exit_program(
                    module_name=_name,
                    msg=f"Invalid input {source} for peptide structure prediction",
                    exc=e,
                )
            predictor = None
            if ":" not in source:
                if self.NSP3_MODEL_PATH:
                    try:
                        from CABS.prediction.secstrpredictor import SecStrPredictor

                        predictor = SecStrPredictor(self.NSP3_MODEL_PATH)
                    except ImportError as e:
                        logger.warning(
                            module_name=_name,
                            msg=f"NetSurfP-3.0 library or its dependencies are missing: {e!s}",
                        )
                    except Exception as e:
                        logger.warning(
                            module_name=_name,
                            msg=f"Cannot load NetSurfP-3.0 model: {e!s}",
                        )
                else:
                    logger.warning(
                        module_name=_name, msg="NSP3 model path not provided."
                    )

                if not predictor:
                    logger.warning(
                        module_name=_name,
                        msg="Secondary structure prediction will not be performed.",
                    )

            if predictor:
                logger.info(
                    module_name=_name,
                    msg="Running secondary structure prediction for the peptide using NetSurfP-3.0.",
                )
                try:
                    sec_str = predictor.predict_q3(sequence_to_predict=source)
                    ss = dict(
                        (a.resid_id(), sec_str[i]) for i, a in enumerate(self.atoms)
                    )
                    logger.info(
                        module_name=_name,
                        msg="Secondary structure prediction for the peptide successful.",
                    )
                except Exception as e:
                    logger.warning(
                        module_name=_name,
                        msg=f"Secondary structure prediction for the peptide failed: {e!s}",
                    )
                    CABS_SS = "CHTE"
                    ss = dict(
                        (a.resid_id(), CABS_SS[int(a.occ) - 1]) for a in self.atoms
                    )
            else:
                CABS_SS = "CHTE"
                ss = dict((a.resid_id(), CABS_SS[int(a.occ) - 1]) for a in self.atoms)

        else:
            try:
                self.atoms = randinit.RandomInitialStructure(source).pdb
                CABS_SS = "CHTE"
                ss = dict((a.resid_id(), CABS_SS[int(a.occ) - 1]) for a in self.atoms)

            except Exception as e:
                logger.debug(module_name=_name, msg=f"RandomInitialStructure failed or bypassed. Loading via Pdb class. Error: {e}")
                pdb = Pdb(source=source, pdb_cache=pdb_cache)
                ss = pdb.dssp(work_dir=work_dir)
                logger.debug(module_name=_name, msg=f"save_initial_pdb flag: {save_initial_pdb}")
                if save_initial_pdb:
                    logger.debug(module_name=_name, msg=f"Calling pdb.save_initial_pdb(work_dir={work_dir}) with sc={sc}")
                    if sc:
                        # Create a copy with SC for saving to start_all.pdb
                        sc_atoms = pdb.atoms.select("name CA")
                        sc_atoms.add_side_chain_centers()
                        header = pdb.mk_ss_header(work_dir=work_dir)
                        path = os.path.join(work_dir, "output_pdbs", "start_all.pdb")
                        os.makedirs(os.path.dirname(path), exist_ok=True)
                        sc_atoms.save_to_pdb(path, header=header)
                    else:
                        pdb.save_initial_pdb(work_dir=work_dir)
                pdb.atoms = pdb.atoms.select("name CA")
                if not pdb.atoms:
                    raise Exception(
                        "No protein alpha carbon (CA) atoms were found in the input structure. "
                        "Please make sure the input PDB/mmCIF contains a valid protein structure."
                    )
                self.atoms = pdb.atoms.models()[0]

        if receptor_ss:
            logger.info("Running manual assignment of receptor's II structure.")
            try:
                ss = ReceptorSS(current_ss=ss, receptor_ss=receptor_ss).ss
            except InvalidReceptorSS:
                logger.warning(msg="Invalid data for --receptor-ss option")

        # setup plddt
        if plddt:
            if plddt.lower() == "pdb" or plddt.lower() == "bf":
                for a in self.atoms:
                    a.plddt = a.bfac / 100
            elif plddt.lower() == "file":
                try:
                    d, de = self.read_plddt(plddt)
                    self.update_plddt(d, de)
                except Exception as e:
                    try:
                        protein_work_dir = os.path.dirname(source)
                        d, de = self.read_plddt(
                            os.path.join(protein_work_dir, "plddt.config")
                        )
                        self.update_plddt(d, de)
                    except Exception as e:
                        logger.warning(
                            _name, "Using default plddt(1.0) for all residues."
                        )
                        self.set_plddt(1.0)
            else:
                try:
                    d, de = self.read_plddt(plddt)
                    self.update_plddt(d, de)
                except OSError:
                    logger.warning(_name, f"Could not read pLDDT file: {plddt}")
                    logger.warning(_name, "Using default plddt(1.0) for all residues.")
                    self.set_plddt(1.0)
                except Exception as e:
                    logger.warning(_name, f"{e}")
                    logger.warning(_name, "Using default plddt(1.0) for all residues.")
                    self.set_plddt(1.0)
        else:
            self.set_plddt(1.0)

        # setup flexibility
        if flexibility:
            try:
                flex = float(flexibility)
                self.set_flexibility(flex)
            except ValueError:
                if flexibility.lower() == "bf":
                    for a in self.atoms:
                        a.flexibility = a.bfac
                elif flexibility.lower() == "bfi":
                    for a in self.atoms:
                        if a.bfac > 1.0:
                            a.flexibility = 0.0
                        elif a.bfac < 0.0:
                            a.flexibility = 1.0
                        else:
                            a.flexibility = 1.0 - a.flexibility
                elif flexibility.lower() == "bfg":
                    for a in self.atoms:
                        if a.bfac < 0.0:
                            a.flexibility = 1.0
                        else:
                            a.flexibility = exp(-0.5 * a.bfac**2)
                else:
                    try:
                        d, de = self.read_flexibility(flexibility)
                        self.update_flexibility(d, de)
                    except OSError:
                        logger.warning(
                            _name, f"Could not read flexibility file: {flexibility}"
                        )
                        logger.warning(
                            _name, "Using default flexibility(1.0) for all residues."
                        )
                        self.set_flexibility(1.0)
        else:
            self.set_flexibility(1.0)

        # setup excluding
        self.exclude = {}
        if exclude:
            for s in exclude:
                words = s.split("@")
                if len(words) == 1:
                    key = "ALL"
                else:
                    key = utils.pep2pep1(words[-1])
                if key in self.exclude:
                    self.exclude[key] += "+" + words[0]
                else:
                    self.exclude[key] = words[0]

            for k, v in self.exclude.items():
                self.exclude[k] = []
                for word in v.split("+"):
                    if ":" in word:
                        if "-" in word:
                            beg, end = word.split("-")
                            self.exclude[k].extend(self.atom_range(beg, end))
                        else:
                            self.exclude[k].append(word)
                    else:
                        chains = re.sub(r"[^%s]*" % word, "", ascii_uppercase)
                        self.exclude[k].extend(
                            a.resid_id() for a in self.select("chain %s" % chains)
                        )

        self.update_sec(ss)

        if work_dir and logger.output_ss():
            output_ss = os.path.join(work_dir, "output_data", "ss.txt")
            odir = os.path.dirname(output_ss)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            logger.to_file(
                filename=output_ss,
                content=str("".join(ss.values())),
                msg=f"Saving secondary structure output to {output_ss}",
            )

        # setup categories
        if category:
            try:
                d, de = self.read_category(category)
                self.update_category(d, de)
            except OSError:
                logger.warning(_name, f"Could not read category file: {category}")
                logger.warning(
                    _name,
                    "Using default categories based on pLDDT and SS for all residues.",
                )
                self.determine_category(mode=mode)
        else:
            self.determine_category(mode=mode)

        self.old_ids = self.fix_broken_chains()

        todrop = [c for c, l in self.list_chains().items() if l < 4]
        if todrop:
            self.atoms = self.atoms.drop("chain " + ",".join(todrop))

        self.new_ids = {v: k for k, v in self.old_ids.items()}

        for key, val in self.exclude.items():
            self.exclude[key] = [self.new_ids[r] for r in val]

        self.cg2all_env_prefix = cg2all_env_prefix
        # setup rmsd_weights
        self.weights = None
        if weights and weights.lower() == "flex":
            self.weights = [a.flexibility for a in self.atoms]
        elif weights and weights.lower() == "ss":
            self.weights = [(a.occ + 1.0) % 2 for a in self.atoms]
        elif weights and (weights.lower() == "off" or weights.lower() == "gauss"):
            self.weights = [1.0] * len(self.atoms)
        elif weights:
            try:
                default = 1.0
                self.weights = []
                weights_dict = {}
                with open(weights) as _file:
                    for line in _file:
                        k, v = line.split()[:2]
                        weights_dict[k] = v

                if "default" in weights_dict:
                    default = float(weights_dict["default"])

                for a in self.atoms:
                    w = weights_dict.get(a.resid_id())
                    w = float(w) if w else default
                    self.weights.append(w)

            except (OSError, ValueError):
                logger.warning(_name, f"Could not read weights file: {weights}")
                logger.warning(_name, "Using default weights(1.0) for all atoms.")
                self.weights = [1.0] * len(self.atoms)
        else:
             self.weights = [1.0] * len(self.atoms)

        self.center = self.cent_of_mass()
        self.dimension = self.max_dimension()
        self.patches = {}

    def convert_patch(self, location):
        if location not in self.patches:
            chains = {}
            for res in [self.new_ids[r] for r in location.split("+")]:
                num, chid = res.split(":")
                if chid in chains:
                    chains[chid].append(num)
                else:
                    chains[chid] = [num]
            s = " or ".join(
                [
                    "(chain " + ch + " and resnum " + ",".join(chains[ch]) + ")"
                    for ch in chains
                ]
            )
            patch = self.select(s)
            self.patches[location] = (patch.cent_of_mass() - self.center).norm()
        return self.patches[location]

    @staticmethod
    def read_flexibility(filename):
        key = r"[0-9A-Z]+:[A-Z]"
        val = r"[0-9.]+"

        patt_range = re.compile(f"({key}) *-* *({key}) +({val})")
        patt_single = re.compile(f"({key}) +({val})")

        with open(filename) as f:
            d = {}
            def_val = 1.0
            for line in f:
                if re.search("default", line):
                    def_val = float(line.split()[-1])
                else:
                    match = re.search(patt_range, line)
                    if match:
                        n1, c1 = match.group(1).split(":")
                        n2, c2 = match.group(2).split(":")
                        n1 = int(n1)
                        n2 = int(n2)
                        if c1 != c2 or n1 > n2:
                            raise Exception(
                                f"Invalid range: '{line}' in file: {filename}!!!"
                            )
                        for i in range(n1, n2 + 1):
                            d[str(i) + ":" + c1] = float(match.group(3))
                    else:
                        match = re.search(patt_single, line)
                        if match:
                            d[match.group(1)] = float(match.group(2))
                        else:
                            raise Exception("Invalid syntax in flexibility file!!!")
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
            with open(file_path) as file:
                json_dict = json.load(file)
                if "plddt" not in json_dict:
                    raise ValueError(
                        "Validation failed: Missing 'plddt' field in JSON."
                    )
                plddt_values = json_dict["plddt"]
                if not isinstance(plddt_values, list):
                    raise ValueError("Validation failed: 'plddt' field must be a list.")
                if not all(
                    isinstance(value, (int, float)) and 0 <= value <= 100
                    for value in plddt_values
                ):
                    raise ValueError(
                        "Validation failed: All pLDDT scores must be numbers between 0 and 100."
                    )
            return "valid_json"
        except json.JSONDecodeError:
            # Not JSON, try TSV
            try:
                with open(file_path) as file:
                    # Skip the header line if present
                    next(file)
                    for line in file:
                        columns = line.strip().split("\t")
                        if len(columns) < 2:
                            raise ValueError(
                                "Validation failed: TSV file must have at least two columns."
                            )
                        float(columns[1])  # Ensure the second column is a valid number
                return "valid_tsv"
            except Exception as e:
                raise ValueError(
                    "Validation failed: File is neither valid JSON nor TSV."
                ) from e

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
                with open(filename) as file:
                    json_dict = json.load(file)
                    plddt_values = [float(entry) / 100 for entry in json_dict["plddt"]]
            elif file_type == "valid_tsv":
                with open(filename) as file:
                    next(file)
                    for line in file:
                        columns = line.strip().split("\t")
                        plddt_values.append(float(columns[1]) / 100)
        except Exception as e:
            raise ValueError("Error while reading pLDDT values.") from e

        for i, atom in enumerate(self.atoms):
            if i < len(plddt_values):
                d[atom.resid_id()] = plddt_values[i]
            else:
                raise ValueError(
                    "Mismatch between the number of residues and pLDDT values."
                )

        return d, def_val

    @staticmethod
    def read_category(filename):
        key = r"[0-9A-Z]+:[A-Z]"
        val = r"[0-9.]+"

        patt_range = re.compile(f"({key}) *-* *({key}) +({val})")
        patt_single = re.compile(f"({key}) +({val})")

        with open(filename) as f:
            d = {}
            def_val = None
            for line in f:
                if re.search("default", line):
                    def_val = float(line.split()[-1])
                else:
                    match = re.search(patt_range, line)
                    if match:
                        n1, c1 = match.group(1).split(":")
                        n2, c2 = match.group(2).split(":")
                        n1 = int(n1)
                        n2 = int(n2)
                        if c1 != c2 or n1 > n2:
                            raise Exception(
                                f"Invalid range: '{line}' in file: {filename}!!!"
                            )
                        for i in range(n1, n2 + 1):
                            d[str(i) + ":" + c1] = float(match.group(3))
                    else:
                        match = re.search(patt_single, line)
                        if match:
                            d[match.group(1)] = float(match.group(2))
                        else:
                            raise Exception("Invalid syntax in category file!!!")
            return d, def_val

    def generate_restraints(self, mode, gap, min_d, max_d):
        gap = int(gap)
        min_d = float(min_d)
        max_d = float(max_d)
        restr = []
        _len = len(self.atoms)

        if mode in ["manual", "plddt"]:
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
                            restr.append(f"{a1.resid_id()} {a2.resid_id()} {d} {w}")

        else:
            for i in range(_len):
                a1 = self.atoms[i]
                ssi = int(a1.occ) % 2
                if mode == "flexible" and ssi:
                    continue
                for j in range(i + gap, _len):
                    a2 = self.atoms[j]
                    ssj = int(a2.occ) % 2
                    if mode == "flexible" and ssj:
                        continue
                    if mode == "ss1" and (ssi * ssj):
                        continue
                    d = (a1.coord - a2.coord).length()
                    if min_d < d < max_d:
                        if a1.flexibility < a2.flexibility:
                            w = a1.flexibility
                        else:
                            w = a2.flexibility
                        if w:
                            restr.append(f"{a1.resid_id()} {a2.resid_id()} {d} {w}")
        return restr

    def generate_ca_restraints(self):
        """
        Generates restraints for the CA-CA bonds.
        :return: list of restraints
        """
        restraints = []
        for i in range(len(self.atoms) - 1):
            if self.atoms[i].chid == self.atoms[i + 1].chid:
                restraints.append(
                    f"{self.atoms[i].resid_id()} {self.atoms[i + 1].resid_id()} 3.8 1.0"
                )
        return restraints

    def calculate_distances(self):
        """
        Generate a matrix of distances between each C-alpha in the protein (server uses this)
        :return: NxN matrix of distances (dict of dicts {'1:A': {'1:A' :20, '2:A': 30}, ...} )
        """
        out = {}
        _len = len(self.atoms)
        for i in range(_len):
            a1 = self.atoms[i]
            out[a1.resid_id()] = {}
            for j in range(_len):
                a2 = self.atoms[j]
                d = (a1.coord - a2.coord).length()
                out[a1.resid_id()][a2.resid_id()] = f"{d:6.4f}"
        return out


class Peptide(Atoms):
    """
    Class for the peptides.
    """

    def __init__(self, source, conformation, location, work_dir=".", pdb_cache=None, cg2all_env_prefix=None, predict_peptide_structure=False):
        logger.info(
            module_name=_name,
            msg=f"Loading ligand: {source}, conformation - {conformation}, location - {location}",
        )
        try:
            # OPTIMIZATION: Heuristic to detect sequences vs PDBs/files
            # If not a file, not 4-char ID, and no special chars commonly used in IDs/files, treat as sequence.
            identifier = source.split(":")[0]
            if (not os.path.isfile(identifier) and 
                len(identifier) != 4 and 
                "_" not in identifier and 
                "." not in identifier):
                raise Pdb.InvalidPdbInput("Input looks like a sequence (no file/ID match)")

            pdb = Pdb(
                source=source, pdb_cache=pdb_cache, no_exit=True
            )
            ss = pdb.dssp(work_dir=work_dir)
            pdb.atoms = pdb.atoms.select("name CA")
            if not pdb.atoms:
                raise Pdb.InvalidPdbInput("No protein alpha carbon (CA) atoms found in structure.")
            atoms = pdb.atoms.models()[0]
            atoms.update_sec(ss)
        except Pdb.InvalidPdbInput:
            if ":" not in source and (predict_peptide_structure or Protein.NSP3_MODEL_PATH):
                # Try prediction for sequence without SS
                try:
                    from CABS.prediction.secstrpredictor import SecStrPredictor
                    # If model path is not set, SecStrPredictor will use propensity fallback
                    predictor = SecStrPredictor(Protein.NSP3_MODEL_PATH)
                    logger.info(_name, f"Predicting SS for peptide sequence: {source}")
                    sec_str = predictor.predict_q3(source)
                    atoms = Atoms(f"{source}:{sec_str}")
                except Exception as e:
                    logger.warning(_name, f"Peptide SS prediction failed: {e}. Defaulting to Coil.")
                    atoms = Atoms(source)
            else:
                atoms = Atoms(source)

        atoms.set_bfac(0.0)
        self.conformation = conformation
        self.location = location
        Atoms.__init__(self, atoms)
        self.cg2all_env_prefix = cg2all_env_prefix


class ProteinComplex(Atoms):
    """
    Class that assembles the initial complex.
    """

    def __init__(
        self,
        protein,
        flexibility,
        exclude,
        weights,
        plddt,
        category,
        mode,
        peptides,
        replicas,
        separation,
        insertion_attempts,
        insertion_clash,
        work_dir,
        receptor_ss,
        pdb_cache,
        save_initial_pdb,
        json_output=False,
        cg2all_env_prefix=None,
        predict_peptide_structure=False,
        sc=False,
    ):
        logger.debug(module_name=_name, msg="Preparing the complex")
        Atoms.__init__(self)

        self.protein = Protein(
            protein,
            flexibility=flexibility,
            exclude=exclude,
            weights=weights,
            plddt=plddt,
            category=category,
            mode=mode,
            work_dir=work_dir,
            receptor_ss=receptor_ss,
            pdb_cache=pdb_cache,
            save_initial_pdb=save_initial_pdb,
            predict_peptide_structure=predict_peptide_structure,
            cg2all_env_prefix=cg2all_env_prefix,
            sc=sc,
        )
        self.chain_list = self.protein.list_chains()
        self.protein_chains = list(self.chain_list.keys())
        self.weights = list(self.protein.weights) if self.protein.weights is not None else []
        self.old_ids = deepcopy(self.protein.old_ids)

        self.peptides = []
        self.peptide_chains = []
        if peptides:
            taken_chains = "".join(self.protein_chains) + "X"
            for num, p in enumerate(peptides):
                peptide = Peptide(*p, work_dir=work_dir, pdb_cache=pdb_cache,
                                  cg2all_env_prefix=cg2all_env_prefix,
                                  predict_peptide_structure=predict_peptide_structure)
                if peptide[0].chid in taken_chains:
                    peptide.change_chid(
                        peptide[0].chid, utils.next_letter(taken_chains)
                    )
                taken_chains += peptide[0].chid
                self.peptide_chains.append(peptide[0].chid)
                self.peptides.append(peptide)
                self.weights.extend([1.0] * len(peptide))
                update_dict = {}
                i = 1
                for atom in peptide:
                    update_dict[atom.resid_id()] = f"{i}:PEP{num + 1}"
                    i += 1
                self.old_ids.update(update_dict)
                self.chain_list.update(peptide.list_chains())
        self.new_ids = {v: k for k, v in self.old_ids.items()}

        exclude = []
        for key, value in self.protein.exclude.items():
            if key == "ALL":
                kword = "PEP"
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
                    raise Exception(
                        f"Maximum number of attempts to insert peptide {peptide} reached!!!"
                    )
            self.atoms.extend(model)
        logger.debug(module_name=_name, msg="Complex successfully created")

        if json_output:
            complex_to_save = deepcopy(self)
            complex_to_save.update_ids(self.old_ids)
            json_file = os.path.join(work_dir, "output_data", "atoms.json")
            odir = os.path.dirname(json_file)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            complex_to_save.save_to_json(json_file)
            logger.debug(module_name=_name, msg="Atoms saved to JSON file")

    def update_ids(self, ids, pedantic=True):
        super().update_ids(ids, pedantic)
        all_chains = list(self.select("model 1").list_chains().keys())
        self.protein_chains = [c for c in all_chains if not c.startswith("PEP")]
        self.peptide_chains = [c for c in all_chains if c.startswith("PEP")]
        # If heuristic failed to find any peptides but there are chains, and we know we had peptides...
        if not self.peptide_chains and self.peptides and all_chains:
            # If there's no protein, everything must be peptides
            if not self.protein or not self.protein.atoms:
                self.peptide_chains = all_chains
            else:
                # This case is complex, but startswith('PEP') should have worked if we followed our own convention
                pass

    def map_user_chain(self, user_chain):
        """Maps user-friendly chain IDs (PEP1) to internal ones (A) during setup."""
        if user_chain.startswith("PEP"):
            try:
                idx = int(user_chain[3:]) - 1
                if 0 <= idx < len(self.peptides):
                    return self.peptides[idx][0].chid
            except (ValueError, IndexError):
                pass
        return user_chain

    def map_user_residue(self, user_residue):
        """Maps user-friendly residue IDs (e.g., '1:PEP1') to internal ones (e.g., '1:A')."""
        if ":" in user_residue:
            resnum, chid = user_residue.split(":", 1)
            return f"{resnum}:{self.map_user_chain(chid)}"
        return user_residue



    @staticmethod
    def insert_peptide(protein, peptide, separation):
        radius = 0.5 * protein.dimension + separation

        if peptide.location == "keep":
            location = peptide.cent_of_mass()
        elif peptide.location == "random":
            peptide.rotate_in_place(utils.random_rotation_matrix())
            location = Vector3d().random() * radius + protein.center
        else:
            location = protein.convert_patch(peptide.location) * radius + protein.center

        if peptide.conformation == "random":
            peptide.random_conformation()

        peptide.move_to(location)


class InvalidReceptorSS(Exception):
    pass


class ReceptorSS:
    def __init__(self, current_ss, receptor_ss):
        self.ss = {}
        chains = {}

        try:
            with open(receptor_ss) as f:
                for line in f:
                    chid, ss = line.replace(" ", "").split(":")
                    chains[chid] = ss
        except OSError:
            chains = dict(ch.split(":") for ch in receptor_ss.split("+"))

        if not current_ss:
            tmp = []
            for c in chains:
                s = list(chains[c])
                for l in range(1, len(s) + 1):
                    st = f"{l}:{c}"
                    tmp.append((st, s[l - 1]))
            current_ss = dict([(a[0], a[1]) for a in tmp])

        if not chains:
            raise InvalidReceptorSS

        current_chains = {}
        for k, v in current_ss.items():
            chid = k.split(":")[-1]
            if chid not in current_chains:
                current_chains[chid] = {}
            current_chains[chid][k] = v

        for chid in current_chains:
            if len(current_chains[chid]) == len(chains[chid]):
                self.ss.update(dict(zip(current_chains[chid].keys(), chains[chid])))
            else:
                raise InvalidReceptorSS


if __name__ == "__main__":
    ss = dict(
        [
            ("1:A", "C"),
            ("2:A", "H"),
            ("3:A", "H"),
            ("4:A", "H"),
            ("5:A", "C"),
            ("2:B", "C"),
            ("3:B", "E"),
            ("4:B", "E"),
            ("5:B", "C"),
            ("6:B", "H"),
            ("7:B", "H"),
            ("8:B", "C"),
        ]
    )
