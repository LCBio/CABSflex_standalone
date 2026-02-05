"""
Classes for representing atoms and collections of atoms from PDB structures.
"""

from collections import defaultdict
from copy import deepcopy
from itertools import combinations
import json
from math import sqrt
from random import randint
import re
from string import ascii_uppercase
from typing import (
    Any,
    Iterator,
    List,
    Literal,
    Optional,
    Union,
)

import numpy as np

from CABS.constants import (
    CABS_SS,
    CABS_SS_REVERSE,
)
from CABS.io.logger import ProgressBar
from CABS.structures.vector3d import Vector3d
from CABS.utils import utils
from CABS.utils.utils import (
    aa_to_long,
    aa_to_short,
    check_peptide_sequence,
    kabsch,
    smart_flatten,
)


class Atom:
    """
    Class for representation of a single atom.
    """

    # pattern used to decompose return value of resid_id() to (resnum, icode, chid)
    RES_ID_PATT = re.compile(r"(-?[0-9]{1,4})([^0-9]?):([A-Z0-9])")

    def __init__(
        self, line: Optional[str] = None, model: int = 0, **kwargs: Any
    ) -> None:
        """
        Constructor. Creates an Atom object from a PDB line or keyword arguments.
        """
        # 1. Initialize with basic defaults
        self.model = model
        self.hetatm = True
        self.serial = 0
        self.name = "XXXX"
        self.alt = ""
        self.resname = "XXX"
        self.chid = "X"
        self.resnum = 0
        self.icode = ""
        self.coord = Vector3d()
        self.occ = 0.0
        self.bfac = 0.0
        self.tail = ""

        # 2. If a raw PDB line is provided, parse it (Legacy support)
        if line:
            self.hetatm = line[:6] == "HETATM"
            self.serial = int(line[6:11])
            self.name = line[11:16].strip()
            self.alt = line[16]
            self.resname = line[17:21].strip()
            self.chid = line[21]
            self.resnum = int(line[22:26])
            self.icode = line[26]
            self.coord = Vector3d(line[30:38], line[38:46], line[46:54])
            self.occ = float(line[54:60])
            self.bfac = float(line[60:66])
            self.tail = " " * 11 + line[77:].replace("\n", "")

        # 3. OVERRIDE with kwargs (This is where data from the new pdblib.py arrives)
        for arg, value in kwargs.items():
            if hasattr(self, arg) or arg in ['plddt', 'flexibility', 'category', 'ss']:
                setattr(self, arg, value)

        # 4. CALCULATE derived metadata now that all data is present
        self.ss = getattr(self, 'ss', self.occ)
        self.flexibility = getattr(self, 'flexibility', self.bfac)
        self.plddt = getattr(self, 'plddt', self.bfac)
        self.category = getattr(self, 'category', 0.0)

        # Finally, set the category based on the final plddt/ss values
        self.set_category()

    def __str__(self) -> str:
        line = "ATOM  "
        if self.hetatm:
            line = "HETATM"
        fmt_name = f" {self.name:<3s}"
        if len(self.name) == 4:
            fmt_name = self.name
        line += f"{self.serial:5d} {fmt_name:4s}{self.alt:1s}{self.resname:<4s}{self.chid:1s}{self.resnum:4d}{self.icode:1s}   {self.coord:24s}{self.occ:6.2f}{self.bfac:6.2f} {self.tail}"
        return line

    def __repr__(self) -> str:
        return f"<Atom: {self.fmt()} {self.resname}>"

    def fmt(self) -> str:
        return f"{self.chid}{self.resnum}{self.icode.strip()}"

    def same_model(self, other: "Atom") -> bool:
        """
        Returns True if both atoms belong to the same model. False otherwise.
        :param other: Atom
        :return: Bool
        """
        return self.model == other.model

    def same_chain(self, other: "Atom") -> bool:
        """
        Returns True if both atoms belong to the same chain and model. False otherwise.
        :param other: Atom
        :return: Bool
        """
        return self.same_model(other) and self.chid == other.chid

    def same_residue(self, other: "Atom") -> bool:
        """
        Returns True if both atoms belong to the same residue, chain and model. False otherwise.
        :param other: Atom
        :return: Bool
        """
        return (
            self.same_chain(other)
            and self.resnum == other.resnum
            and self.icode == other.icode
        )

    def is_hydrogen(self) -> bool:
        """
        Returns true if Atom is hydrogen, false otherwise.
        Determined by the first non-digit character in atom's name. If "H" then hydrogen.
        :return: Bool
        """
        m = re.search("([A-Z])", self.name)
        return m and m.group(0) == "H"

    def dist2(self, other: "Atom") -> float:
        """
        Returns squared distance between two atoms.
        :param other: Atom
        :return: float
        """
        return (self.coord - other.coord).mod2()

    def distance(self, other: "Atom") -> float:
        """
        Returns distance in Angstroms between two atoms.
        :param other: Atom
        :return: float
        """
        return sqrt(self.dist2(other))

    def min_distance(self, other: "Atoms") -> float:
        """
        Returns minimal distance between atom and group of atoms.
        :param other: Atoms collection
        :return: float
        """
        return min(self.distance(atom) for atom in other)

    def match_token(self, token: str) -> bool:
        """
        Returns True if Atom matches selection token. False otherwise.
        :param token: Selection token string
        :return: Bool
        """
        words = token.split()
        if len(words) == 1:
            keyword = token.upper()
            if keyword == "HETERO":
                return self.hetatm
            else:
                raise Exception("Invalid selection keyword: " + keyword)
        elif len(words) > 1:
            keyword = words[0].upper()
            args = "".join(words[1:]).split(",")
            if keyword == "MODEL":
                return self.model in smart_flatten(args)
            elif keyword == "CHAIN":
                return self.chid in args
            elif keyword == "RESNUM":
                return self.resnum in smart_flatten(args)
            elif keyword == "RESNAME":
                return any([re.match("^%s$" % a, self.resname) for a in args])
            elif keyword == "NAME":
                return any([re.match("^%s$" % a, self.name) for a in args])
            else:
                raise Exception("Invalid selection keyword: " + keyword)
        else:
            raise Exception("Invalid selection syntax: " + token)

    def match(self, selection: "Selection") -> bool:
        """
         Returns True if Atom matches selection pattern. False otherwise.
        :param selection: Selection object
        :return: Bool
        """
        pattern = deepcopy(selection.tokens)
        for i, t in enumerate(pattern):
            if t.upper() not in Selection.JOINTS:
                pattern[i] = str(self.match_token(t))
        return eval(" ".join(pattern))

    def resid_id(self) -> str:
        """
        Returns a string with residue identification i.e. 123:A
        :return: str
        """
        return (str(self.resnum) + self.icode).strip() + ":" + self.chid

    def update_id(self, res_id: str) -> "Atom":
        """
        Updates resnum, chid and icode(when necessary) with values taken from dictionary res_id(i.e. 123A:B or 123:C).
        :param res_id: Residue ID string
        """

        match = re.match(Atom.RES_ID_PATT, res_id)
        if not match:
            raise Exception("Invalid res_id format: %s" % res_id)
        else:
            if match.group(2):
                self.icode = match.group(2)
            else:
                self.icode = " "
            self.resnum = int(match.group(1))
            self.chid = match.group(3)
        return self

    def set_category(
        self,
        mode: Literal[
            "rigid", "flexible", "no-protein-restraints", "unleashed", "none"
        ] = "rigid",
    ) -> None:
        category = 0

        if mode == "flexible":
            if self.occ == 2 or self.occ == 4:
                category = 3
        elif mode == "no-protein-restraints" or mode == "unleashed" or mode == "none":
            pass
        else:
            if self.plddt < 0.5:
                category -= 1
            elif self.plddt < 0.7:
                pass
            elif self.plddt < 0.9:
                category += 1
            else:
                category += 2

            if self.occ == 1:
                pass
            elif self.occ == 3:
                category += 1
            else:
                category += 2

            if category < 0:
                category = 0
            elif category > 3:
                category = 3

        self.category = category
        return self


class Atoms:
    """
    Container for atoms. Has most methods of a list. Also has methods common
    for all multi-atom objects: move, rotate etc.
    """

    def __init__(self, arg: Union[None, List[Atom], Any, str, int] = None) -> None:
        """
        Constructor. arg should be either:
            - None
            - list[Atom]
            - any object that has attribute named atoms, which is a list of objects named 'Atom'
            - string 'SEQUENCE' or 'SEQUENCE:SECONDARY' [SECONDARY is H/E/C/T for helix/sheet/turn/coil]
            - int for poly-alanine
        """
        if type(arg) is list:
            self.atoms = arg
        elif hasattr(arg, "atoms"):
            self.atoms = arg.atoms
        elif type(arg) is str:
            self.atoms = []
            if ":" in arg:
                seq, sec = arg.split(":")
                if len(sec) != len(seq):
                    raise Exception(
                        "Sequence length != secondary structure in " + arg + " !!!"
                    )
            else:
                seq = arg
                sec = "C" * len(seq)

            check_peptide_sequence(seq)

            for i, ch in enumerate(seq):
                self.atoms.append(
                    Atom(
                        hetatm=False,
                        serial=i + 1,
                        name="CA",
                        resname=aa_to_long(ch),
                        resnum=i + 1,
                        occ=CABS_SS.get(sec[i], 1),
                    )
                )

        elif type(arg) is int:
            self.atoms = []
            for i in range(arg):
                self.atoms.append(
                    Atom(
                        hetatm=False,
                        serial=i + 1,
                        name="CA",
                        resname="ALA",
                        resnum=i + 1,
                    )
                )

        else:
            self.atoms: List[Atom] = []

    def __len__(self) -> int:
        return len(self.atoms)

    def __iter__(self) -> Iterator[Atom]:
        return iter(self.atoms)

    def __getitem__(self, index: Union[int, slice]) -> Union[Atom, List[Atom]]:
        return self.atoms[index]

    def __setitem__(self, index: int, atom: Atom) -> None:
        self.atoms[index] = atom

    def __delitem__(self, index: int) -> None:
        del self.atoms[index]

    def append(self, atom: Atom) -> None:
        self.atoms.append(atom)

    def extend(self, other: "Atoms") -> None:
        self.atoms.extend(other.atoms)

    def __str__(self) -> str:
        return "\n".join(str(atom) for atom in self.atoms)

    def __repr__(self) -> str:
        return "<Atoms: %i>" % len(self.atoms)

    def __eq__(self, other: "Atoms") -> bool:
        if len(self.atoms) != len(other.atoms):
            return False
        for atom1, atom2 in zip(self.atoms, other.atoms):
            if atom1 != atom2:
                return False
        return True

    def __ne__(self, other: "Atoms") -> bool:
        return not (self == other)

    def residues(self) -> List["Atoms"]:
        """
        Returns a list of Atoms objects representing residues.
        :return: List of Atoms objects
        """
        res = []
        residue = Atoms()
        residue.append(self.atoms[0])
        for atom in self[1:]:
            if atom.same_residue(residue.atoms[-1]):
                residue.append(atom)
            else:
                if residue not in res:
                    res.append(residue)
                for r in res:
                    if atom.same_residue(r[0]):
                        r.append(atom)
                        break
                else:
                    residue = Atoms()
                    residue.append(atom)
        if residue not in res:
            res.append(residue)
        return res

    def chains(self):
        """
        Returns a list of Atoms objects representing chains.
        :return: [Atoms]
        """
        chn = []
        chain = Atoms()
        chain.append(self.atoms[0])
        for atom in self[1:]:
            if atom.same_chain(chain.atoms[-1]):
                chain.append(atom)
            else:
                if chain not in chn:
                    chn.append(chain)
                for ch in chn:
                    if atom.same_chain(ch[0]):
                        ch.append(atom)
                        break
                else:
                    chain = Atoms()
                    chain.append(atom)
        if chain not in chn:
            chn.append(chain)
        return chn

    def models(self):
        """
        Returns a list of Atoms objects representing models.
        :return: [Atoms]
        """
        mdl = []
        model = Atoms()
        model.append(self.atoms[0])
        for atom in self[1:]:
            if atom.same_model(model.atoms[-1]):
                model.append(atom)
            else:
                if model not in mdl:
                    mdl.append(model)
                for m in mdl:
                    if atom.same_model(m[0]):
                        m.append(atom)
                        break
                else:
                    model = Atoms()
                    model.append(atom)
        if model not in mdl:
            mdl.append(model)
        return mdl

    def to_numpy(self):
        """ "
        Returns np.array(N, 3) with coordinates, where N is the number of Atoms.
        """
        return np.array([a.coord.to_numpy() for a in self.atoms])

    def from_numpy(self, matrix):
        """
        Sets Atoms' coordinates from np.array(3,N) or (N,3).
        """
        if matrix.shape == (len(self), 3):
            for index, atom in enumerate(self.atoms):
                atom.coord = Vector3d(matrix[index])
        elif matrix.shape == (3, len(self)):
            self.from_numpy(matrix.T)
        else:
            raise Exception("Invalid matrix shape: " + str(matrix.shape))
        return self

    def move(self, v):
        """
        Move atoms by vector.
        :param v: Vector3d
        :return: Atoms
        """
        if v is not None:
            for a in self.atoms:
                a.coord += v
        return self

    def rotate(self, matrix):
        """
        Rotate atoms by rotation matrix.
        :param matrix: np.array(3, 3)
        :return: Atoms
        """
        if matrix.shape != (3, 3):
            raise Exception("Invalid matrix shape: " + matrix.shape)
        self.from_numpy(np.dot(matrix, self.to_numpy().T))
        return self

    def rotate_in_place(self, matrix):
        if matrix.shape != (3, 3):
            raise Exception("Invalid matrix shape: " + matrix.shape)
        com = self.cent_of_mass()
        self.move(-com)
        self.from_numpy(np.dot(matrix, self.to_numpy().T))
        self.move(com)
        return self

    def cent_of_mass(self):
        """
        Returns a vector of the geometrical center of atoms.
        :return: Vector3d
        """
        com = Vector3d()
        for atom in self.atoms:
            com += atom.coord
        com_x = com.x / len(self)
        com_y = com.y / len(self)
        com_z = com.z / len(self)

        return Vector3d(com_x, com_y, com_z)  # changed  to_numpy()  and it fucks

    def center_at_origin(self):
        """
        Moves atoms so that their geometrical center is in [0, 0, 0].
        :return: Atoms
        """
        self.move(-self.cent_of_mass())
        return self

    def move_to(self, v):
        """
        Moves atoms so that their geometrical center is in [vx, vy, vz].
        :param v: Vector3d
        :return: Atoms
        """
        self.move(v - self.cent_of_mass())
        return self

    def compute_rotation(self, other, concentric=False):
        """
        Computes the rotation matrix for best fit between two sets of Atoms.
        :param other: Atoms
        :param concentric: Bool
        :return: np.array(3, 3)
        """
        if len(self) != len(other):
            raise Exception(
                "Atom sets have different length: %i != %i" % (len(self), len(other))
            )
        t = other.to_numpy()
        q = self.to_numpy()
        return kabsch(t, q, concentric=concentric)

    def str_align(self, other):
        """
        Aligns structurally set of atoms to another set.
        :param other: Atoms
        :return: Atoms
        """
        r = self.compute_rotation(other)
        self.center_at_origin().rotate(r).move(other.cent_of_mass())
        return self

    def rmsd(self, other):
        """
        Calculates rmsd between two sets of atoms.
        :param other: Atoms
        :return: float
        """
        if len(self) != len(other):
            raise Exception(
                "Atom sets have different length: %i != %i" % (len(self), len(other))
            )
        r = 0
        for a1, a2 in zip(self, other):
            r += a1.dist2(a2)
        return sqrt(r / len(self))

    def min_distance(self, other):
        """
        Calculates minimal distance between two sets of atoms.
        :param other: Atoms
        :return: float
        """
        return min(a.min_distance(other) for a in self.atoms)

    def change_chid(self, old, new):
        """
        Changes chain ID. Can do multiple changes at once.
        :param old: str
        :param new: str
        :return: Atoms
        """
        if len(old) == len(new) and len(old) > 0:
            d = {}
            for o, n in zip(old, new):
                d[o] = n
            for atom in self:
                if atom.chid in d:
                    atom.chid = d[atom.chid]
        return self

    def model_count(self):
        """
        Returns number of models in Atoms object.
        :return: int
        """
        return len(self.models())

    def chain_count(self):
        """
        Returns number of chains in Atoms object.
        :return: int
        """
        return len(self.chains())

    def residue_count(self):
        """
        Returns number of residues in Atoms object.
        :return: int
        """
        return len(self.residues())

    def list_chains(self):
        """
        Returns a dictionary [chain ID] = chain_residue_count
        :return: {str: int}
        """
        d = {}
        for ch in self.chains():
            d[ch[0].chid] = Atoms(ch).residue_count()
        return d

    def max_dimension(self):
        """
        Returns maximal distance between any two atoms from the Atoms object.
        :return: float
        """
        return sqrt(max([p[0].dist2(p[1]) for p in combinations(self.atoms, 2)]))

    def make_pdb(self, bar_msg=""):
        """
        Returns a pdb-like formatted string. bar_msg is a string with message to show at ProgressBar initialization.
        bar_msg = '' disables the bar.
        :param bar_msg: str
        :return: str
        """
        models = self.models()
        if bar_msg:
            bar = ProgressBar(len(models), bar_msg)
        else:
            bar = None
        if len(models) == 1:
            s = str(self)
            s += "\n"
        else:
            s = ""
            for m in models:
                s += f"MODEL{m[0].model:9d}\n"
                s += str(m)
                s += "\nENDMDL\n"
                if bar:
                    bar.update()
        if bar:
            bar.done(False)
        return s

    def save_to_pdb(self, filename, bar_msg="", header=""):
        """
        Saves atoms to a file in the pdb format. Calls Atoms.make_pdb(). bar_msg is a string with message to show
        at ProgressBar initialization. bar_msg = '' disables the bar.
        :param filename: str
        :param bar_msg: str
        :param header: str
        :return: None
        """
        with open(filename, "w") as f:
            f.write(header)
            f.write(self.make_pdb(bar_msg=bar_msg))

    def save_to_json(self, filename):
        """
        Saves atoms to a file in the json format.
        :param filename:
        :return:
        """

        chain_indices = {
            chain_id: id_num
            for id_num, chain_id in enumerate(self.list_chains().keys())
        }

        chain_residues = defaultdict(dict)

        residue_counter = {}
        for atom in self.atoms:
            chain_id = atom.chid
            chain_num = chain_indices[chain_id]
            if chain_id not in residue_counter:
                residue_counter[chain_id] = 0

            if atom.resnum not in chain_residues[chain_num]:
                chain_residues[chain_num][atom.resnum] = residue_counter[chain_id]
                residue_counter[chain_id] += 1

        data = {chain_idx: res_dict for chain_idx, res_dict in chain_residues.items()}

        with open(filename, "w") as out:
            json.dump(data, out)

    def select(self, selection):
        """
        Selects subset of atoms defined by selection sentence.
        :param selection: str or Selection
        :return: Atoms
        """
        if type(selection) is str:
            s = Selection(selection)
        else:
            s = selection

        return Atoms([a for a in self if a.match(s)])

    def drop(self, selection):
        """
        Removes subset of atoms defined by selection sentence.
        :param selection: str or Selection
        :return: Atoms
        """
        if type(selection) is str:
            s = Selection(selection)
        else:
            s = selection
        return self.select(~s)

    def get_resname(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = resname.
        :return: {str: str}
        """
        res = {}
        for a in self.atoms:
            res[a.resid_id()] = a.resname
        return res

    def get_resname_short(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = resname.
        :return: {str: str}
        """
        res = {}
        for a in self.atoms:
            res[a.resid_id()] = aa_to_short(a.resname)
        return res

    def get_coordinates(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = coordinates.
        :return: {str: Vector3d}
        """
        res = {}
        for a in self.atoms:
            res[a.resid_id()] = a.coord
        return res

    def update_sec(self, sec):
        """
        Reads secondary structure dictionary sec[] with Atoms.resid_id() as keys
        and puts it into Atom.occ in CABS code:
        Helix - > 2.0, Sheet -> 4.0, Turn -> 3.0, Coil -> 1.0
        :param sec: {str: str}
        :return: Atoms
        """
        if sec:
            for a in self.atoms:
                a.occ = CABS_SS[sec.get(a.resid_id(), "C")]
        return self

    def get_sec(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = secondary structure.
        :return: {str: str}
        """
        sec = {}
        for a in self.atoms:
            sec[a.resid_id()] = CABS_SS_REVERSE[a.occ]
        return sec

    def update_occ(self, occ):
        """
        Reads dictionary with keys = Atom.resid_id() and values = occupancy and puts it into Atom.occ
        if key is found, default otherwise.
        :param occ: {str: float}
        :return: Atoms
        """
        for a in self.atoms:
            a.occ = occ.get(a.resid_id(), 1.0)
        return self

    def get_occ(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = beta factors.
        :return: {str: float}
        """
        sec = {}
        for a in self.atoms:
            sec[a.resid_id()] = a.occ
        return sec

    def update_bfac(self, bfac, default=0.0):
        """
        Reads dictionary with keys = Atom.resid_id() and values = beta factors and puts it into Atom.bfac
        if key is found, default otherwise.
        :param bfac: {str: float}
        :param default: float
        :return: Atoms
        """
        for a in self.atoms:
            a.bfac = bfac.get(a.resid_id(), default)
        return self

    def set_bfac(self, bfac=0.0):
        """
        Sets beta factor of all atoms to bfac.
        :param bfac: float
        :return: Atoms
        """
        for atom in self.atoms:
            atom.bfac = bfac
        return self

    def get_bfac(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = beta factors.
        :return: {str: float}
        """
        bfac = {}
        for a in self.atoms:
            bfac[a.resid_id()] = a.bfac
        return bfac

    def update_flexibility(self, flexibility, default=1.0):
        """
        Reads dictionary with keys = Atom.resid_id() and values = flexibility and puts it into Atom.flexibility
        if key is found, default otherwise.
        :param flexibility: {str: float}
        :param default: float
        :return: Atoms
        """
        for a in self.atoms:
            a.flexibility = flexibility.get(a.resid_id(), default)
        return self

    def set_flexibility(self, flexibility=1.0):
        """
        Sets beta factor of all atoms to bfac.
        :param flexibility: float
        :return: Atoms
        """
        for atom in self.atoms:
            atom.flexibility = flexibility
        return self

    def update_plddt(self, plddt, default=1.0):
        """
        Reads dictionary with keys = Atom.resid_id() and values = plddt and puts it into Atom.plddt
        if key is found, default otherwise.
        :param plddt: {str: float}
        :param default: float
        :return: Atoms
        """
        for a in self.atoms:
            a.plddt = plddt.get(a.resid_id(), default)
        return self

    def set_plddt(self, plddt=1.0):
        """
        Sets pLDDT of all atoms to plddt.
        :param plddt: float
        :return: Atoms
        """
        for atom in self.atoms:
            atom.plddt = plddt
        return self

    def get_plddt(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = plddt.
        :return: {str: float}
        """
        plddt = {}
        for a in self.atoms:
            plddt[a.resid_id()] = a.plddt
        return plddt

    def update_category(self, category, default=None):
        """
        Reads dictionary with keys = Atom.resid_id() and values = beta factors and puts it into Atom.category
        if key is found, default otherwise.
        :param category: {str: float}
        :param default: float
        :return: Atoms
        """
        if default:
            for a in self.atoms:
                a.category = category.get(a.resid_id(), default)
        else:
            for a in self.atoms:
                if a.resid_id() in category.keys():
                    a.category = category[a.resid_id()]
                else:
                    a.set_category()
        return self

    def determine_category(self, mode="rigid"):
        """
        Sets category of all atoms according to plddt and secondary structure.
        :return: Atoms
        """
        for atom in self.atoms:
            atom.set_category(mode=mode)
        return self

    def get_category(self):
        """
        Returns dictionary with keys = Atom.resid_id() and values = category.
        :return: {str: float}
        """
        category = {}
        for a in self.atoms:
            category[a.resid_id()] = a.category
        return category

    def valid_residues(self, must_have="CA, N, C, O"):
        """
        Returns only those residues that have atoms specified in "must_have" parameter.
        TODO: This is just temporary and it will be replaced by conditional selection class.
        :param must_have: str
        :return: Atoms
        """
        valid = Atoms()
        mh = [word.strip() for word in must_have.split(",")]
        for residue in self.residues():
            keep = True
            for nm in mh:
                if len(residue.select("name " + nm)) == 0:
                    keep = False
                    break
            if keep:
                valid += residue
        return valid

    def remove_alternative_locations(self):
        """
        Removes atoms with alternative locations other than ' ' or 'A'
        :return: Atoms
        """
        self.atoms = [
            atom for atom in self.atoms if (atom.alt == " " or atom.alt == "A")
        ]
        return self

    def set_model_number(self, number):
        """
        Sets model number to [number] for all atoms.
        :param number: int
        :return: Atoms
        """
        for a in self.atoms:
            a.model = number
        return self

    def fix_broken_chains(self, cut_off=4.5, used_letters=""):
        """
        Checks for gaps in protein chains (Ca-Ca distance > cut_off). Splits broken chains
        on gaps taking next available letter for the new chain, except for those in used_letters.
        Returns a dictionary with residue ids (new -> old).
        :param cut_off: float
        :param used_letters: str
        :return: Atoms
        """

        used_letters += "".join(self.list_chains().keys())

        old_ids = {}
        for chain in self.chains():
            prev = None
            for residue in chain.residues():
                ca = residue.select("name CA")[0]
                res_id = ca.resid_id()
                if prev:
                    chid = prev.chid
                    d = (ca.coord - prev.coord).length()
                    if d > cut_off:
                        chid = sorted(
                            re.sub("[" + used_letters + "]", "", ascii_uppercase)
                        )[0]
                        used_letters += chid
                    for a in residue:
                        a.chid = chid
                old_ids[ca.resid_id()] = res_id
                prev = ca
        return old_ids

    def update_ids(self, ids, pedantic=True):
        """
        Updates resnum, icode, chid from dictionary ids with pairs old, new resid_id. Pedantic controls behaviour
        if key not found in ids. pedantic = True raises Exception, pedantic = False does nothing
        :param ids: {str: str}
        :param pedantic: Bool
        :return: Atoms
        """
        for a in self.atoms:
            r_id = ids.get(a.resid_id())
            if not r_id:
                if pedantic:
                    raise Exception(
                        "%s not found in %s" % (a.resid_id(), sorted(ids.keys()))
                    )
            else:
                a.update_id(r_id)
        return self

    def erase_tail(self):
        """
        Clears atoms' tail field.
        :return: Atoms
        """
        for a in self.atoms:
            a.tail = ""
        return self

    def atom_range(self, first, last):
        atoms = []
        add = False
        for a in self.atoms:
            ar = a.resid_id()
            if ar == first:
                add = True
            if add:
                atoms.append(ar)
            if ar == last:
                return atoms
        return []

    def random_conformation(self, lib=utils.RANDOM_LIGAND_LIBRARY):
        length = len(self)
        models, max_length, dim = lib.shape
        if length > max_length:
            raise Exception(
                "Cannot generate random coordinates for peptide length = %d (max is %d)"
                % (length, max_length)
            )
        model = randint(0, models - 1)
        index = randint(0, max_length - length)
        self.from_numpy(lib[model][index : index + length])
        return self


class Selection:
    """
    Class representing atom selection.
    TODO: aliases
    """

    ARGS_KEYWORDS = ["MODEL", "CHAIN", "RESNUM", "RESNAME", "NAME"]
    NO_ARGS_KEYWORDS = ["HETERO"]
    OPERATORS = ["NOT", "AND", "OR"]
    PARENTHESIS = ["(", ")"]
    JOINTS = OPERATORS + PARENTHESIS
    KEYWORDS = ARGS_KEYWORDS + NO_ARGS_KEYWORDS + JOINTS

    def __init__(self, s=""):
        """
        Takes a string as input and parses it into selection tokens
        :param s: str
        :return: [str]
        """
        self.tokens = []
        if s is not None:
            for word in s.replace("(", " ( ").replace(")", " ) ").split():
                if word.upper() in self.KEYWORDS:
                    self.tokens.append(word)
                elif len(self.tokens) != 0:  # changed to != 0
                    self.tokens[-1] += " " + word

    def __invert__(self):
        return Selection("not( " + repr(self) + " )")

    def __repr__(self):
        return " ".join(self.tokens)


if __name__ == "__main__":
    pass
