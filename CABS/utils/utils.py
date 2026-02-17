"""Utility functions for CABS with type annotations."""

from contextlib import closing
import os
from pathlib import Path
import re
from string import ascii_uppercase
import subprocess
from tempfile import NamedTemporaryFile
from typing import (
    Any,
    List,
    Literal,
    Optional,
    TextIO,
    Tuple,
    Union,
)
import warnings

import numpy as np
import numpy.typing as npt

# Set numpy to handle denormals more gracefully - this is the proper way
# to handle IEEE_DENORMAL warnings without masking legitimate small values
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # This tells numpy to flush denormals to zero automatically in operations
    np.seterr(under="ignore")  # Ignore underflow warnings (denormals)

try:
    from biopandas.pdb import PandasPdb
    import pandas as pd

    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    from importlib.resources import as_file, files
except ImportError:
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=DeprecationWarning, module="pkg_resources"
        )
        from pkg_resources import resource_filename

from CABS.constants import (
    AA_NAMES,
    AA_SUB_NAMES,
    AA_SUB_NAMES_EXTENDED,
    CABS_SS,
    CABS_SS_REVERSE,
    SIDECNT,
    AminoAcid,
    AminoAcidCode,
    SecondaryStructure,
)
from CABS.io import logger

__all__ = [
    "AA_NAMES",
    "AA_SUB_NAMES",
    "CABS_SS",
    "CABS_SS_REVERSE",
    "RANDOM_LIGAND_LIBRARY",
    "SIDECNT",
    "AminoAcid",
    "EprGenerator",
    "InvalidAAName",
    "SCModeler",
    "SecondaryStructure",
    "aa_to_long",
    "aa_to_short",
    "check_peptide_sequence",
    "convert_cg_to_all",
    "dynamic_kabsch",
    "fix_residue",
    "kabsch",
    "line_count",
    "next_letter",
    "random_rotation_matrix",
    "ranges",
    "rmsd",
    "smart_flatten",
]

# Type aliases
RotationMatrix = npt.NDArray[np.float64]
CoordinateArray = npt.NDArray[np.float64]  # Shape (N, 3)
TrajectoryArray = npt.NDArray[np.float64]  # Shape (frames, atoms, 3)
WeightArray = Optional[npt.NDArray[np.float64]]


# Load random ligand library
def _load_random_ligand_library() -> npt.NDArray[np.float64]:
    """Load the random ligand library from data file."""
    try:
        # Try modern importlib.resources first
        try:
            with as_file(files("CABS") / "data" / "data2.dat") as data_file:
                data = np.loadtxt(str(data_file))
        except (ImportError, AttributeError, NameError):
            # Fallback to pkg_resources
            data_path = resource_filename("CABS", "data/data2.dat")
            data = np.loadtxt(data_path)
    except Exception:
        # Return empty array if data file cannot be loaded
        return np.array([], dtype=float).reshape(-1, 3)

    return data.reshape(1, -1, 3)


RANDOM_LIGAND_LIBRARY: npt.NDArray[np.float64] = _load_random_ligand_library()


class SCModeler:
    """Side Chain Modeler.

    Rebuilds side chains in CABS representation on C-alpha trace vector.
    """

    def __init__(self, nms: Any) -> None:
        """Side Chain Modeler initialization.

        Arguments:
            nms: sequence of CABS.Atom representing subsequent mers.
        """
        self.nms = nms

    @staticmethod
    def _mk_local_system(
        c1: npt.NDArray[np.float64],
        c2: npt.NDArray[np.float64],
        c3: npt.NDArray[np.float64],
    ) -> Tuple[
        npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64], float
    ]:
        """Return local system base for given CA.

        Arguments:
            c1, c2, c3: subsequent CA position vectors.

        Base will be calculated for c2.
        Returns 3 vectors and distance between c1 and c3.
        """
        rdif = c3 - c1
        rdnorm = np.linalg.norm(rdif)
        rsum = (c3 - c2) + (c1 - c2)
        rsum_norm = np.linalg.norm(rsum)
        if rsum_norm < 1e-12:  # Avoid division by very small numbers
            z = np.array([0.0, 0.0, 1.0])  # Default z direction
        else:
            z = -1 * rsum / rsum_norm
        x = rdif / rdnorm
        y = np.cross(z, x)
        return x, y, z, rdnorm

    @staticmethod
    def _calc_nodes_line(
        old_v: npt.NDArray[np.float64], new_v: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Return versor unchanged during transformation of one versor into another."""
        w = np.cross(old_v, new_v)
        norm_w = np.linalg.norm(w)
        if norm_w < 1e-12:  # Avoid division by very small numbers
            return np.array([1.0, 0.0, 0.0])  # Return default unit vector
        return w / norm_w

    @staticmethod
    def _calc_trig_fnc(
        ve1: npt.NDArray[np.float64],
        ve2: npt.NDArray[np.float64],
        axis: npt.NDArray[np.float64],
    ) -> Tuple[float, float]:
        """Return cos and sin between given versors."""
        cos = np.dot(ve1, ve2)
        perp_vec = np.cross(ve1, ve2)
        sin = np.linalg.norm(perp_vec) * np.sign(np.dot(perp_vec, axis))
        return cos, sin

    def _calc_rot_mtx(
        self,
        c1: npt.NDArray[np.float64],
        c2: npt.NDArray[np.float64],
        c3: npt.NDArray[np.float64],
        dbg: bool = False,
    ) -> Tuple[npt.NDArray[np.float64], float]:
        """Return rotation matrix transforming Cartesian system to system
        of alpha carbon c2 in sequence of alpha carbons c1, c2, c3.

        Arguments:
            c1, c2, c3: position vectors of subsequent alpha carbons.

        Returns matrix and distance between c1 and c3.
        """
        if dbg:
            # Debug functionality removed for production code
            pass
        x, y, z, rdnorm = self._mk_local_system(c1, c2, c3)

        setting = np.geterr()
        np.seterr(all="raise")
        try:
            w = self._calc_nodes_line(np.array((0, 0, 1)), z)
        except FloatingPointError:
            w = x
        np.seterr(**setting)

        cph, sph = self._calc_trig_fnc(np.array((1, 0, 0)), w, np.array((0, 0, 1)))
        # phi angle trig fncs -- rotation around z axis so x -> w

        cps, sps = self._calc_trig_fnc(w, x, z)
        # psi angle -- rotation around z' so x -> x'

        cth, sth = self._calc_trig_fnc(np.array((0, 0, 1)), z, w)
        # theta angle -- rotation around nodes line to transform z on z'

        rot = np.matrix(
            [
                [cps * cph - sps * sph * cth, sph * cps + sps * cth * cph, sps * sth],
                [
                    -1 * sps * cph - sph * cps * cth,
                    -1 * sps * sph + cps * cth * cph,
                    cps * sth,
                ],
                [sth * sph, -sth * cph, cth],
            ]
        )

        return rot, rdnorm

    @staticmethod
    def _calc_scatter_coef(dist: float) -> float:
        """Calculate scatter coefficient based on distance."""
        if dist < 5.3:
            return 1.0
        if dist > 6.4:
            return 0.0
        return float((dist - 5.3) * -(1 / 1.1) + 1)

    def rebuild_one(
        self, vec: npt.NDArray[np.float64], sc: bool = True
    ) -> npt.NDArray[np.float64]:
        """Takes vector of C alpha coords and residue names and returns vector of C beta coords."""
        vec = np.insert(vec, 0, vec[0] - (vec[2] - vec[1]), axis=0)
        vec = np.append(vec, np.array([vec[-1] + (vec[-2] - vec[-3])]), axis=0)

        nvec = np.zeros((len(vec) - 2, 3))
        nms = (lambda x: self.nms[x].resname) if sc else (lambda x: "ALA")

        for i in range(len(vec) - 2):
            rot, casdist = self._calc_rot_mtx(*vec[i : i + 3])
            coef = self._calc_scatter_coef(casdist)
            comp = np.array(SIDECNT[nms(i)][:3]) * coef
            scat = np.array(SIDECNT[nms(i)][3:]) * (1 - coef)
            rbld = np.array(comp + scat)
            nvec[i] = rbld.dot(rot).A1 + vec[i + 1]

        return nvec

    def _calculate_traj(
        self, traj: TrajectoryArray, sc: bool = False
    ) -> TrajectoryArray:
        """Calculate trajectory with side chain modeling."""
        return np.array([np.array([self.rebuild_one(j, sc) for j in i]) for i in traj])

    def calculate_cb_traj(self, traj: TrajectoryArray) -> TrajectoryArray:
        """Calculate C-beta trajectory."""
        return self._calculate_traj(traj, False)

    def calculate_sc_traj(self, traj: TrajectoryArray) -> TrajectoryArray:
        """Calculate side chain trajectory."""
        return self._calculate_traj(traj, True)


class InvalidAAName(Exception):
    """Exception raised when invalid amino acid name is used"""

    def __init__(self, name: str, length: int) -> None:
        self.name = (name, length)

    def __str__(self) -> str:
        return f"{self.name[0]} is not a valid {self.name[1]}-letter amino acid code"


def aa_to_long(seq: str) -> str:
    """Converts short amino acid name to long."""
    s = seq.upper()
    if s in AA_NAMES:
        return AA_NAMES[s]
    else:
        raise InvalidAAName(seq, 1)


def aa_to_short(seq: str) -> AminoAcidCode:
    """Converts long amino acid name to short."""
    s = seq.upper()
    for short, full in AA_NAMES.items():
        if full == s:
            return short
    else:
        raise InvalidAAName(seq, 3)


def next_letter(taken_letters: str) -> str:
    """Returns next available letter for new protein chain."""
    return re.sub("[" + taken_letters + "]", "", ascii_uppercase)[0]


def line_count(filename: Union[str, Path]) -> int:
    """Count lines in a file."""
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f, 1):
            pass
    return i


def ranges(data: List[int]) -> List[Tuple[int, int]]:
    """Convert list of integers to ranges."""
    result = []
    if not data:
        return result
    idata = iter(data)
    first = prev = next(idata)
    for following in idata:
        if following - prev == 1:
            prev = following
        else:
            result.append((first, prev + 1))
            first = prev = following
    result.append((first, prev + 1))
    return result


def kabsch(
    target: CoordinateArray,
    query: CoordinateArray,
    weights: WeightArray = None,
    concentric: bool = False,
) -> RotationMatrix:
    """
    Function for the calculation of the best fit rotation.

    Arguments:
        target: a N x 3 np.array with coordinates of the reference structure
        query: a N x 3 np.array with coordinates of the fitted structure
        weights: a N-length list with weights - floats from [0:1]
        concentric: True/False specifying if target and query are centered at origin

    IMPORTANT: If weights are not None centering at origin should account for them.
    proper procedure: A -= np.average(A, 0, WEIGHTS)

    Returns rotation matrix as 3 x 3 np.array
    """

    if not concentric:
        t = target - np.average(target, axis=0, weights=weights)
        q = query - np.average(query, axis=0, weights=weights)
    else:
        t = target
        q = query

    c = np.dot(weights * t.T, q) if weights is not None else np.dot(t.T, q)
    v, s, w = np.linalg.svd(c)
    d = np.identity(3)
    if np.linalg.det(c) < 0:
        d[2, 2] = -1

    return np.dot(np.dot(w.T, d), v.T)


_LARGE = 1000.0
_TINY = 0.001


def rmsd(
    target: CoordinateArray,
    query: Optional[CoordinateArray] = None,
    weights: WeightArray = None,
) -> float:
    """Calculate RMSD between target and query coordinates."""
    _diff = target if query is None else query - target
    _rmsd_sq = np.average(np.sum(_diff**2, axis=1), weights=weights)
    _rmsd = np.sqrt(_rmsd_sq)
    return _rmsd if _rmsd > _TINY else 0.0


GAUSS_MAX_ITER = 100


def dynamic_kabsch(
    target: CoordinateArray, query: CoordinateArray
) -> Tuple[float, RotationMatrix, npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Dynamic Kabsch algorithm with iterative weighting."""
    _rmsd = _LARGE
    w = [1.0] * len(target)
    for i in range(GAUSS_MAX_ITER):
        t_com = np.average(target, 0, w)
        q_com = np.average(query, 0, w)
        t = target - t_com
        q = query - q_com
        r = kabsch(target=t, query=q, weights=w, concentric=True)
        q = np.dot(q, r)
        _diff = q - t
        _current = rmsd(_diff, weights=w)
        if np.abs(_current - _rmsd) < _TINY:
            return _rmsd, r, t_com, q_com
        _rmsd = _current
        w = np.exp(-np.sum(_diff**2, axis=1) / max(_rmsd, 2.0)).tolist()
    else:
        raise Exception(f"Dynamic Kabsch did not converge in {GAUSS_MAX_ITER} steps.")


def smart_flatten(l: List[str]) -> List[int]:
    """
    Function which expands and flattens a list of integers.
    m-n -> m, m+1, ..., n
    """
    fl = []
    for i in l:
        if "-" in i:
            j = i.split("-")
            if len(j) != 2:
                raise Exception(f"Invalid range syntax: {l}")

            # temporary fix, nie bedzie wiecej bledow z tego powodu,
            # ale nie do konca to dziala tak jak powinno
            beg = int(re.sub(r"[^0-9]", "", j[0]))
            end = int(re.sub(r"[^0-9]", "", j[1]))

            if beg > end:
                raise Exception(
                    f"The left index({beg}) is greater than the right({end})"
                )
            for k in range(beg, end + 1):
                fl.append(k)
        else:
            fl.append(int(i))
    return fl


def check_peptide_sequence(sequence: str) -> bool:
    """
    Checks the peptide sequence for non-standard AAs.

    Arguments:
        sequence: the peptide sequence.

    Returns:
        True if the sequence does not contain non-standard AAs.

    Raises:
        Exception: if the sequence contains non-standard AAs.
    """
    standard_one_letter_residues = AA_NAMES.keys()
    wrong_residues = set(sequence) - set(standard_one_letter_residues)
    if wrong_residues:
        raise Exception(
            f"The input peptide sequence contains a non-standard residue(s): {' '.join(wrong_residues)}. "
            f"Only the 20 standard amino acid symbols are allowed."
        )
    return True


def fix_residue(residue: str) -> str:
    """
    Fixes non-standard AA residues in the receptor.

    Arguments:
        residue: three-letter residue code

    Returns:
        If the residue is non-standard the method returns three-letter code
        of the appropriate substitution.

    Raises:
        Exception: if the residue is non-standard and there is no substitution available.
    """
    standard_three_letter_residues = AA_NAMES.values()
    known_non_standard_three_letter_residues = AA_SUB_NAMES_EXTENDED.keys()
    if residue in standard_three_letter_residues:
        return residue
    elif residue in known_non_standard_three_letter_residues:
        modified_single = AA_SUB_NAMES_EXTENDED[residue]
        modified = AA_NAMES[modified_single]
        warnings.warn(
            f'In the current version residue "{residue}" is not supported.'
            f'"{residue}" was replaced with "{modified}" to perform the simulation.',
            UserWarning,
        )
        return modified
    else:
        raise Exception(f'The PDB file contains unknown residue "{residue}"')


def _chunk_lst(
    lst: List[Any], sl_len: int, extend_last: Optional[Any] = None
) -> List[List[Any]]:
    """Slices given list for slices of given len.

    Arguments:
        lst: list to be sliced.
        sl_len: len of one slice.
        extend_last: value to be put in last slice in order to extend it to proper length.
    """
    slists = []
    while lst:
        slists.append(lst[:sl_len])
        lst = lst[sl_len:]
    if extend_last is not None:
        _extend_last(slists, sl_len, extend_last)
    return slists


def _extend_last(sseries: List[List[Any]], slen: int, token: Any) -> None:
    """Extend the last sublist to the specified length."""
    try:
        sseries[-1].extend([token] * (slen - len(sseries[-1])))
    except IndexError:
        sseries.append([token] * slen)


def _fmt_res_name(atom: Any) -> str:
    """Format residue name from atom."""
    return (atom.chid + str(atom.resnum) + atom.icode).strip()


def pep2pep1(_id: str) -> str:
    """Convert peptide ID."""
    if re.search("PEP$", _id):
        return _id + "1"
    else:
        return _id


class EprGenerator:
    """Class that generates EPR coefficients for peptide residues based on external output (e.g. SASA)."""

    def __init__(self, file: Optional[str] = None) -> None:
        super().__init__()
        self.external_file = file
        self.variant_dictionary = {
            "easy": self.easy,
            "normalized": self.normalized,
            "multiplied": self.multiplied,
        }
        self.residues: List[str] = []
        self.external_weights: List[float] = []
        self.external_parameters: List[float] = []
        self.outfile: str = ""

    def load_file(self, file: str = "") -> None:
        """Load external parameters from file."""
        self.residues = []
        self.external_weights = []
        self.external_parameters = []
        with open(file) as file_handle:
            for line in file_handle:
                residue, external_weight, external_parameter = line.split()
                self.residues.append(residue)
                self.external_weights.append(float(external_weight))
                self.external_parameters.append(float(external_parameter))

    def generate_epr(
        self,
        variant: Literal["easy", "normalized", "multiplied"] = "easy",
        outfile: str = "outfile.txt",
    ) -> None:
        """Generate EPR coefficients using specified variant."""
        self.outfile = outfile
        self.variant_dictionary[variant]()

    def save_weights_file(
        self, outfile: str = "", new_weights: Optional[List[float]] = None
    ) -> None:
        """Save weights to file."""
        if new_weights is None:
            new_weights = []
        with open(self.outfile, "w") as outfile_handle:
            for residue, new_weight, external_parameter in zip(
                self.residues, new_weights, self.external_parameters
            ):
                outfile_handle.write(f"{residue}\t{new_weight}\t{external_parameter}\n")

    def easy(self) -> None:
        """Easy variant - use external weights as is."""
        new_weights = self.external_weights
        self.save_weights_file(outfile=self.outfile, new_weights=new_weights)

    def normalized(self) -> None:
        """Normalized variant - normalize weights."""
        new_weights = (
            len(self.external_weights)
            * np.array(self.external_weights)
            / sum(self.external_weights)
        )
        self.save_weights_file(outfile=self.outfile, new_weights=new_weights.tolist())

    def multiplied(self) -> None:
        """Multiplied variant - double the weights."""
        new_weights = 2 * np.array(self.external_weights)
        self.save_weights_file(outfile=self.outfile, new_weights=new_weights.tolist())


def random_rotation_matrix() -> RotationMatrix:
    """Generate a random rotation matrix."""
    cos = np.cos
    sin = np.sin

    x = np.random.uniform(0, np.pi * 2)
    y = np.random.uniform(0, np.pi * 2)
    z = np.random.uniform(0, np.pi * 2)

    rx = np.array([[1, 0, 0], [0, cos(x), -sin(x)], [0, sin(x), cos(x)]])
    ry = np.array([[cos(y), 0, sin(y)], [0, 1, 0], [-sin(y), 0, cos(y)]])
    rz = np.array([[cos(z), -sin(z), 0], [sin(z), cos(z), 0], [0, 0, 1]])

    return (rx.dot(ry)).dot(rz)


def convert_cg_to_all(
    filename: Union[str, TextIO],
    work_dir: str = ".",
    iter: int = 0,
    reference_pdb: Optional[str] = None,
    renumber_flag: bool = False,
    env_prefix: Optional[str] = None
) -> str:
    """
    Convert coarse-grained model to all-atom
    """
    with NamedTemporaryFile(
        prefix=".", suffix=".pdb", dir=work_dir, mode="w", delete=False
    ) as tmp_file:
        pdb = tmp_file.name

        atoms = []
        pattern = re.compile("ATOM.{9}CA .([A-Z]{3}) ([A-Z ])(.{5}).{27}(.{12}).*")
        with closing(filename) as f:
            for line in f:
                if line.startswith("ENDMDL"):
                    break
                else:
                    match = re.match(pattern, line)
                    if match:
                        atoms.append(match.groups())
                        tmp_file.write(line)

    if renumber_flag:
        sync_residues(input_pdb_path=Path(reference_pdb), output_pdb_path=Path(pdb))

    output_dir = Path(work_dir) / "output_pdbs"
    input_pdb = Path(pdb)
    fout = f"model_{iter}.pdb"
    # Modify the subprocess call to use the specified environment's python executable
    if env_prefix:
        # Use the specific executable from the isolated environment's 'bin' directory
        executable_path = os.path.join(env_prefix, "bin", "convert_cg2all")
    else:
        # Fallback to the default path if no specific env is passed
        executable_path = "convert_cg2all"

    command_parts = [
        executable_path,
        "-p", str(input_pdb),
        "-o", str(output_dir / fout),
        "--device", "cpu"
    ]

    try:
        result = subprocess.run(
            command_parts,
            shell=False,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.critical(
            module_name="CG2ALL",
            msg=f"CG2ALL failed with exit code: {e.returncode} and error: {e.stderr}",
        )
        raise Exception("CG2ALL failed to convert CG model to all-atom model")
    except Exception as e:
        logger.warning(module_name="CG2ALL", msg=f"CG2ALL failed with error: {e}")
        raise Exception("CG2ALL failed to convert CG model to all-atom model")
    finally:
        os.remove(pdb)


def sync_residues(input_pdb_path: Path, output_pdb_path: Path) -> str:
    """Synchronize residue numbering between input and output PDB files."""
    if not HAS_PANDAS:
        raise ImportError(
            "pandas and biopandas are required for sync_residues function"
        )

    input_pdb = PandasPdb().read_pdb(str(input_pdb_path))
    output_pdb = PandasPdb().read_pdb(str(output_pdb_path))
    input_atom_df = input_pdb.df["ATOM"]
    hetatm_df = input_pdb.df["HETATM"]
    hetatm_df["residue_name"] = hetatm_df["residue_name"].map(AA_SUB_NAMES)
    input_atom_df = pd.concat([input_atom_df, hetatm_df]).dropna(
        subset=["residue_name"]
    )
    output_atom_df = output_pdb.df["ATOM"]
    output_atom_df.loc[output_atom_df["atom_name"] == "CA", "residue_number"] = (
        input_atom_df.loc[input_atom_df["atom_name"] == "CA", "residue_number"].values
    )
    output_pdb.df["ATOM"] = output_atom_df
    output_pdb.to_pdb(str(output_pdb_path))
    return "Residues synchronized"
