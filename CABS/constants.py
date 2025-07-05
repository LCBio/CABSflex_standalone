"""
Constants and enums for the CABS package.

This module contains constant data including amino acid mappings,
secondary structure codes, side chain coordinates, and configuration templates.
"""

from enum import Enum
from typing import Dict, Final, List, Literal, Tuple, Union

import numpy as np
import numpy.typing as npt

try:
    from importlib.resources import as_file, files
except ImportError:
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=DeprecationWarning, module="pkg_resources"
        )
        from pkg_resources import resource_filename

# Import configuration loader
from CABS.config_loader import (
    get_aa_names,
    get_aa_sub_names,
    get_aa_sub_names_extended,
    get_allowed_aa_methods,
    get_bfac_output_options,
    get_cabs_files,
    get_cabs_lattice_defaults,
    get_cabs_ss,
    get_cabs_ss_reverse,
    get_calculation_constants,
    get_config_header,
    get_contact_map_constants,
    get_csv_output_options,
    get_default_colors,
    get_default_filenames,
    get_default_peptide_ss,
    get_default_values,
    get_dssp_ss_mapping,
    get_error_messages,
    get_file_extensions,
    get_math_constants,
    get_model_constants,
    get_nsp3_constants,
    get_output_directories,
    get_pdb_output_options,
    get_peptide_replacements,
    get_protein_category_modes,
    get_protein_restraints_modes,
    get_sc_modeling_thresholds,
    get_sidecnt,
    get_string_patterns,
    get_system_constants,
    get_unleashed_aliases,
    get_valid_letters,
)

# Type aliases for better type hints
AminoAcidCode = Literal[
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]
SecondaryStructureCode = Literal["C", "H", "T", "E", "c", "h", "t", "e"]
ColorHex = str  # Hex color string like '#ffffff'


class SecondaryStructure(Enum):
    """Secondary structure types with CABS encoding."""

    COIL = 1
    HELIX = 2
    TURN = 3
    STRAND = 4


class AminoAcid(Enum):
    """Standard amino acids with their properties."""

    ALA = ("A", "Alanine", "Ala")
    CYS = ("C", "Cysteine", "Cys")
    ASP = ("D", "Aspartic acid", "Asp")
    GLU = ("E", "Glutamic acid", "Glu")
    PHE = ("F", "Phenylalanine", "Phe")
    GLY = ("G", "Glycine", "Gly")
    HIS = ("H", "Histidine", "His")
    ILE = ("I", "Isoleucine", "Ile")
    LYS = ("K", "Lysine", "Lys")
    LEU = ("L", "Leucine", "Leu")
    MET = ("M", "Methionine", "Met")
    ASN = ("N", "Asparagine", "Asn")
    PRO = ("P", "Proline", "Pro")
    GLN = ("Q", "Glutamine", "Gln")
    ARG = ("R", "Arginine", "Arg")
    SER = ("S", "Serine", "Ser")
    THR = ("T", "Threonine", "Thr")
    VAL = ("V", "Valine", "Val")
    TRP = ("W", "Tryptophan", "Trp")
    TYR = ("Y", "Tyrosine", "Tyr")

    def __init__(self, single: str, full_name: str, three_letter: str) -> None:
        self.single = single
        self.full_name = full_name
        self.three_letter = three_letter


# Dictionary for conversion of secondary structure from DSSP to CABS
CABS_SS: Final[Dict[SecondaryStructureCode, int]] = get_cabs_ss()

CABS_SS_REVERSE: Final[Dict[int, str]] = {
    int(k): v for k, v in get_cabs_ss_reverse().items()
}

# Side chain relative coordinates
SIDECNT: Final[Dict[str, Tuple[float, ...]]] = {
    k: tuple(v) for k, v in get_sidecnt().items()
}

# Create amino acid lookup dictionaries
AA_NAMES: Final[Dict[AminoAcidCode, str]] = get_aa_names()
AA_SUB_NAMES: Final[Dict[str, AminoAcidCode]] = get_aa_sub_names()


# Load random ligand library
def _load_random_ligand_library() -> npt.NDArray[np.float64]:
    """Load the random ligand library from data file."""
    try:
        # Try modern importlib.resources first
        try:
            with as_file(files("CABS") / "data" / "data2.dat") as data_file:
                return np.reshape(np.fromfile(str(data_file), sep=" "), (1000, 50, 3))
        except (ImportError, AttributeError):
            # Fallback to pkg_resources
            data_file = resource_filename("CABS", "data/data2.dat")
            return np.reshape(np.fromfile(data_file, sep=" "), (1000, 50, 3))
    except Exception:
        # Return zeros if data file cannot be loaded
        return np.zeros((1000, 50, 3))


RANDOM_LIGAND_LIBRARY: Final[npt.NDArray[np.float64]] = _load_random_ligand_library()

# Extended amino acid substitution dictionary (non-standard amino acids)
AA_SUB_NAMES_EXTENDED: Final[Dict[str, AminoAcidCode]] = get_aa_sub_names_extended()

# File extensions and formats
_file_ext_config = get_file_extensions()
PDB_EXTENSIONS: Final[List[str]] = _file_ext_config["pdb"]
IMAGE_FORMATS: Final[List[str]] = _file_ext_config["image_formats"]
DEFAULT_IMAGE_FORMAT: Final[str] = _file_ext_config["default_image_format"]

# CABS-specific constants
CABS_LATTICE_DEFAULTS: Final[Dict[str, Union[float, Tuple[float, float]]]] = (
    get_cabs_lattice_defaults()
)

# Side chain modeling constants
SC_MODELING_THRESHOLDS: Final[Dict[str, float]] = get_sc_modeling_thresholds()

# PDB output options
PDB_OUTPUT_OPTIONS: Final[Dict[str, str]] = get_pdb_output_options()

# Beta factor output options
BFAC_OUTPUT_OPTIONS: Final[Dict[str, str]] = get_bfac_output_options()

# CSV output options
CSV_OUTPUT_OPTIONS: Final[Dict[str, str]] = get_csv_output_options()

# Valid letters for different output types
_valid_letters = get_valid_letters()
VALID_PDB_OUTPUT_LETTERS: Final[str] = _valid_letters["pdb_output"]
VALID_BFAC_OUTPUT_LETTERS: Final[str] = _valid_letters["bfac_output"]
VALID_CSV_OUTPUT_LETTERS: Final[str] = _valid_letters["csv_output"]

# Protein restraints modes
PROTEIN_RESTRAINTS_MODES: Final[List[str]] = get_protein_restraints_modes()
PROTEIN_CATEGORY_MODES: Final[List[str]] = get_protein_category_modes()

# Special protein restraints aliases
UNLEASHED_ALIASES: Final[List[str]] = get_unleashed_aliases()

# Amino acid reconstruction methods
ALLOWED_AA_METHODS: Final[List[str]] = get_allowed_aa_methods()

# CABS files that are generated during simulation
CABS_FILES: Final[List[str]] = get_cabs_files()

# DSSP secondary structure mapping
DSSP_SS_MAPPING: Final[Dict[str, str]] = get_dssp_ss_mapping()

# Default colors for contact maps and plots
DEFAULT_COLORS: Final[List[ColorHex]] = get_default_colors()

# Contact map histogram constants
CONTACT_MAP_CONSTANTS: Final[Dict[str, Union[int, float]]] = get_contact_map_constants()

# Peptide replacement patterns
PEPTIDE_REPLACEMENTS: Final[Dict[str, str]] = get_peptide_replacements()

# File paths and directories
OUTPUT_DIRECTORIES: Final[Dict[str, str]] = get_output_directories()

# Default file names
DEFAULT_FILENAMES: Final[Dict[str, str]] = get_default_filenames()

# Calculation thresholds and limits
CALCULATION_CONSTANTS: Final[Dict[str, Union[int, float]]] = get_calculation_constants()

# Legacy constants for backward compatibility
_calc_constants = CALCULATION_CONSTANTS
_LARGE: Final[float] = _calc_constants["large_value"]
_TINY: Final[float] = _calc_constants["tiny_value"]
GAUSS_MAX_ITER: Final[int] = _calc_constants["gauss_max_iter"]

# Mathematical constants
MATH_CONSTANTS: Final[Dict[str, Union[int, float]]] = get_math_constants()

# String patterns and templates
STRING_PATTERNS: Final[Dict[str, str]] = get_string_patterns()

# System and platform constants
SYSTEM_CONSTANTS: Final[Dict[str, Union[str, List[str]]]] = get_system_constants()

# Model file extensions and paths
MODEL_CONSTANTS: Final[Dict[str, Union[str, List[str]]]] = get_model_constants()

# Configuration header template
CONFIG_HEADER: Final[str] = get_config_header()

# Error messages and warnings
ERROR_MESSAGES: Final[Dict[str, str]] = get_error_messages()

# Default values for various parameters
DEFAULT_VALUES: Final[Dict[str, Union[int, float, str, bool]]] = get_default_values()

# NetSurfP-3.0 related constants
NSP3_CONSTANTS: Final[Dict[str, str]] = get_nsp3_constants()

# Default secondary structure for unknown peptides
DEFAULT_PEPTIDE_SS: Final[str] = get_default_peptide_ss()
