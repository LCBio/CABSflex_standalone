"""
Constants and enums for the CABS package.

This module contains all constant data extracted from various CABS modules,
including amino acid mappings, secondary structure codes, side chain coordinates,
and configuration templates.
"""

from enum import Enum
from typing import Dict, List, Tuple, Union, Final
from typing_extensions import Literal
import numpy as np
import numpy.typing as npt

try:
    from importlib.resources import files, as_file
except ImportError:
    # Fallback for Python < 3.9
    from pkg_resources import resource_filename

# Type aliases for better type hints
AminoAcidCode = Literal['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
SecondaryStructureCode = Literal['C', 'H', 'T', 'E', 'c', 'h', 't', 'e']
ColorHex = str  # Hex color string like '#ffffff'

class SecondaryStructure(Enum):
    """Secondary structure types with CABS encoding."""
    COIL = 1
    HELIX = 2
    TURN = 3
    STRAND = 4


class AminoAcid(Enum):
    """Standard amino acids with their properties."""
    ALA = ('A', 'Alanine', 'Ala')
    CYS = ('C', 'Cysteine', 'Cys')
    ASP = ('D', 'Aspartic acid', 'Asp')
    GLU = ('E', 'Glutamic acid', 'Glu')
    PHE = ('F', 'Phenylalanine', 'Phe')
    GLY = ('G', 'Glycine', 'Gly')
    HIS = ('H', 'Histidine', 'His')
    ILE = ('I', 'Isoleucine', 'Ile')
    LYS = ('K', 'Lysine', 'Lys')
    LEU = ('L', 'Leucine', 'Leu')
    MET = ('M', 'Methionine', 'Met')
    ASN = ('N', 'Asparagine', 'Asn')
    PRO = ('P', 'Proline', 'Pro')
    GLN = ('Q', 'Glutamine', 'Gln')
    ARG = ('R', 'Arginine', 'Arg')
    SER = ('S', 'Serine', 'Ser')
    THR = ('T', 'Threonine', 'Thr')
    VAL = ('V', 'Valine', 'Val')
    TRP = ('W', 'Tryptophan', 'Trp')
    TYR = ('Y', 'Tyrosine', 'Tyr')

    def __init__(self, single: str, full_name: str, three_letter: str) -> None:
        self.single = single
        self.full_name = full_name
        self.three_letter = three_letter


# Dictionary for conversion of secondary structure from DSSP to CABS
CABS_SS: Final[Dict[SecondaryStructureCode, int]] = {
    'C': SecondaryStructure.COIL.value,
    'H': SecondaryStructure.HELIX.value,
    'T': SecondaryStructure.TURN.value,
    'E': SecondaryStructure.STRAND.value,
    'c': SecondaryStructure.COIL.value,
    'h': SecondaryStructure.HELIX.value,
    't': SecondaryStructure.TURN.value,
    'e': SecondaryStructure.STRAND.value
}

CABS_SS_REVERSE: Final[Dict[int, str]] = {
    SecondaryStructure.COIL.value: 'C',
    SecondaryStructure.HELIX.value: 'H',
    SecondaryStructure.TURN.value: 'T',
    SecondaryStructure.STRAND.value: 'E'
}

# Side chain relative coordinates
SIDECNT: Final[Dict[str, Tuple[float, ...]]] = {
    'ALA': (-0.464, -0.464, 1.073, 0.0, 0.0, 0.0),
    'ARG': (2.103, -0.479, 3.266, 0.931, -0.479, 4.281),
    'ASN': (-0.402, -1.237, 2.111, 0.132, -0.863, 2.328),
    'ASP': (-0.391, -1.358, 1.927, 0.151, -0.914, 1.927),
    'CYS': (-0.464, -0.866, 1.545, 0.0, 0.0, 0.0),
    'GLN': (0.481, -0.766, 2.795, 1.015, -0.392, 3.810),
    'GLU': (0.470, -0.887, 2.611, 1.004, -0.443, 3.625),
    'GLY': (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    'HIS': (0.394, -0.688, 2.598, 0.928, -0.314, 3.613),
    'ILE': (-0.175, -0.175, 2.120, 0.289, 0.289, 2.120),
    'LEU': (0.708, -0.175, 2.120, 1.197, 0.289, 2.120),
    'LYS': (1.620, -0.105, 3.450, 2.154, 0.269, 4.465),
    'MET': (1.101, -0.392, 2.314, 1.635, 0.082, 3.329),
    'PHE': (0.394, -0.688, 2.598, 0.928, -0.314, 3.613),
    'PRO': (-0.175, -0.175, 1.547, 0.289, 0.289, 1.547),
    'SER': (-0.464, -0.866, 1.545, 0.0, 0.0, 0.0),
    'THR': (-0.175, -0.577, 1.895, 0.289, -0.103, 1.895),
    'TRP': (0.394, -0.688, 2.598, 0.928, -0.314, 3.613),
    'TYR': (0.308, -1.387, 3.492, -0.618, -0.799, 3.634),
    'VAL': (-0.175, -0.175, 1.547, 0.289, 0.289, 1.547)
}

# Create amino acid lookup dictionaries
AA_NAMES: Final[Dict[AminoAcidCode, str]] = {aa.single: aa.three_letter for aa in AminoAcid}
AA_SUB_NAMES: Final[Dict[str, AminoAcidCode]] = {aa.three_letter: aa.single for aa in AminoAcid}

# Load random ligand library
def _load_random_ligand_library() -> npt.NDArray[np.float64]:
    """Load the random ligand library from data file."""
    try:
        # Try modern importlib.resources first
        try:
            with as_file(files('CABS') / 'data' / 'data2.dat') as data_file:
                return np.reshape(np.fromfile(str(data_file), sep=' '), (1000, 50, 3))
        except (ImportError, AttributeError):
            # Fallback to pkg_resources
            data_file = resource_filename('CABS', 'data/data2.dat')
            return np.reshape(np.fromfile(data_file, sep=' '), (1000, 50, 3))
    except Exception:
        # Return zeros if data file cannot be loaded
        return np.zeros((1000, 50, 3))

RANDOM_LIGAND_LIBRARY: Final[npt.NDArray[np.float64]] = _load_random_ligand_library()

# Extended amino acid substitution dictionary (non-standard amino acids)
AA_SUB_NAMES_EXTENDED: Final[Dict[str, AminoAcidCode]] = {
    **AA_SUB_NAMES,
    '0CS': 'A',  # 0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
    '1AB': 'P',  # 1AB PRO  1,4-DIDEOXY-1,4-IMINO-D-ARABINITOL
    '1LU': 'L',  # 1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
    '1PA': 'F',  # 1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
    '1TQ': 'W',  # 1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
    '1TY': 'Y',  # 1TY TYR
    '23F': 'F',  # 23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
    '23S': 'W',  # 23S TRP  MODIFIED TRYPTOPHAN
    '2BU': 'A',  # 2BU ADE
    '2ML': 'L',  # 2ML LEU  2-METHYLLEUCINE
    '2MR': 'R',  # 2MR ARG  N3, N4-DIMETHYLARGININE
    '2MT': 'P',  # 2MT PRO
    '2OP': 'A',  # 2OP (2S  2-HYDROXYPROPANAL
    '2TY': 'Y',  # 2TY TYR
    '32S': 'W',  # 32S TRP  MODIFIED TRYPTOPHAN
    '32T': 'W',  # 32T TRP  MODIFIED TRYPTOPHAN
    '3AH': 'H',  # 3AH HIS
    '3MD': 'D',  # 3MD ASP  2S,3S-3-METHYLASPARTIC ACID
    '3TY': 'Y',  # 3TY TYR  MODIFIED TYROSINE
    '4DP': 'W',  # 4DP TRP
    '4F3': 'A',  # 4F3 ALA  CYCLIZED
    '4FB': 'P',  # 4FB PRO  (4S)-4-FLUORO-L-PROLINE
    '4FW': 'W',  # 4FW TRP  4-FLUOROTRYPTOPHANE
    '4HT': 'W',  # 4HT TRP  4-HYDROXYTRYPTOPHAN
    '4IN': 'W',  # 4IN TRP  4-AMINO-L-TRYPTOPHAN
    '4PH': 'F',  # 4PH PHE  4-METHYL-L-PHENYLALANINE
    '5CS': 'C',  # 5CS CYS
    '6CL': 'K',  # 6CL LYS  6-CARBOXYLYSINE
    '6CW': 'W',  # 6CW TRP  6-CHLORO-L-TRYPTOPHAN
    'A0A': 'D',  # A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
    'AA4': 'A',  # AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
    'AAR': 'R',  # AAR ARG  ARGININEAMIDE
    'AB7': 'E',  # AB7 GLU  ALPHA-AMINOBUTYRIC ACID
    'ABA': 'A',  # ABA ALA  ALPHA-AMINOBUTYRIC ACID
    'ACB': 'D',  # ACB ASP  3-METHYL-ASPARTIC ACID
    'ACL': 'R',  # ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
    'ACY': 'G',  # ACY GLY  POST-TRANSLATIONAL MODIFICATION
    'AEI': 'T',  # AEI THR  ACYLATED THR
    'AFA': 'N',  # AFA ASN  N-[7-METHYL-OCT-2,4-DIENOYL]ASPARAGINE
    'AGM': 'R',  # AGM ARG  4-METHYL-ARGININE
    'AGT': 'C',  # AGT CYS  AGMATINE-CYSTEINE ADDUCT
    'AHB': 'N',  # AHB ASN  BETA-HYDROXYASPARAGINE
    'AHO': 'A',  # AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
    'AHP': 'A',  # AHP ALA  2-AMINO-HEPTANOIC ACID
    'AIB': 'A',  # AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
    'AKL': 'D',  # AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
    'ALC': 'A',  # ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
    'ALG': 'R',  # ALG ARG  GUANIDINOBUTYRYL GROUP
    'ALM': 'A',  # ALM ALA  1-METHYL-ALANINAL
    'ALN': 'A',  # ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
    'ALO': 'T',  # ALO THR  ALLO-THREONINE
    'ALS': 'A',  # ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
    'ALT': 'A',  # ALT ALA  THIOALANINE
    'ALY': 'K',  # ALY LYS  N(6)-ACETYLLYSINE
    'AME': 'M',  # AME MET  ACETYLATED METHIONINE
    'AP7': 'A',  # AP7 ADE
    'APH': 'A',  # APH ALA  P-AMIDINOPHENYL-3-ALANINE
    'API': 'K',  # API LYS  2,6-DIAMINOPIMELIC ACID
    'APK': 'K',  # APK LYS
    'AR2': 'R',  # AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
    'AR4': 'E',  # AR4 GLU
    'ARM': 'R',  # ARM ARG  DEOXY-METHYL-ARGININE
    'ARO': 'R',  # ARO ARG  C-GAMMA-HYDROXY ARGININE
    'ASA': 'D',  # ASA ASP  ASPARTIC ALDEHYDE
    'ASB': 'D',  # ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
    'ASI': 'D',  # ASI ASP  L-ISO-ASPARTATE
    'ASK': 'D',  # ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
    'ASL': 'D',  # ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
    'AYA': 'A',  # AYA ALA  N-ACETYLALANINE
    'AYG': 'A',  # AYG ALA
    'AZK': 'K',  # AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
    'B2A': 'A',  # B2A ALA  ALANINE BORONIC ACID
    'B2F': 'F',  # B2F PHE  PHENYLALANINE BORONIC ACID
    'B2I': 'I',  # B2I ILE  ISOLEUCINE BORONIC ACID
    'B2V': 'V',  # B2V VAL  VALINE BORONIC ACID
    'B3A': 'A',  # B3A ALA  (3S)-3-AMINOBUTANOIC ACID
    'B3D': 'D',  # B3D ASP  3-AMINOPENTANEDIOIC ACID
    'B3E': 'E',  # B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
    'B3K': 'K',  # B3K LYS  (3S)-3,7-DIAMINOHEPTANOIC ACID
    'B3S': 'S',  # B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
    'B3X': 'N',  # B3X ASN  (3S)-3,5-DIAMINO-5-OXOPENTANOIC ACID
    'B3Y': 'Y',  # B3Y TYR
    'BAL': 'A',  # BAL ALA  BETA-ALANINE
    'BBC': 'C',  # BBC CYS
    'BCS': 'C',  # BCS CYS  BENZYLCYSTEINE
    'BCX': 'C',  # BCX CYS  BETA-3-CYSTEINE
    'BFD': 'D',  # BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
    'BG1': 'S',  # BG1 SER
    'BHD': 'D',  # BHD ASP  BETA-HYDROXYASPARTIC ACID
    'BIF': 'F',  # BIF PHE
    'BLE': 'L',  # BLE LEU  LEUCINE BORONIC ACID
    'BLY': 'K',  # BLY LYS  LYSINE BORONIC ACID
    'BMT': 'T',  # BMT THR
    'BNN': 'A',  # BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
    'BOR': 'R',  # BOR ARG
    'BPE': 'C',  # BPE CYS
    'BTR': 'W',  # BTR TRP  6-BROMO-TRYPTOPHAN
    'BUC': 'C',  # BUC CYS  S,S-BUTYLTHIOCYSTEINE
    'BUG': 'L',  # BUG LEU  TERT-LEUCYL AMINE
    'C12': 'A',  # C12 ALA
    'C1X': 'K',  # C1X LYS  MODIFIED LYSINE
    'C3Y': 'C',  # C3Y CYS  MODIFIED CYSTEINE
    'C5C': 'C',  # C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
    'C6C': 'C',  # C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
    'C99': 'A',  # C99 ALA
    'CAB': 'A',  # CAB ALA  4-CARBOXY-4-AMINOBUTANAL
    'CAF': 'C',  # CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
    'CAS': 'C',  # CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
    'CCS': 'C',  # CCS CYS  CARBOXYMETHYLATED CYSTEINE
    'CGU': 'E',  # CGU GLU  CARBOXYLATION OF THE CG ATOM
    'CH6': 'A',  # CH6 ALA
    'CH7': 'A',  # CH7 ALA
    'CHG': 'G',  # CHG GLY  CYCLOHEXYL GLYCINE
    'CHP': 'G',  # CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
    'CHS': 'F',  # CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
    'CIR': 'R',  # CIR ARG  CITRULLINE
    'CLB': 'A',  # CLB ALA
    'CLD': 'A',  # CLD ALA
    'CLE': 'L',  # CLE LEU  LEUCINE AMIDE
    'CLG': 'K',  # CLG LYS
    'CLH': 'K',  # CLH LYS
    'CLV': 'A',  # CLV ALA
    'CME': 'C',  # CME CYS  MODIFIED CYSTEINE
    'CML': 'C',  # CML CYS
    'CMT': 'C',  # CMT CYS  O-METHYLCYSTEINE
    'CQR': 'A',  # CQR ALA
    'CR2': 'A',  # CR2 ALA  POST-TRANSLATIONAL MODIFICATION
    'CR5': 'A',  # CR5 ALA
    'CR7': 'A',  # CR7 ALA
    'CR8': 'A',  # CR8 ALA
    'CRK': 'A',  # CRK ALA
    'CRO': 'T',  # CRO THR  CYCLIZED
    'CRQ': 'Y',  # CRQ TYR
    'CRW': 'A',  # CRW ALA
    'CRX': 'A',  # CRX ALA
    'CS1': 'C',  # CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
    'CS3': 'C',  # CS3 CYS
    'CS4': 'C',  # CS4 CYS
    'CSA': 'C',  # CSA CYS  S-ACETONYLCYSTEIN
    'CSB': 'C',  # CSB CYS  CYS BOUND TO LEAD ION
    'CSD': 'C',  # CSD CYS  3-SULFINOALANINE
    'CSE': 'C',  # CSE CYS  SELENOCYSTEINE
    'CSI': 'A',  # CSI ALA
    'CSO': 'C',  # CSO CYS  INE S-HYDROXYCYSTEINE
    'CSR': 'C',  # CSR CYS  S-ARSONOCYSTEINE
    'CSS': 'C',  # CSS CYS  1,3-THIAZOLE-4-CARBOXYLIC ACID
    'CSU': 'C',  # CSU CYS  CYSTEINE-S-SULFONIC ACID
    'CSW': 'C',  # CSW CYS  CYSTEINE-S-DIOXIDE
    'CSX': 'C',  # CSX CYS  OXOCYSTEINE
    'CSY': 'A',  # CSY ALA  MODIFIED TYROSINE COMPLEX
    'CSZ': 'C',  # CSZ CYS  S-SELANYL CYSTEINE
    'CTH': 'T',  # CTH THR  4-CHLOROTHREONINE
    'CWR': 'A',  # CWR ALA
    'CXM': 'M',  # CXM MET  N-CARBOXYMETHIONINE
    'CY0': 'C',  # CY0 CYS  MODIFIED CYSTEINE
    'CY1': 'C',  # CY1 CYS  ACETAMIDOMETHYLCYSTEINE
    'CY3': 'C',  # CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
    'CY4': 'C',  # CY4 CYS  S-BUTYRYL-CYSTEIN
    'CY7': 'C',  # CY7 CYS  MODIFIED CYSTEINE
    'CYD': 'C',  # CYD CYS
    'CYF': 'C',  # CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
    'CYG': 'C',  # CYG CYS
    'CYJ': 'K',  # CYJ LYS  MODIFIED LYSINE
    'CYQ': 'C',  # CYQ CYS
    'CYR': 'C',  # CYR CYS
    'CZ2': 'C',  # CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
    'CZZ': 'C',  # CZZ CYS  THIARSAHYDROXY-CYSTEINE
    'DA2': 'R',  # DA2 ARG  MODIFIED ARGININE
    'DAB': 'A',  # DAB ALA  2,4-DIAMINOBUTYRIC ACID
    'DAH': 'F',  # DAH PHE  3,4-DIHYDROXYDAHNYLALANINE
    'DAL': 'A',  # DAL ALA  D-ALANINE
    'DAM': 'A',  # DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
    'DAR': 'R',  # DAR ARG  D-ARGININE
    'DAS': 'D',  # DAS ASP  D-ASPARTIC ACID
    'DBU': 'A',  # DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
    'DBY': 'Y',  # DBY TYR  3,5 DIBROMOTYROSINE
    'DBZ': 'A',  # DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
    'DCL': 'L',  # DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
    'DCY': 'C',  # DCY CYS  D-CYSTEINE
    'DDE': 'H',  # DDE HIS
    'DGL': 'E',  # DGL GLU  D-GLU
    'DGN': 'Q',  # DGN GLN  D-GLUTAMINE
    'DHA': 'A',  # DHA ALA  2-AMINO-ACRYLIC ACID
    'DHI': 'H',  # DHI HIS  D-HISTIDINE
    'DHL': 'S',  # DHL SER  POST-TRANSLATIONAL MODIFICATION
    'DIL': 'I',  # DIL ILE  D-ISOLEUCINE
    'DIV': 'V',  # DIV VAL  D-ISOVALINE
    'DLE': 'L',  # DLE LEU  D-LEUCINE
    'DLS': 'K',  # DLS LYS  DI-ACETYL-LYSINE
    'DLY': 'K',  # DLY LYS  D-LYSINE
    'DMH': 'N',  # DMH ASN  N4,N4-DIMETHYL-ASPARAGINE
    'DMK': 'D',  # DMK ASP  DIMETHYL ASPARTIC ACID
    'DNE': 'L',  # DNE LEU  D-NORLEUCINE
    'DNG': 'L',  # DNG LEU  N-FORMYL-D-NORLEUCINE
    'DNL': 'K',  # DNL LYS  6-AMINO-HEXANAL
    'DNM': 'L',  # DNM LEU  D-N-METHYL NORLEUCINE
    'DPH': 'F',  # DPH PHE  DEAMINO-METHYL-PHENYLALANINE
    'DPL': 'P',  # DPL PRO  4-OXOPROLINE
    'DPN': 'F',  # DPN PHE  D-CONFIGURATION
    'DPP': 'A',  # DPP ALA  DIAMMINOPROPANOIC ACID
    'DPQ': 'Y',  # DPQ TYR  TYROSINE DERIVATIVE
    'DPR': 'P',  # DPR PRO  D-PROLINE
    'DSE': 'S',  # DSE SER  D-SERINE N-METHYLATED
    'DSG': 'N',  # DSG ASN  D-ASPARAGINE
    'DSN': 'S',  # DSN SER  D-SERINE
    'DTH': 'T',  # DTH THR  D-THREONINE
    'DTR': 'W',  # DTR TRP  D-TRYPTOPHAN
    'DTY': 'Y',  # DTY TYR  D-TYROSINE
    'DVA': 'V',  # DVA VAL  D-VALINE
    'DYG': 'A',  # DYG ALA
    'DYS': 'C',  # DYS CYS
    'EFC': 'C',  # EFC CYS  S,S-(2-FLUOROETHYL)THIOCYSTEINE
    'ESB': 'Y',  # ESB TYR
    'ESC': 'M',  # ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
    'FCL': 'F',  # FCL PHE  3-CHLORO-L-PHENYLALANINE
    'FGL': 'A',  # FGL ALA  2-AMINOPROPANEDIOIC ACID
    'FGP': 'S',  # FGP SER
    'FHL': 'K',  # FHL LYS  MODIFIED LYSINE
    'FLE': 'L',  # FLE LEU  FUROYL-LEUCINE
    'FLT': 'Y',  # FLT TYR  FLUOROMALONYL TYROSINE
    'FME': 'M',  # FME MET  FORMYL-METHIONINE
    'FOE': 'C',  # FOE CYS
    'FOG': 'F',  # FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
    'FOR': 'M',  # FOR MET
    'FRF': 'F',  # FRF PHE  PHE FOLLOWED BY REDUCED PHE
    'FTR': 'W',  # FTR TRP  FLUOROTRYPTOPHANE
    'FTY': 'Y',  # FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
    'GHG': 'Q',  # GHG GLN  GAMMA-HYDROXY-GLUTAMINE
    'GHP': 'G',  # GHP GLY  4-HYDROXYPHENYLGLYCINE
    'GL3': 'G',  # GL3 GLY  POST-TRANSLATIONAL MODIFICATION
    'GLH': 'Q',  # GLH GLN
    'GLZ': 'G',  # GLZ GLY  AMINO-ACETALDEHYDE
    'GMA': 'E',  # GMA GLU  1-AMIDO-GLUTAMIC ACID
    'GMU': 'A',  # GMU 5MU
    'GPL': 'K',  # GPL LYS  LYSINE GUANOSINE-5'-MONOPHOSPHATE
    'GT9': 'C',  # GT9 CYS  SG ALKYLATED
    'GVL': 'S',  # GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
    'GYC': 'C',  # GYC CYS
    'GYS': 'G',  # GYS GLY
    'H5M': 'P',  # H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
    'HHK': 'A',  # HHK ALA  (2S)-2,8-DIAMINOOCTANOIC ACID
    'HIA': 'H',  # HIA HIS  L-HISTIDINE AMIDE
    'HIC': 'H',  # HIC HIS  4-METHYL-HISTIDINE
    'HIP': 'H',  # HIP HIS  ND1-PHOSPHONOHISTIDINE
    'HIQ': 'H',  # HIQ HIS  MODIFIED HISTIDINE
    'HLU': 'L',  # HLU LEU  BETA-HYDROXYLEUCINE
    'HMF': 'A',  # HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
    'HMR': 'R',  # HMR ARG  BETA-HOMOARGININE
    'HPE': 'F',  # HPE PHE  HOMOPHENYLALANINE
    'HPH': 'F',  # HPH PHE  PHENYLALANINOL GROUP
    'HPQ': 'F',  # HPQ PHE  HOMOPHENYLALANINYLMETHANE
    'HRG': 'R',  # HRG ARG  L-HOMOARGININE
    'HSE': 'S',  # HSE SER  L-HOMOSERINE
    'HSL': 'S',  # HSL SER  HOMOSERINE LACTONE
    'HSO': 'H',  # HSO HIS  HISTIDINOL
    'HTI': 'C',  # HTI CYS
    'HTR': 'W',  # HTR TRP  BETA-HYDROXYTRYPTOPHANE
    'HY3': 'P',  # HY3 PRO  3-HYDROXYPROLINE
    'HYP': 'P',  # HYP PRO  4-HYDROXYPROLINE
    'IAM': 'A',  # IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
    'IAS': 'D',  # IAS ASP  ASPARTYL GROUP
    'IGL': 'A',  # IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
    'IIL': 'I',  # IIL ILE  ISO-ISOLEUCINE
    'ILG': 'E',  # ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
    'ILX': 'I',  # ILX ILE  4,5-DIHYDROXYISOLEUCINE
    'IML': 'I',  # IML ILE  N-METHYLATED
    'IPG': 'G',  # IPG GLY  N-ISOPROPYL GLYCINE
    'IT1': 'K',  # IT1 LYS
    'IYR': 'Y',  # IYR TYR  3-IODO-TYROSINE
    'KCX': 'K',  # KCX LYS  CARBAMOYLATED LYSINE
    'KGC': 'K',  # KGC LYS
    'KOR': 'C',  # KOR CYS  MODIFIED CYSTEINE
    'KST': 'K',  # KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
    'KYN': 'A',  # KYN ALA  KYNURENINE
    'LA2': 'K',  # LA2 LYS
    'LAL': 'A',  # LAL ALA  N,N-DIMETHYL-L-ALANINE
    'LCK': 'K',  # LCK LYS
    'LCX': 'K',  # LCX LYS  CARBAMYLATED LYSINE
    'LDH': 'K',  # LDH LYS  N~6~-ETHYL-L-LYSINE
    'LED': 'L',  # LED LEU  POST-TRANSLATIONAL MODIFICATION
    'LEF': 'L',  # LEF LEU  2-5-FLUOROLEUCINE
    'LET': 'K',  # LET LYS  ODIFIED LYSINE
    'LLP': 'K',  # LLP LYS
    'LLY': 'K',  # LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
    'LME': 'E',  # LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
    'LNT': 'L',  # LNT LEU
    'LPD': 'P',  # LPD PRO  L-PROLINAMIDE
    'LSO': 'K',  # LSO LYS  MODIFIED LYSINE
    'LYM': 'K',  # LYM LYS  DEOXY-METHYL-LYSINE
    'LYN': 'K',  # LYN LYS  2,6-DIAMINO-HEXANOIC ACID AMIDE
    'LYP': 'K',  # LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
    'LYR': 'K',  # LYR LYS  MODIFIED LYSINE
    'LYX': 'K',  # LYX LYS  N''-(2-COENZYME A)-PROPANOYL-LYSINE
    'LYZ': 'K',  # LYZ LYS  5-HYDROXYLYSINE
    'M0H': 'C',  # M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
    'M2L': 'K',  # M2L LYS
    'M3L': 'K',  # M3L LYS  N-TRIMETHYLLYSINE
    'MAA': 'A',  # MAA ALA  N-METHYLALANINE
    'MAI': 'R',  # MAI ARG  DEOXO-METHYLARGININE
    'MBQ': 'Y',  # MBQ TYR
    'MC1': 'S',  # MC1 SER  METHICILLIN ACYL-SERINE
    'MCL': 'K',  # MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
    'MCS': 'C',  # MCS CYS  MALONYLCYSTEINE
    'MDO': 'A',  # MDO ALA
    'MEA': 'F',  # MEA PHE  N-METHYLPHENYLALANINE
    'MEG': 'E',  # MEG GLU  (2S,3R)-3-METHYL-GLUTAMIC ACID
    'MEN': 'N',  # MEN ASN  GAMMA METHYL ASPARAGINE
    'MEU': 'G',  # MEU GLY  O-METHYL-GLYCINE
    'MFC': 'A',  # MFC ALA  CYCLIZED
    'MGG': 'R',  # MGG ARG  MODIFIED D-ARGININE
    'MGN': 'Q',  # MGN GLN  2-METHYL-GLUTAMINE
    'MHL': 'L',  # MHL LEU  N-METHYLATED, HYDROXY
    'MHO': 'M',  # MHO MET  POST-TRANSLATIONAL MODIFICATION
    'MHS': 'H',  # MHS HIS  1-N-METHYLHISTIDINE
    'MIS': 'S',  # MIS SER  MODIFIED SERINE
    'MLE': 'L',  # MLE LEU  N-METHYLATED
    'MLL': 'L',  # MLL LEU  METHYL L-LEUCINATE
    'MLY': 'K',  # MLY LYS  METHYLATED LYSINE
    'MLZ': 'K',  # MLZ LYS  N-METHYL-LYSINE
    'MME': 'M',  # MME MET  N-METHYL METHIONINE
    'MNL': 'L',  # MNL LEU  4,N-DIMETHYLNORLEUCINE
    'MNV': 'V',  # MNV VAL  N-METHYL-C-AMINO VALINE
    'MPQ': 'G',  # MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
    'MSA': 'G',  # MSA GLY  (2-S-METHYL) SARCOSINE
    'MSE': 'M',  # MSE MET  ELENOMETHIONINE
    'MSO': 'M',  # MSO MET  METHIONINE SULFOXIDE
    'MTY': 'F',  # MTY PHE  3-HYDROXYPHENYLALANINE
    'MVA': 'V',  # MVA VAL  N-METHYLATED
    'N10': 'S',  # N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
    'NAL': 'A',  # NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
    'NAM': 'A',  # NAM ALA  NAM NAPTHYLAMINOALANINE
    'NBQ': 'Y',  # NBQ TYR
    'NC1': 'S',  # NC1 SER  NITROCEFIN ACYL-SERINE
    'NCB': 'A',  # NCB ALA  CHEMICAL MODIFICATION
    'NEP': 'H',  # NEP HIS  N1-PHOSPHONOHISTIDINE
    'NFA': 'F',  # NFA PHE  MODIFIED PHENYLALANINE
    'NIY': 'Y',  # NIY TYR  META-NITRO-TYROSINE
    'NLE': 'L',  # NLE LEU  NORLEUCINE
    'NLN': 'L',  # NLN LEU  NORLEUCINE AMIDE
    'NLO': 'L',  # NLO LEU  O-METHYL-L-NORLEUCINE
    'NMC': 'G',  # NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
    'NMM': 'R',  # NMM ARG  MODIFIED ARGININE
    'NPH': 'C',  # NPH CYS
    'NRQ': 'A',  # NRQ ALA
    'NVA': 'V',  # NVA VAL  NORVALINE
    'NYC': 'A',  # NYC ALA
    'NYS': 'C',  # NYS CYS
    'NZH': 'H',  # NZH HIS
    'OAS': 'S',  # OAS SER  O-ACETYLSERINE
    'OBS': 'K',  # OBS LYS  MODIFIED LYSINE
    'OCS': 'C',  # OCS CYS  CYSTEINE SULFONIC ACID
    'OCY': 'C',  # OCY CYS  HYDROXYETHYLCYSTEINE
    'OHI': 'H',  # OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
    'OHS': 'D',  # OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
    'OLT': 'T',  # OLT THR  O-METHYL-L-THREONINE
    'OMT': 'M',  # OMT MET  METHIONINE SULFONE
    'OPR': 'R',  # OPR ARG  C-(3-OXOPROPYL)ARGININE
    'ORN': 'A',  # ORN ALA  ORNITHINE
    'ORQ': 'R',  # ORQ ARG  N~5~-ACETYL-L-ORNITHINE
    'OSE': 'S',  # OSE SER  O-SULFO-L-SERINE
    'OTY': 'Y',  # OTY TYR
    'OXX': 'D',  # OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
    'P1L': 'C',  # P1L CYS  S-PALMITOYL CYSTEINE
    'P2Y': 'P',  # P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
    'PAQ': 'Y',  # PAQ TYR  SEE REMARK 999
    'PAT': 'W',  # PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
    'PBB': 'C',  # PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
    'PBF': 'F',  # PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
    'PCA': 'P',  # PCA PRO  5-OXOPROLINE
    'PCS': 'F',  # PCS PHE  PHENYLALANYLMETHYLCHLORIDE
    'PEC': 'C',  # PEC CYS  S,S-PENTYLTHIOCYSTEINE
    'PF5': 'F',  # PF5 PHE  2,3,4,5,6-PENTAFLUORO-L-PHENYLALANINE
    'PFF': 'F',  # PFF PHE  4-FLUORO-L-PHENYLALANINE
    'PG1': 'S',  # PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
    'PG9': 'G',  # PG9 GLY  D-PHENYLGLYCINE
    'PHA': 'F',  # PHA PHE  PHENYLALANINAL
    'PHD': 'D',  # PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
    'PHI': 'F',  # PHI PHE  IODO-PHENYLALANINE
    'PHL': 'F',  # PHL PHE  L-PHENYLALANINOL
    'PHM': 'F',  # PHM PHE  PHENYLALANYLMETHANE
    'PIA': 'A',  # PIA ALA  FUSION OF ALA 65, TYR 66, GLY 67
    'PLE': 'L',  # PLE LEU  LEUCINE PHOSPHINIC ACID
    'PM3': 'F',  # PM3 PHE
    'POM': 'P',  # POM PRO  CIS-5-METHYL-4-OXOPROLINE
    'PPH': 'L',  # PPH LEU  PHENYLALANINE PHOSPHINIC ACID
    'PPN': 'F',  # PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
    'PR3': 'C',  # PR3 CYS  INE DTT-CYSTEINE
    'PRQ': 'F',  # PRQ PHE  PHENYLALANINE
    'PRR': 'A',  # PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
    'PRS': 'P',  # PRS PRO  THIOPROLINE
    'PSA': 'F',  # PSA PHE
    'PSH': 'H',  # PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
    'PTH': 'Y',  # PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
    'PTM': 'Y',  # PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
    'PTR': 'Y',  # PTR TYR  O-PHOSPHOTYROSINE
    'PYA': 'A',  # PYA ALA  3-(1,10-PHENANTHROL-2-YL)-L-ALANINE
    'PYC': 'A',  # PYC ALA  PYRROLE-2-CARBOXYLATE
    'PYR': 'S',  # PYR SER  CHEMICALLY MODIFIED
    'PYT': 'A',  # PYT ALA  MODIFIED ALANINE
    'PYX': 'C',  # PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
    'R1A': 'C',  # R1A CYS
    'R1B': 'C',  # R1B CYS
    'R1F': 'C',  # R1F CYS
    'R7A': 'C',  # R7A CYS
    'RC7': 'A',  # RC7 ALA
    'RCY': 'C',  # RCY CYS
    'S1H': 'S',  # S1H SER
    'SAC': 'S',  # SAC SER
    'SAH': 'C',  # SAH CYS
    'SAR': 'G',  # SAR GLY
    'SBD': 'S',  # SBD SER
    'SBG': 'S',  # SBG SER
    'SBL': 'S',  # SBL SER
    'SC2': 'C',  # SC2 CYS  N-ACETYL-L-CYSTEINE
    'SCH': 'C',  # SCH CYS  S-METHYL THIOCYSTEINE GROUP
    'SCS': 'C',  # SCS CYS  MODIFIED CYSTEINE
    'SCY': 'C',  # SCY CYS  CETYLATED CYSTEINE
    'SDP': 'S',  # SDP SER
    'SEB': 'S',  # SEB SER  O-BENZYLSULFONYL-SERINE
    'SEC': 'A',  # SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
    'SEL': 'S',  # SEL SER  2-AMINO-1,3-PROPANEDIOL
    'SEP': 'S',  # SEP SER  E PHOSPHOSERINE
    'SET': 'S',  # SET SER  AMINOSERINE
    'SGB': 'S',  # SGB SER  MODIFIED SERINE
    'SGR': 'S',  # SGR SER  MODIFIED SERINE
    'SHC': 'C',  # SHC CYS  S-HEXYLCYSTEINE
    'SHP': 'G',  # SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
    'SIC': 'A',  # SIC ALA
    'SLZ': 'K',  # SLZ LYS  L-THIALYSINE
    'SMC': 'C',  # SMC CYS  POST-TRANSLATIONAL MODIFICATION
    'SME': 'M',  # SME MET  METHIONINE SULFOXIDE
    'SMF': 'F',  # SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
    'SNC': 'C',  # SNC CYS  S-NITROSO CYSTEINE
    'SNN': 'D',  # SNN ASP  POST-TRANSLATIONAL MODIFICATION
    'SOC': 'C',  # SOC CYS  DIOXYSELENOCYSTEINE
    'SOY': 'S',  # SOY SER  OXACILLOYL-ACYLATED SERINE
    'SUI': 'A',  # SUI ALA
    'SUN': 'S',  # SUN SER  TABUN CONJUGATED SERINE
    'SVA': 'S',  # SVA SER  SERINE VANADATE
    'SVV': 'S',  # SVV SER  MODIFIED SERINE
    'SVX': 'S',  # SVX SER  MODIFIED SERINE
    'SVY': 'S',  # SVY SER  MODIFIED SERINE
    'SVZ': 'S',  # SVZ SER  MODIFIED SERINE
    'SXE': 'S',  # SXE SER  MODIFIED SERINE
    'TBG': 'G',  # TBG GLY  T-BUTYL GLYCINE
    'TBM': 'T',  # TBM THR
    'TCQ': 'Y',  # TCQ TYR  MODIFIED TYROSINE
    'TEE': 'C',  # TEE CYS  POST-TRANSLATIONAL MODIFICATION
    'TH5': 'T',  # TH5 THR  O-ACETYL-L-THREONINE
    'THC': 'T',  # THC THR  N-METHYLCARBONYLTHREONINE
    'TIH': 'A',  # TIH ALA  BETA(2-THIENYL)ALANINE
    'TMD': 'T',  # TMD THR  N-METHYLATED, EPSILON C ALKYLATED
    'TNB': 'C',  # TNB CYS  S-(2,3,6-TRINITROPHENYL)CYSTEINE
    'TOX': 'W',  # TOX TRP
    'TPL': 'W',  # TPL TRP  TRYTOPHANOL
    'TPO': 'T',  # TPO THR  HOSPHOTHREONINE
    'TPQ': 'A',  # TPQ ALA  2,4,5-TRIHYDROXYPHENYLALANINE
    'TQQ': 'W',  # TQQ TRP
    'TRF': 'W',  # TRF TRP  N1-FORMYL-TRYPTOPHAN
    'TRN': 'W',  # TRN TRP  AZA-TRYPTOPHAN
    'TRO': 'W',  # TRO TRP  2-HYDROXY-TRYPTOPHAN
    'TRQ': 'W',  # TRQ TRP
    'TRW': 'W',  # TRW TRP
    'TRX': 'W',  # TRX TRP  6-HYDROXYTRYPTOPHAN
    'TTQ': 'W',  # TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
    'TTS': 'Y',  # TTS TYR
    'TY2': 'Y',  # TY2 TYR  3-AMINO-L-TYROSINE
    'TY3': 'Y',  # TY3 TYR  3-HYDROXY-L-TYROSINE
    'TYB': 'Y',  # TYB TYR  TYROSINAL
    'TYC': 'Y',  # TYC TYR  L-TYROSINAMIDE
    'TYI': 'Y',  # TYI TYR  3,5-DIIODOTYROSINE
    'TYN': 'Y',  # TYN TYR  ADDUCT AT HYDROXY GROUP
    'TYO': 'Y',  # TYO TYR
    'TYQ': 'Y',  # TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
    'TYS': 'Y',  # TYS TYR  INE SULPHONATED TYROSINE
    'TYT': 'Y',  # TYT TYR
    'TYX': 'C',  # TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
    'TYY': 'Y',  # TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
    'TYZ': 'R',  # TYZ ARG  PARA ACETAMIDO BENZOIC ACID
    'UMA': 'A',  # UMA ALA
    'VAD': 'V',  # VAD VAL  DEAMINOHYDROXYVALINE
    'VAF': 'V',  # VAF VAL  METHYLVALINE
    'VDL': 'V',  # VDL VAL  (2R,3R)-2,3-DIAMINOBUTANOIC ACID
    'VLL': 'V',  # VLL VAL  (2S)-2,3-DIAMINOBUTANOIC ACID
    'HSD': 'H',  # HSD HIS  ND1-H TAUTOMER OF HISTIDINE
    'VME': 'V',  # VME VAL  O- METHYLVALINE
    'X9Q': 'A',  # X9Q ALA
    'XX1': 'K',  # XX1 LYS
    'XXY': 'A',  # XXY ALA
    'XYG': 'A',  # XYG ALA
    'YCM': 'C',  # YCM CYS
    'YOF': 'Y',  # YOF TYR  3-FLUOROTYROSINE
    'MSE': 'M',  # Selenomethionine -> Methionine
    'SEC': 'C',  # Selenocysteine -> Cysteine  
    'PYL': 'K',  # Pyrrolysine -> Lysine
}

# Default secondary structure for unknown peptides
DEFAULT_PEPTIDE_SS: Final[str] = 'CHTE'

# File extensions and formats
PDB_EXTENSIONS: Final[List[str]] = ['.pdb', '.pdb.gz']
IMAGE_FORMATS: Final[List[str]] = ['svg', 'png', 'pdf', 'eps']
DEFAULT_IMAGE_FORMAT: Final[str] = 'svg'

# CABS-specific constants
CABS_LATTICE_DEFAULTS: Final[Dict[str, Union[float, Tuple[float, float]]]] = {
    'grid_spacing': 0.61,
    'r12': (3.28, 4.27),
    'r13': (4.1, 7.35)
}

# Side chain modeling constants
SC_MODELING_THRESHOLDS: Final[Dict[str, float]] = {
    'min_distance': 5.3,
    'max_distance': 6.4,
    'slope_factor': 1.1
}

# PDB output options
PDB_OUTPUT_OPTIONS: Final[Dict[str, str]] = {
    'R': 'replicas',
    'F': 'filtered',
    'C': 'clusters',
    'M': 'models',
    'S': 'starting',
    'A': 'all',
    'N': 'none'
}

# Beta factor output options
BFAC_OUTPUT_OPTIONS: Final[Dict[str, str]] = {
    'B': 'bfac',
    'C': 'category',
    'P': 'plddt',
    'R': 'rmsf',
    'S': 'secstr',
    'A': 'all',
    'N': 'none'
}

# CSV output options
CSV_OUTPUT_OPTIONS: Final[Dict[str, str]] = {
    'B': 'bfac',
    'C': 'category',
    'P': 'plddt',
    'S': 'secstr',
    'A': 'all',
    'N': 'none'
}

# Valid letters for different output types
VALID_PDB_OUTPUT_LETTERS: Final[str] = 'RFCMSAN'
VALID_BFAC_OUTPUT_LETTERS: Final[str] = 'ABCPRSN'
VALID_CSV_OUTPUT_LETTERS: Final[str] = 'ABCPSN'

# Protein restraints modes
PROTEIN_RESTRAINTS_MODES: Final[List[str]] = ['rigid', 'plddt', 'manual', 'flexible']
PROTEIN_CATEGORY_MODES: Final[List[str]] = ['rigid', 'flexible', 'unleashed']

# Special protein restraints aliases
UNLEASHED_ALIASES: Final[List[str]] = ['none', 'unleashed', 'no-protein-restraints']

# Amino acid reconstruction methods
ALLOWED_AA_METHODS: Final[List[str]] = ['modeller', 'cg2all']

# CABS files that are generated during simulation
CABS_FILES: Final[List[str]] = [
    'TRAF', 'SEQ', 'INP', 'FCHAINS', 'EPAIRMOD'
]

# DSSP secondary structure mapping
DSSP_SS_MAPPING: Final[Dict[str, str]] = {
    'HGIP': 'H',  # Helix
    'BE': 'E',    # Extended/Beta
    'T': 'T',     # Turn
    # Default is 'C' (Coil)
}

# Default colors for contact maps and plots
DEFAULT_COLORS: Final[List[ColorHex]] = [
    '#ffffff',  # White
    '#f2d600',  # Yellow
    '#4b8f24',  # Green
    '#666666',  # Gray
    '#e80915',  # Red
    '#000000'   # Black
]

# Contact map histogram constants
CONTACT_MAP_CONSTANTS: Final[Dict[str, Union[int, float]]] = {
    'max_bars': 50,
    'width_const': 50,
    'default_max_y': 0.05
}

# Peptide replacement patterns
PEPTIDE_REPLACEMENTS: Final[Dict[str, str]] = {
    'PEP ': 'PEP1 '
}

# File paths and directories
OUTPUT_DIRECTORIES: Final[Dict[str, str]] = {
    'pdbs': 'output_pdbs',
    'data': 'output_data',
    'plots': 'plots',
    'contact_maps': 'contact_maps'
}

# Default file names
DEFAULT_FILENAMES: Final[Dict[str, str]] = {
    'log': 'CABS.log',
    'config': 'config.ini',
    'restraints': 'restraints.txt',
    'plddt_config': 'plddt.config'
}

# Calculation thresholds and limits
CALCULATION_CONSTANTS: Final[Dict[str, Union[int, float]]] = {
    'large_value': 1000.0,
    'tiny_value': 0.001,
    'gauss_max_iter': 100,
    'max_dimension_check': 5,
    'progress_width': 65,
    'line_break': 76,
    'default_log_level': 2,
    'zero_threshold': 0,
    'one_threshold': 1
}

# Legacy constants for backward compatibility  
_LARGE: Final[float] = 1000.0  # sort of ...
_TINY: Final[float] = 0.001   # useful only for rmsd/rmsf calc
GAUSS_MAX_ITER: Final[int] = 100

# Mathematical constants
MATH_CONSTANTS: Final[Dict[str, Union[int, float]]] = {
    'multiplier_variant': 2,
    'percent_scale': 100,
    'diagonal_fill': 0
}

# String patterns and templates
STRING_PATTERNS: Final[Dict[str, str]] = {
    'pep_pattern': r'PEP$',
    'pep_replacement': 'PEP1',
    'random_placement': 'random'
}

# System and platform constants
SYSTEM_CONSTANTS: Final[Dict[str, Union[str, List[str]]]] = {
    'windows': 'Windows',
    'dev_null': '/dev/null',
    'stdin': '/dev/stdin',
    'python2_pdb_path': '/usr/lib/python2.7/pdb.py'
}

# Model file extensions and paths
MODEL_CONSTANTS: Final[Dict[str, Union[str, List[str]]]] = {
    'pdb_extension': '.pdb',
    'cbs_extension': '.cbs',
    'gz_extension': '.gz',
    'json_extension': '.json',
    'txt_extension': '.txt',
    'csv_extension': '.csv'
}

# Configuration header template
CONFIG_HEADER: Final[str] = """# CABS configuration file
# Generated automatically - modify with caution
# Lines starting with # are comments and will be ignored

"""

# Error messages and warnings
ERROR_MESSAGES: Final[Dict[str, str]] = {
    'invalid_pdb_output': "Contains letters outside of 'RFCMSAN'.",
    'invalid_bfac_output': "Contains letters outside of 'ABCPRSN'.",
    'invalid_csv_output': "Contains letters outside of 'ABCPSN'.",
    'no_peptide': 'No peptide given',
    'modeller_not_found': 'Modeller not found. Skipping AA rebuild.',
    'unknown_aa_method': 'Unknown AA method: %s. Skipping AA rebuild.'
}

# Default values for various parameters
DEFAULT_VALUES: Final[Dict[str, Union[int, float, str, bool]]] = {
    'iterations': 1,
    'modeller_iterations': 3,
    'flexibility': 1.0,
    'plddt': 1.0,
    'bfac': 0.0,
    'temperature': 300.0,
    'random_seed': None,
    'replicas': 10,
    'mc_cycles': 50,
    'mc_annealing': 20,
    'work_dir': '.',
    'verbosity': 2
}

# NetSurfP-3.0 related constants
NSP3_CONSTANTS: Final[Dict[str, str]] = {
    'model_path_config_key': 'nsp3_model_path',
    'cpu_device': 'cpu',
    'peptide_id': '>peptide'
}

# Legacy constants for backward compatibility
CONFIG_HEADER: Final[str] = """############### CABSdock CONFIGURATION FILE ################

; Options available from the command line can be set here.
; Run CABSdock with -c <config file name> option
;
; Options set from the command line overwrite these set from
; the config file, unless option supports accumulation of
; the arguments. In such case arguments are first accumula-
; ted in order they appear in the config file or on the com-
; mand line. Finally arguments from the command line are ap-
; pended to those from the config file.

########################## SYNTAX ##########################

# this is a comment
; this is also a comment

################### ONE-ARGUMENT OPTIONS ###################

; option = value             OK
; option : value             OK
; option value               NO

################ MULTIPLE ARGUMENT OPTIONS #################

; option = value1 value2     OK
; option : value1 value2     OK
; option = value1, value2    NO

########################## FLAGS ###########################

; flag                       OK
; flag = 1                   NO
; flag = True                NO
; set flag                   NO

############################################################
"""

# For legacy compatibility - deprecated aliases
CABS_SS_reverse = CABS_SS_REVERSE  # Deprecated: use CABS_SS_REVERSE
_CABS_files = CABS_FILES  # Deprecated: use CABS_FILES
_allowed_aa_methods = ALLOWED_AA_METHODS  # Deprecated: use ALLOWED_AA_METHODS
