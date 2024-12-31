import numpy as np
import re
import warnings
from string import ascii_uppercase
from pkg_resources import resource_filename
from random import uniform
from pathlib import Path
import subprocess
from CABS import logger
import os
from tempfile import mkstemp
from contextlib import closing
from biopandas.pdb import PandasPdb
import pandas as pd


# Dictionary for conversion of secondary structure from DSSP to CABS
CABS_SS = {'C': 1, 'H': 2, 'T': 3, 'E': 4, 'c': 1, 'h': 2, 't': 3, 'e': 4}

# sidechains relative coords
SIDECNT = {
    'CYS': (-0.139, -1.265, 1.619, 0.019, -0.813, 1.897),
    'GLN': (-0.095, -1.674, 2.612, 0.047, -0.886, 2.991),
    'ILE': (0.094, -1.416, 1.836, -0.105, -0.659, 2.219),
    'SER': (0.121, -1.476, 1.186, 0.223, -1.042, 1.571),
    'VAL': (0.264, -1.194, 1.531, 0.077, -0.631, 1.854),
    'MET': (-0.04, -1.446, 2.587, 0.072, -0.81, 2.883),
    'PRO': (-0.751, -1.643, 0.467, -1.016, -1.228, 0.977),
    'LYS': (0.032, -1.835, 2.989, 0.002, -0.882, 3.405),
    'THR': (0.075, -1.341, 1.398, 0.051, -0.909, 1.712),
    'PHE': (0.151, -1.256, 3.161, -0.448, -0.791, 3.286),
    'ALA': (0.253, -1.133, 0.985, 0.119, -0.763, 1.312),
    'HIS': (-0.301, -1.405, 2.801, -0.207, -0.879, 3.019),
    'GLY': (0.0, -0.111, -0.111, 0.0, -0.111, -0.111),
    'ASP': (-0.287, -1.451, 1.989, 0.396, -0.798, 2.313),
    'GLU': (-0.028, -1.774, 2.546, 0.096, -0.923, 3.016),
    'LEU': (-0.069, -1.247, 2.292, 0.002, -0.462, 2.579),
    'ARG': (-0.057, -2.522, 3.639, -0.057, -1.21, 3.986),
    'TRP': (0.558, -1.694, 3.433, -0.06, -0.574, 3.834),
    'ASN': (-0.402, -1.237, 2.111, 0.132, -0.863, 2.328),
    'TYR': (0.308, -1.387, 3.492, -0.618, -0.799, 3.634)}

# Library of 1000 random peptides with up to 50 amino acids each
RANDOM_LIGAND_LIBRARY = np.reshape(np.fromfile(
    resource_filename('CABS', 'data/data2.dat'), sep=' '), (1000, 50, 3))

CONFIG_HEADER = """############### CABSdock CONFIGURATION FILE ################

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

AA_NAMES = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

AA_SUB_NAMES = {
    '0CS': 'ALA',  # 0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
    '1AB': 'PRO',  # 1AB PRO  1,4-DIDEOXY-1,4-IMINO-D-ARABINITOL
    '1LU': 'LEU',  # 1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
    '1PA': 'PHE',  # 1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
    '1TQ': 'TRP',  # 1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
    '1TY': 'TYR',  # 1TY TYR
    '23F': 'PHE',  # 23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
    '23S': 'TRP',  # 23S TRP  MODIFIED TRYPTOPHAN
    '2BU': 'ALA',  # 2BU ADE
    '2ML': 'LEU',  # 2ML LEU  2-METHYLLEUCINE
    '2MR': 'ARG',  # 2MR ARG  N3, N4-DIMETHYLARGININE
    '2MT': 'PRO',  # 2MT PRO
    '2OP': 'ALA',  # 2OP (2S  2-HYDROXYPROPANAL
    '2TY': 'TYR',  # 2TY TYR
    '32S': 'TRP',  # 32S TRP  MODIFIED TRYPTOPHAN
    '32T': 'TRP',  # 32T TRP  MODIFIED TRYPTOPHAN
    '3AH': 'HIS',  # 3AH HIS
    '3MD': 'ASP',  # 3MD ASP  2S,3S-3-METHYLASPARTIC ACID
    '3TY': 'TYR',  # 3TY TYR  MODIFIED TYROSINE
    '4DP': 'TRP',  # 4DP TRP
    '4F3': 'ALA',  # 4F3 ALA  CYCLIZED
    '4FB': 'PRO',  # 4FB PRO  (4S)-4-FLUORO-L-PROLINE
    '4FW': 'TRP',  # 4FW TRP  4-FLUOROTRYPTOPHANE
    '4HT': 'TRP',  # 4HT TRP  4-HYDROXYTRYPTOPHAN
    '4IN': 'TRP',  # 4IN TRP  4-AMINO-L-TRYPTOPHAN
    '4PH': 'PHE',  # 4PH PHE  4-METHYL-L-PHENYLALANINE
    '5CS': 'CYS',  # 5CS CYS
    '6CL': 'LYS',  # 6CL LYS  6-CARBOXYLYSINE
    '6CW': 'TRP',  # 6CW TRP  6-CHLORO-L-TRYPTOPHAN
    'A0A': 'ASP',  # A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
    'AA4': 'ALA',  # AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
    'AAR': 'ARG',  # AAR ARG  ARGININEAMIDE
    'AB7': 'GLU',  # AB7 GLU  ALPHA-AMINOBUTYRIC ACID
    'ABA': 'ALA',  # ABA ALA  ALPHA-AMINOBUTYRIC ACID
    'ACB': 'ASP',  # ACB ASP  3-METHYL-ASPARTIC ACID
    'ACL': 'ARG',  # ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
    'ACY': 'GLY',  # ACY GLY  POST-TRANSLATIONAL MODIFICATION
    'AEI': 'THR',  # AEI THR  ACYLATED THR
    'AFA': 'ASN',  # AFA ASN  N-[7-METHYL-OCT-2,4-DIENOYL]ASPARAGINE
    'AGM': 'ARG',  # AGM ARG  4-METHYL-ARGININE
    'AGT': 'CYS',  # AGT CYS  AGMATINE-CYSTEINE ADDUCT
    'AHB': 'ASN',  # AHB ASN  BETA-HYDROXYASPARAGINE
    'AHO': 'ALA',  # AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
    'AHP': 'ALA',  # AHP ALA  2-AMINO-HEPTANOIC ACID
    'AIB': 'ALA',  # AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
    'AKL': 'ASP',  # AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
    'ALA': 'ALA',  # ALA ALA
    'ALC': 'ALA',  # ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
    'ALG': 'ARG',  # ALG ARG  GUANIDINOBUTYRYL GROUP
    'ALM': 'ALA',  # ALM ALA  1-METHYL-ALANINAL
    'ALN': 'ALA',  # ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
    'ALO': 'THR',  # ALO THR  ALLO-THREONINE
    'ALS': 'ALA',  # ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
    'ALT': 'ALA',  # ALT ALA  THIOALANINE
    'ALY': 'LYS',  # ALY LYS  N(6)-ACETYLLYSINE
    'AME': 'MET',  # AME MET  ACETYLATED METHIONINE
    'AP7': 'ALA',  # AP7 ADE
    'APH': 'ALA',  # APH ALA  P-AMIDINOPHENYL-3-ALANINE
    'API': 'LYS',  # API LYS  2,6-DIAMINOPIMELIC ACID
    'APK': 'LYS',  # APK LYS
    'AR2': 'ARG',  # AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
    'AR4': 'GLU',  # AR4 GLU
    'ARG': 'ARG',  # ARG ARG
    'ARM': 'ARG',  # ARM ARG  DEOXY-METHYL-ARGININE
    'ARO': 'ARG',  # ARO ARG  C-GAMMA-HYDROXY ARGININE
    'ASA': 'ASP',  # ASA ASP  ASPARTIC ALDEHYDE
    'ASB': 'ASP',  # ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
    'ASI': 'ASP',  # ASI ASP  L-ISO-ASPARTATE
    'ASK': 'ASP',  # ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
    'ASL': 'ASP',  # ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
    'ASN': 'ASN',  # ASN ASN
    'ASP': 'ASP',  # ASP ASP
    'AYA': 'ALA',  # AYA ALA  N-ACETYLALANINE
    'AYG': 'ALA',  # AYG ALA
    'AZK': 'LYS',  # AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
    'B2A': 'ALA',  # B2A ALA  ALANINE BORONIC ACID
    'B2F': 'PHE',  # B2F PHE  PHENYLALANINE BORONIC ACID
    'B2I': 'ILE',  # B2I ILE  ISOLEUCINE BORONIC ACID
    'B2V': 'VAL',  # B2V VAL  VALINE BORONIC ACID
    'B3A': 'ALA',  # B3A ALA  (3S)-3-AMINOBUTANOIC ACID
    'B3D': 'ASP',  # B3D ASP  3-AMINOPENTANEDIOIC ACID
    'B3E': 'GLU',  # B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
    'B3K': 'LYS',  # B3K LYS  (3S)-3,7-DIAMINOHEPTANOIC ACID
    'B3S': 'SER',  # B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
    'B3X': 'ASN',  # B3X ASN  (3S)-3,5-DIAMINO-5-OXOPENTANOIC ACID
    'B3Y': 'TYR',  # B3Y TYR
    'BAL': 'ALA',  # BAL ALA  BETA-ALANINE
    'BBC': 'CYS',  # BBC CYS
    'BCS': 'CYS',  # BCS CYS  BENZYLCYSTEINE
    'BCX': 'CYS',  # BCX CYS  BETA-3-CYSTEINE
    'BFD': 'ASP',  # BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
    'BG1': 'SER',  # BG1 SER
    'BHD': 'ASP',  # BHD ASP  BETA-HYDROXYASPARTIC ACID
    'BIF': 'PHE',  # BIF PHE
    'BLE': 'LEU',  # BLE LEU  LEUCINE BORONIC ACID
    'BLY': 'LYS',  # BLY LYS  LYSINE BORONIC ACID
    'BMT': 'THR',  # BMT THR
    'BNN': 'ALA',  # BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
    'BOR': 'ARG',  # BOR ARG
    'BPE': 'CYS',  # BPE CYS
    'BTR': 'TRP',  # BTR TRP  6-BROMO-TRYPTOPHAN
    'BUC': 'CYS',  # BUC CYS  S,S-BUTYLTHIOCYSTEINE
    'BUG': 'LEU',  # BUG LEU  TERT-LEUCYL AMINE
    'C12': 'ALA',  # C12 ALA
    'C1X': 'LYS',  # C1X LYS  MODIFIED LYSINE
    'C3Y': 'CYS',  # C3Y CYS  MODIFIED CYSTEINE
    'C5C': 'CYS',  # C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
    'C6C': 'CYS',  # C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
    'C99': 'ALA',  # C99 ALA
    'CAB': 'ALA',  # CAB ALA  4-CARBOXY-4-AMINOBUTANAL
    'CAF': 'CYS',  # CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
    'CAS': 'CYS',  # CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
    'CCS': 'CYS',  # CCS CYS  CARBOXYMETHYLATED CYSTEINE
    'CGU': 'GLU',  # CGU GLU  CARBOXYLATION OF THE CG ATOM
    'CH6': 'ALA',  # CH6 ALA
    'CH7': 'ALA',  # CH7 ALA
    'CHG': 'GLY',  # CHG GLY  CYCLOHEXYL GLYCINE
    'CHP': 'GLY',  # CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
    'CHS': 'PHE',  # CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
    'CIR': 'ARG',  # CIR ARG  CITRULLINE
    'CLB': 'ALA',  # CLB ALA
    'CLD': 'ALA',  # CLD ALA
    'CLE': 'LEU',  # CLE LEU  LEUCINE AMIDE
    'CLG': 'LYS',  # CLG LYS
    'CLH': 'LYS',  # CLH LYS
    'CLV': 'ALA',  # CLV ALA
    'CME': 'CYS',  # CME CYS  MODIFIED CYSTEINE
    'CML': 'CYS',  # CML CYS
    'CMT': 'CYS',  # CMT CYS  O-METHYLCYSTEINE
    'CQR': 'ALA',  # CQR ALA
    'CR2': 'ALA',  # CR2 ALA  POST-TRANSLATIONAL MODIFICATION
    'CR5': 'ALA',  # CR5 ALA
    'CR7': 'ALA',  # CR7 ALA
    'CR8': 'ALA',  # CR8 ALA
    'CRK': 'ALA',  # CRK ALA
    'CRO': 'THR',  # CRO THR  CYCLIZED
    'CRQ': 'TYR',  # CRQ TYR
    'CRW': 'ALA',  # CRW ALA
    'CRX': 'ALA',  # CRX ALA
    'CS1': 'CYS',  # CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
    'CS3': 'CYS',  # CS3 CYS
    'CS4': 'CYS',  # CS4 CYS
    'CSA': 'CYS',  # CSA CYS  S-ACETONYLCYSTEIN
    'CSB': 'CYS',  # CSB CYS  CYS BOUND TO LEAD ION
    'CSD': 'CYS',  # CSD CYS  3-SULFINOALANINE
    'CSE': 'CYS',  # CSE CYS  SELENOCYSTEINE
    'CSI': 'ALA',  # CSI ALA
    'CSO': 'CYS',  # CSO CYS  INE S-HYDROXYCYSTEINE
    'CSR': 'CYS',  # CSR CYS  S-ARSONOCYSTEINE
    'CSS': 'CYS',  # CSS CYS  1,3-THIAZOLE-4-CARBOXYLIC ACID
    'CSU': 'CYS',  # CSU CYS  CYSTEINE-S-SULFONIC ACID
    'CSW': 'CYS',  # CSW CYS  CYSTEINE-S-DIOXIDE
    'CSX': 'CYS',  # CSX CYS  OXOCYSTEINE
    'CSY': 'ALA',  # CSY ALA  MODIFIED TYROSINE COMPLEX
    'CSZ': 'CYS',  # CSZ CYS  S-SELANYL CYSTEINE
    'CTH': 'THR',  # CTH THR  4-CHLOROTHREONINE
    'CWR': 'ALA',  # CWR ALA
    'CXM': 'MET',  # CXM MET  N-CARBOXYMETHIONINE
    'CY0': 'CYS',  # CY0 CYS  MODIFIED CYSTEINE
    'CY1': 'CYS',  # CY1 CYS  ACETAMIDOMETHYLCYSTEINE
    'CY3': 'CYS',  # CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
    'CY4': 'CYS',  # CY4 CYS  S-BUTYRYL-CYSTEIN
    'CY7': 'CYS',  # CY7 CYS  MODIFIED CYSTEINE
    'CYD': 'CYS',  # CYD CYS
    'CYF': 'CYS',  # CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
    'CYG': 'CYS',  # CYG CYS
    'CYJ': 'LYS',  # CYJ LYS  MODIFIED LYSINE
    'CYQ': 'CYS',  # CYQ CYS
    'CYR': 'CYS',  # CYR CYS
    'CYS': 'CYS',  # CYS CYS
    'CZ2': 'CYS',  # CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
    'CZZ': 'CYS',  # CZZ CYS  THIARSAHYDROXY-CYSTEINE
    'DA2': 'ARG',  # DA2 ARG  MODIFIED ARGININE
    'DAB': 'ALA',  # DAB ALA  2,4-DIAMINOBUTYRIC ACID
    'DAH': 'PHE',  # DAH PHE  3,4-DIHYDROXYDAHNYLALANINE
    'DAL': 'ALA',  # DAL ALA  D-ALANINE
    'DAM': 'ALA',  # DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
    'DAR': 'ARG',  # DAR ARG  D-ARGININE
    'DAS': 'ASP',  # DAS ASP  D-ASPARTIC ACID
    'DBU': 'ALA',  # DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
    'DBY': 'TYR',  # DBY TYR  3,5 DIBROMOTYROSINE
    'DBZ': 'ALA',  # DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
    'DCL': 'LEU',  # DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
    'DCY': 'CYS',  # DCY CYS  D-CYSTEINE
    'DDE': 'HIS',  # DDE HIS
    'DGL': 'GLU',  # DGL GLU  D-GLU
    'DGN': 'GLN',  # DGN GLN  D-GLUTAMINE
    'DHA': 'ALA',  # DHA ALA  2-AMINO-ACRYLIC ACID
    'DHI': 'HIS',  # DHI HIS  D-HISTIDINE
    'DHL': 'SER',  # DHL SER  POST-TRANSLATIONAL MODIFICATION
    'DIL': 'ILE',  # DIL ILE  D-ISOLEUCINE
    'DIV': 'VAL',  # DIV VAL  D-ISOVALINE
    'DLE': 'LEU',  # DLE LEU  D-LEUCINE
    'DLS': 'LYS',  # DLS LYS  DI-ACETYL-LYSINE
    'DLY': 'LYS',  # DLY LYS  D-LYSINE
    'DMH': 'ASN',  # DMH ASN  N4,N4-DIMETHYL-ASPARAGINE
    'DMK': 'ASP',  # DMK ASP  DIMETHYL ASPARTIC ACID
    'DNE': 'LEU',  # DNE LEU  D-NORLEUCINE
    'DNG': 'LEU',  # DNG LEU  N-FORMYL-D-NORLEUCINE
    'DNL': 'LYS',  # DNL LYS  6-AMINO-HEXANAL
    'DNM': 'LEU',  # DNM LEU  D-N-METHYL NORLEUCINE
    'DPH': 'PHE',  # DPH PHE  DEAMINO-METHYL-PHENYLALANINE
    'DPL': 'PRO',  # DPL PRO  4-OXOPROLINE
    'DPN': 'PHE',  # DPN PHE  D-CONFIGURATION
    'DPP': 'ALA',  # DPP ALA  DIAMMINOPROPANOIC ACID
    'DPQ': 'TYR',  # DPQ TYR  TYROSINE DERIVATIVE
    'DPR': 'PRO',  # DPR PRO  D-PROLINE
    'DSE': 'SER',  # DSE SER  D-SERINE N-METHYLATED
    'DSG': 'ASN',  # DSG ASN  D-ASPARAGINE
    'DSN': 'SER',  # DSN SER  D-SERINE
    'DTH': 'THR',  # DTH THR  D-THREONINE
    'DTR': 'TRP',  # DTR TRP  D-TRYPTOPHAN
    'DTY': 'TYR',  # DTY TYR  D-TYROSINE
    'DVA': 'VAL',  # DVA VAL  D-VALINE
    'DYG': 'ALA',  # DYG ALA
    'DYS': 'CYS',  # DYS CYS
    'EFC': 'CYS',  # EFC CYS  S,S-(2-FLUOROETHYL)THIOCYSTEINE
    'ESB': 'TYR',  # ESB TYR
    'ESC': 'MET',  # ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
    'FCL': 'PHE',  # FCL PHE  3-CHLORO-L-PHENYLALANINE
    'FGL': 'ALA',  # FGL ALA  2-AMINOPROPANEDIOIC ACID
    'FGP': 'SER',  # FGP SER
    'FHL': 'LYS',  # FHL LYS  MODIFIED LYSINE
    'FLE': 'LEU',  # FLE LEU  FUROYL-LEUCINE
    'FLT': 'TYR',  # FLT TYR  FLUOROMALONYL TYROSINE
    'FME': 'MET',  # FME MET  FORMYL-METHIONINE
    'FOE': 'CYS',  # FOE CYS
    'FOG': 'PHE',  # FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
    'FOR': 'MET',  # FOR MET
    'FRF': 'PHE',  # FRF PHE  PHE FOLLOWED BY REDUCED PHE
    'FTR': 'TRP',  # FTR TRP  FLUOROTRYPTOPHANE
    'FTY': 'TYR',  # FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
    'GHG': 'GLN',  # GHG GLN  GAMMA-HYDROXY-GLUTAMINE
    'GHP': 'GLY',  # GHP GLY  4-HYDROXYPHENYLGLYCINE
    'GL3': 'GLY',  # GL3 GLY  POST-TRANSLATIONAL MODIFICATION
    'GLH': 'GLN',  # GLH GLN
    'GLN': 'GLN',  # GLN GLN
    'GLU': 'GLU',  # GLU GLU
    'GLY': 'GLY',  # GLY GLY
    'GLZ': 'GLY',  # GLZ GLY  AMINO-ACETALDEHYDE
    'GMA': 'GLU',  # GMA GLU  1-AMIDO-GLUTAMIC ACID
    'GMU': 'ALA',  # GMU 5MU
    'GPL': 'LYS',  # GPL LYS  LYSINE GUANOSINE-5'-MONOPHOSPHATE
    'GT9': 'CYS',  # GT9 CYS  SG ALKYLATED
    'GVL': 'SER',  # GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
    'GYC': 'CYS',  # GYC CYS
    'GYS': 'GLY',  # GYS GLY
    'H5M': 'PRO',  # H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
    'HHK': 'ALA',  # HHK ALA  (2S)-2,8-DIAMINOOCTANOIC ACID
    'HIA': 'HIS',  # HIA HIS  L-HISTIDINE AMIDE
    'HIC': 'HIS',  # HIC HIS  4-METHYL-HISTIDINE
    'HIP': 'HIS',  # HIP HIS  ND1-PHOSPHONOHISTIDINE
    'HIQ': 'HIS',  # HIQ HIS  MODIFIED HISTIDINE
    'HIS': 'HIS',  # HIS HIS
    'HLU': 'LEU',  # HLU LEU  BETA-HYDROXYLEUCINE
    'HMF': 'ALA',  # HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
    'HMR': 'ARG',  # HMR ARG  BETA-HOMOARGININE
    'HPE': 'PHE',  # HPE PHE  HOMOPHENYLALANINE
    'HPH': 'PHE',  # HPH PHE  PHENYLALANINOL GROUP
    'HPQ': 'PHE',  # HPQ PHE  HOMOPHENYLALANINYLMETHANE
    'HRG': 'ARG',  # HRG ARG  L-HOMOARGININE
    'HSE': 'SER',  # HSE SER  L-HOMOSERINE
    'HSL': 'SER',  # HSL SER  HOMOSERINE LACTONE
    'HSO': 'HIS',  # HSO HIS  HISTIDINOL
    'HTI': 'CYS',  # HTI CYS
    'HTR': 'TRP',  # HTR TRP  BETA-HYDROXYTRYPTOPHANE
    'HY3': 'PRO',  # HY3 PRO  3-HYDROXYPROLINE
    'HYP': 'PRO',  # HYP PRO  4-HYDROXYPROLINE
    'IAM': 'ALA',  # IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
    'IAS': 'ASP',  # IAS ASP  ASPARTYL GROUP
    'IGL': 'ALA',  # IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
    'IIL': 'ILE',  # IIL ILE  ISO-ISOLEUCINE
    'ILE': 'ILE',  # ILE ILE
    'ILG': 'GLU',  # ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
    'ILX': 'ILE',  # ILX ILE  4,5-DIHYDROXYISOLEUCINE
    'IML': 'ILE',  # IML ILE  N-METHYLATED
    'IPG': 'GLY',  # IPG GLY  N-ISOPROPYL GLYCINE
    'IT1': 'LYS',  # IT1 LYS
    'IYR': 'TYR',  # IYR TYR  3-IODO-TYROSINE
    'KCX': 'LYS',  # KCX LYS  CARBAMOYLATED LYSINE
    'KGC': 'LYS',  # KGC LYS
    'KOR': 'CYS',  # KOR CYS  MODIFIED CYSTEINE
    'KST': 'LYS',  # KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
    'KYN': 'ALA',  # KYN ALA  KYNURENINE
    'LA2': 'LYS',  # LA2 LYS
    'LAL': 'ALA',  # LAL ALA  N,N-DIMETHYL-L-ALANINE
    'LCK': 'LYS',  # LCK LYS
    'LCX': 'LYS',  # LCX LYS  CARBAMYLATED LYSINE
    'LDH': 'LYS',  # LDH LYS  N~6~-ETHYL-L-LYSINE
    'LED': 'LEU',  # LED LEU  POST-TRANSLATIONAL MODIFICATION
    'LEF': 'LEU',  # LEF LEU  2-5-FLUOROLEUCINE
    'LET': 'LYS',  # LET LYS  ODIFIED LYSINE
    'LEU': 'LEU',  # LEU LEU
    'LLP': 'LYS',  # LLP LYS
    'LLY': 'LYS',  # LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
    'LME': 'GLU',  # LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
    'LNT': 'LEU',  # LNT LEU
    'LPD': 'PRO',  # LPD PRO  L-PROLINAMIDE
    'LSO': 'LYS',  # LSO LYS  MODIFIED LYSINE
    'LYM': 'LYS',  # LYM LYS  DEOXY-METHYL-LYSINE
    'LYN': 'LYS',  # LYN LYS  2,6-DIAMINO-HEXANOIC ACID AMIDE
    'LYP': 'LYS',  # LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
    'LYR': 'LYS',  # LYR LYS  MODIFIED LYSINE
    'LYS': 'LYS',  # LYS LYS
    'LYX': 'LYS',  # LYX LYS  N''-(2-COENZYME A)-PROPANOYL-LYSINE
    'LYZ': 'LYS',  # LYZ LYS  5-HYDROXYLYSINE
    'M0H': 'CYS',  # M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
    'M2L': 'LYS',  # M2L LYS
    'M3L': 'LYS',  # M3L LYS  N-TRIMETHYLLYSINE
    'MAA': 'ALA',  # MAA ALA  N-METHYLALANINE
    'MAI': 'ARG',  # MAI ARG  DEOXO-METHYLARGININE
    'MBQ': 'TYR',  # MBQ TYR
    'MC1': 'SER',  # MC1 SER  METHICILLIN ACYL-SERINE
    'MCL': 'LYS',  # MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
    'MCS': 'CYS',  # MCS CYS  MALONYLCYSTEINE
    'MDO': 'ALA',  # MDO ALA
    'MEA': 'PHE',  # MEA PHE  N-METHYLPHENYLALANINE
    'MEG': 'GLU',  # MEG GLU  (2S,3R)-3-METHYL-GLUTAMIC ACID
    'MEN': 'ASN',  # MEN ASN  GAMMA METHYL ASPARAGINE
    'MET': 'MET',  # MET MET
    'MEU': 'GLY',  # MEU GLY  O-METHYL-GLYCINE
    'MFC': 'ALA',  # MFC ALA  CYCLIZED
    'MGG': 'ARG',  # MGG ARG  MODIFIED D-ARGININE
    'MGN': 'GLN',  # MGN GLN  2-METHYL-GLUTAMINE
    'MHL': 'LEU',  # MHL LEU  N-METHYLATED, HYDROXY
    'MHO': 'MET',  # MHO MET  POST-TRANSLATIONAL MODIFICATION
    'MHS': 'HIS',  # MHS HIS  1-N-METHYLHISTIDINE
    'MIS': 'SER',  # MIS SER  MODIFIED SERINE
    'MLE': 'LEU',  # MLE LEU  N-METHYLATED
    'MLL': 'LEU',  # MLL LEU  METHYL L-LEUCINATE
    'MLY': 'LYS',  # MLY LYS  METHYLATED LYSINE
    'MLZ': 'LYS',  # MLZ LYS  N-METHYL-LYSINE
    'MME': 'MET',  # MME MET  N-METHYL METHIONINE
    'MNL': 'LEU',  # MNL LEU  4,N-DIMETHYLNORLEUCINE
    'MNV': 'VAL',  # MNV VAL  N-METHYL-C-AMINO VALINE
    'MPQ': 'GLY',  # MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
    'MSA': 'GLY',  # MSA GLY  (2-S-METHYL) SARCOSINE
    'MSE': 'MET',  # MSE MET  ELENOMETHIONINE
    'MSO': 'MET',  # MSO MET  METHIONINE SULFOXIDE
    'MTY': 'PHE',  # MTY PHE  3-HYDROXYPHENYLALANINE
    'MVA': 'VAL',  # MVA VAL  N-METHYLATED
    'N10': 'SER',  # N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
    'NAL': 'ALA',  # NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
    'NAM': 'ALA',  # NAM ALA  NAM NAPTHYLAMINOALANINE
    'NBQ': 'TYR',  # NBQ TYR
    'NC1': 'SER',  # NC1 SER  NITROCEFIN ACYL-SERINE
    'NCB': 'ALA',  # NCB ALA  CHEMICAL MODIFICATION
    'NEP': 'HIS',  # NEP HIS  N1-PHOSPHONOHISTIDINE
    'NFA': 'PHE',  # NFA PHE  MODIFIED PHENYLALANINE
    'NIY': 'TYR',  # NIY TYR  META-NITRO-TYROSINE
    'NLE': 'LEU',  # NLE LEU  NORLEUCINE
    'NLN': 'LEU',  # NLN LEU  NORLEUCINE AMIDE
    'NLO': 'LEU',  # NLO LEU  O-METHYL-L-NORLEUCINE
    'NMC': 'GLY',  # NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
    'NMM': 'ARG',  # NMM ARG  MODIFIED ARGININE
    'NPH': 'CYS',  # NPH CYS
    'NRQ': 'ALA',  # NRQ ALA
    'NVA': 'VAL',  # NVA VAL  NORVALINE
    'NYC': 'ALA',  # NYC ALA
    'NYS': 'CYS',  # NYS CYS
    'NZH': 'HIS',  # NZH HIS
    'OAS': 'SER',  # OAS SER  O-ACETYLSERINE
    'OBS': 'LYS',  # OBS LYS  MODIFIED LYSINE
    'OCS': 'CYS',  # OCS CYS  CYSTEINE SULFONIC ACID
    'OCY': 'CYS',  # OCY CYS  HYDROXYETHYLCYSTEINE
    'OHI': 'HIS',  # OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
    'OHS': 'ASP',  # OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
    'OLT': 'THR',  # OLT THR  O-METHYL-L-THREONINE
    'OMT': 'MET',  # OMT MET  METHIONINE SULFONE
    'OPR': 'ARG',  # OPR ARG  C-(3-OXOPROPYL)ARGININE
    'ORN': 'ALA',  # ORN ALA  ORNITHINE
    'ORQ': 'ARG',  # ORQ ARG  N~5~-ACETYL-L-ORNITHINE
    'OSE': 'SER',  # OSE SER  O-SULFO-L-SERINE
    'OTY': 'TYR',  # OTY TYR
    'OXX': 'ASP',  # OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
    'P1L': 'CYS',  # P1L CYS  S-PALMITOYL CYSTEINE
    'P2Y': 'PRO',  # P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
    'PAQ': 'TYR',  # PAQ TYR  SEE REMARK 999
    'PAT': 'TRP',  # PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
    'PBB': 'CYS',  # PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
    'PBF': 'PHE',  # PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
    'PCA': 'PRO',  # PCA PRO  5-OXOPROLINE
    'PCS': 'PHE',  # PCS PHE  PHENYLALANYLMETHYLCHLORIDE
    'PEC': 'CYS',  # PEC CYS  S,S-PENTYLTHIOCYSTEINE
    'PF5': 'PHE',  # PF5 PHE  2,3,4,5,6-PENTAFLUORO-L-PHENYLALANINE
    'PFF': 'PHE',  # PFF PHE  4-FLUORO-L-PHENYLALANINE
    'PG1': 'SER',  # PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
    'PG9': 'GLY',  # PG9 GLY  D-PHENYLGLYCINE
    'PHA': 'PHE',  # PHA PHE  PHENYLALANINAL
    'PHD': 'ASP',  # PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
    'PHE': 'PHE',  # PHE PHE
    'PHI': 'PHE',  # PHI PHE  IODO-PHENYLALANINE
    'PHL': 'PHE',  # PHL PHE  L-PHENYLALANINOL
    'PHM': 'PHE',  # PHM PHE  PHENYLALANYLMETHANE
    'PIA': 'ALA',  # PIA ALA  FUSION OF ALA 65, TYR 66, GLY 67
    'PLE': 'LEU',  # PLE LEU  LEUCINE PHOSPHINIC ACID
    'PM3': 'PHE',  # PM3 PHE
    'POM': 'PRO',  # POM PRO  CIS-5-METHYL-4-OXOPROLINE
    'PPH': 'LEU',  # PPH LEU  PHENYLALANINE PHOSPHINIC ACID
    'PPN': 'PHE',  # PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
    'PR3': 'CYS',  # PR3 CYS  INE DTT-CYSTEINE
    'PRO': 'PRO',  # PRO PRO
    'PRQ': 'PHE',  # PRQ PHE  PHENYLALANINE
    'PRR': 'ALA',  # PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
    'PRS': 'PRO',  # PRS PRO  THIOPROLINE
    'PSA': 'PHE',  # PSA PHE
    'PSH': 'HIS',  # PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
    'PTH': 'TYR',  # PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
    'PTM': 'TYR',  # PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
    'PTR': 'TYR',  # PTR TYR  O-PHOSPHOTYROSINE
    'PYA': 'ALA',  # PYA ALA  3-(1,10-PHENANTHROL-2-YL)-L-ALANINE
    'PYC': 'ALA',  # PYC ALA  PYRROLE-2-CARBOXYLATE
    'PYR': 'SER',  # PYR SER  CHEMICALLY MODIFIED
    'PYT': 'ALA',  # PYT ALA  MODIFIED ALANINE
    'PYX': 'CYS',  # PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
    'R1A': 'CYS',  # R1A CYS
    'R1B': 'CYS',  # R1B CYS
    'R1F': 'CYS',  # R1F CYS
    'R7A': 'CYS',  # R7A CYS
    'RC7': 'ALA',  # RC7 ALA
    'RCY': 'CYS',  # RCY CYS
    'S1H': 'SER',  # S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
    'SAC': 'SER',  # SAC SER  N-ACETYL-SERINE
    'SAH': 'CYS',  # SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
    'SAR': 'GLY',  # SAR GLY  SARCOSINE
    'SBD': 'SER',  # SBD SER
    'SBG': 'SER',  # SBG SER  MODIFIED SERINE
    'SBL': 'SER',  # SBL SER
    'SC2': 'CYS',  # SC2 CYS  N-ACETYL-L-CYSTEINE
    'SCH': 'CYS',  # SCH CYS  S-METHYL THIOCYSTEINE GROUP
    'SCS': 'CYS',  # SCS CYS  MODIFIED CYSTEINE
    'SCY': 'CYS',  # SCY CYS  CETYLATED CYSTEINE
    'SDP': 'SER',  # SDP SER
    'SEB': 'SER',  # SEB SER  O-BENZYLSULFONYL-SERINE
    'SEC': 'ALA',  # SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
    'SEL': 'SER',  # SEL SER  2-AMINO-1,3-PROPANEDIOL
    'SEP': 'SER',  # SEP SER  E PHOSPHOSERINE
    'SER': 'SER',  # SER SER
    'SET': 'SER',  # SET SER  AMINOSERINE
    'SGB': 'SER',  # SGB SER  MODIFIED SERINE
    'SGR': 'SER',  # SGR SER  MODIFIED SERINE
    'SHC': 'CYS',  # SHC CYS  S-HEXYLCYSTEINE
    'SHP': 'GLY',  # SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
    'SIC': 'ALA',  # SIC ALA
    'SLZ': 'LYS',  # SLZ LYS  L-THIALYSINE
    'SMC': 'CYS',  # SMC CYS  POST-TRANSLATIONAL MODIFICATION
    'SME': 'MET',  # SME MET  METHIONINE SULFOXIDE
    'SMF': 'PHE',  # SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
    'SNC': 'CYS',  # SNC CYS  S-NITROSO CYSTEINE
    'SNN': 'ASP',  # SNN ASP  POST-TRANSLATIONAL MODIFICATION
    'SOC': 'CYS',  # SOC CYS  DIOXYSELENOCYSTEINE
    'SOY': 'SER',  # SOY SER  OXACILLOYL-ACYLATED SERINE
    'SUI': 'ALA',  # SUI ALA
    'SUN': 'SER',  # SUN SER  TABUN CONJUGATED SERINE
    'SVA': 'SER',  # SVA SER  SERINE VANADATE
    'SVV': 'SER',  # SVV SER  MODIFIED SERINE
    'SVX': 'SER',  # SVX SER  MODIFIED SERINE
    'SVY': 'SER',  # SVY SER  MODIFIED SERINE
    'SVZ': 'SER',  # SVZ SER  MODIFIED SERINE
    'SXE': 'SER',  # SXE SER  MODIFIED SERINE
    'TBG': 'GLY',  # TBG GLY  T-BUTYL GLYCINE
    'TBM': 'THR',  # TBM THR
    'TCQ': 'TYR',  # TCQ TYR  MODIFIED TYROSINE
    'TEE': 'CYS',  # TEE CYS  POST-TRANSLATIONAL MODIFICATION
    'TH5': 'THR',  # TH5 THR  O-ACETYL-L-THREONINE
    'THC': 'THR',  # THC THR  N-METHYLCARBONYLTHREONINE
    'THR': 'THR',  # THR THR
    'TIH': 'ALA',  # TIH ALA  BETA(2-THIENYL)ALANINE
    'TMD': 'THR',  # TMD THR  N-METHYLATED, EPSILON C ALKYLATED
    'TNB': 'CYS',  # TNB CYS  S-(2,3,6-TRINITROPHENYL)CYSTEINE
    'TOX': 'TRP',  # TOX TRP
    'TPL': 'TRP',  # TPL TRP  TRYTOPHANOL
    'TPO': 'THR',  # TPO THR  HOSPHOTHREONINE
    'TPQ': 'ALA',  # TPQ ALA  2,4,5-TRIHYDROXYPHENYLALANINE
    'TQQ': 'TRP',  # TQQ TRP
    'TRF': 'TRP',  # TRF TRP  N1-FORMYL-TRYPTOPHAN
    'TRN': 'TRP',  # TRN TRP  AZA-TRYPTOPHAN
    'TRO': 'TRP',  # TRO TRP  2-HYDROXY-TRYPTOPHAN
    'TRP': 'TRP',  # TRP TRP
    'TRQ': 'TRP',  # TRQ TRP
    'TRW': 'TRP',  # TRW TRP
    'TRX': 'TRP',  # TRX TRP  6-HYDROXYTRYPTOPHAN
    'TTQ': 'TRP',  # TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
    'TTS': 'TYR',  # TTS TYR
    'TY2': 'TYR',  # TY2 TYR  3-AMINO-L-TYROSINE
    'TY3': 'TYR',  # TY3 TYR  3-HYDROXY-L-TYROSINE
    'TYB': 'TYR',  # TYB TYR  TYROSINAL
    'TYC': 'TYR',  # TYC TYR  L-TYROSINAMIDE
    'TYI': 'TYR',  # TYI TYR  3,5-DIIODOTYROSINE
    'TYN': 'TYR',  # TYN TYR  ADDUCT AT HYDROXY GROUP
    'TYO': 'TYR',  # TYO TYR
    'TYQ': 'TYR',  # TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
    'TYR': 'TYR',  # TYR TYR
    'TYS': 'TYR',  # TYS TYR  INE SULPHONATED TYROSINE
    'TYT': 'TYR',  # TYT TYR
    'TYX': 'CYS',  # TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
    'TYY': 'TYR',  # TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
    'TYZ': 'ARG',  # TYZ ARG  PARA ACETAMIDO BENZOIC ACID
    'UMA': 'ALA',  # UMA ALA
    'VAD': 'VAL',  # VAD VAL  DEAMINOHYDROXYVALINE
    'VAF': 'VAL',  # VAF VAL  METHYLVALINE
    'VAL': 'VAL',  # VAL VAL
    'VDL': 'VAL',  # VDL VAL  (2R,3R)-2,3-DIAMINOBUTANOIC ACID
    'VLL': 'VAL',  # VLL VAL  (2S)-2,3-DIAMINOBUTANOIC ACID
    'HSD': 'HIS',  # up his
    'VME': 'VAL',  # VME VAL  O- METHYLVALINE
    'X9Q': 'ALA',  # X9Q ALA
    'XX1': 'LYS',  # XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
    'XXY': 'ALA',  # XXY ALA
    'XYG': 'ALA',  # XYG ALA
    'YCM': 'CYS',  # YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
    'YOF': 'TYR'}  # YOF TYR  3-FLUOROTYROSINE


class SCModeler(object):

    """Side Chain Modeler.

    Rebuilds side chains in CABS representation on C-alpha trace vector.
    """

    def __init__(self, nms):
        """Side Chain Modeler initialization.

        Argument:
        nms -- sequence of CABS.Atom representing subsequent mers.
        """
        self.nms = nms

    @staticmethod
    def _mk_local_system(c1, c2, c3):
        """Return local system base for given CA.

        Arguments:
        c1, c2, c3 -- subsequent CA position vectors.

        Base will be calculated for c2.
        Returns 3 vectors and distance between c1 nad c3.
        """
        rdif = c3 - c1
        rdnorm = np.linalg.norm(rdif)
        rsum = (c3 - c2) + (c1 - c2)
        z = -1 * rsum / np.linalg.norm(rsum)
        x = rdif / rdnorm
        y = np.cross(z, x)
        return x, y, z, rdnorm

    @staticmethod
    def _calc_nodes_line(old_v, new_v):
        """Return versor unchanged during transformation of one versor into another."""
        w = np.cross(old_v, new_v)
        return w / np.linalg.norm(w)

    @staticmethod
    def _calc_trig_fnc(ve1, ve2, axis):
        """Return cos and sin between given versors."""
        cos = np.dot(ve1, ve2)
        perp_vec = np.cross(ve1, ve2)
        sin = np.linalg.norm(perp_vec) * np.sign(np.dot(perp_vec, axis))
        return cos, sin

    def _calc_rot_mtx(self, c1, c2, c3, dbg=False):
        """Return rotation matrix transforming Cartesian system to system
        of alpha carbon c2 in sequence of alpha carbons c1, c2, c3.

        Arguments:
        c1, c2, c3 -- position vectors of subsequent alpha carbons.

        Returns matrix and distance detween c1 and c3.
        """
        if dbg:
            import imp
            pdbx = imp.load_source('test', '/usr/lib/python2.7/pdb.py')
            pdbx.set_trace()
        x, y, z, rdnorm = self._mk_local_system(c1, c2, c3)

        setting = np.geterr()
        np.seterr(all='raise')
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
                [-1 * sps * cph - sph * cps * cth, -1 * sps * sph + cps * cth * cph, cps * sth],
                [sth * sph, -sth * cph, cth]
            ]
        )

        return rot, rdnorm

    @staticmethod
    def _calc_scatter_coef(dist):
        if dist < 5.3:
            return 1.
        if dist > 6.4:
            return 0.
        return float((dist - 5.3) * -(1 / 1.1) + 1)

    def rebuild_one(self, vec, sc=True):
        """Takes vector of C alpha coords and residue names and returns vector of C beta coords."""
        vec = np.insert(vec, 0, vec[0] - (vec[2] - vec[1]), axis=0)
        vec = np.append(vec, np.array([vec[-1] + (vec[-2] - vec[-3])]), axis=0)

        nvec = np.zeros((len(vec) - 2, 3))
        nms = (lambda x: self.nms[x].resname) if sc else (lambda x: 'ALA')

        for i in range(len(vec) - 2):
            rot, casdist = self._calc_rot_mtx(*vec[i: i + 3])
            coef = self._calc_scatter_coef(casdist)
            comp = np.array(SIDECNT[nms(i)][:3]) * coef
            scat = np.array(SIDECNT[nms(i)][3:]) * (1 - coef)
            rbld = np.array(comp + scat)
            nvec[i] = rbld.dot(rot).A1 + vec[i + 1]

        return nvec

    def _calculate_traj(self, traj, sc=False):
        return np.array([np.array([self.rebuild_one(j, sc) for j in i]) for i in traj])

    def calculate_cb_traj(self, traj):
        return self._calculate_traj(traj, False)

    def calculate_sc_traj(self, traj):
        return self._calculate_traj(traj, True)


class InvalidAAName(Exception):
    """Exception raised when invalid amino acid name is used"""

    def __init__(self, name, l):
        self.name = (name, l)

    def __str__(self):
        return '%s is not a valid %d-letter amino acid code' % self.name


def aa_to_long(seq):
    """Converts short amino acid name to long."""
    s = seq.upper()
    if s in AA_NAMES:
        return AA_NAMES[s]
    else:
        raise InvalidAAName(seq, 1)


def aa_to_short(seq):
    """Converts long amino acid name to short."""
    s = seq.upper()
    for short, full in AA_NAMES.items():
        if full == s:
            return short
    else:
        raise InvalidAAName(seq, 3)


def next_letter(taken_letters):
    """Returns next available letter for new protein chain."""
    return re.sub('[' + taken_letters + ']', '', ascii_uppercase)[0]


def line_count(filename):
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f, 1):
            pass
    return i


def ranges(data):
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


def kabsch(target, query, weights=None, concentric=False):
    """
    Function for the calculation of the best fit rotation.

    target     - a N x 3 np.array with coordinates of the reference structure
    query      - a N x 3 np.array with coordinates of the fitted structure
    weights    - a N-length list with weights - floats from [0:1]
    concentric - True/False specifying if target and query are centered at origin

    IMPORTANT: If weights are not None centering at origin should account for them.
    proper procedure: A -= np.average(A, 0, WEIGHTS)

    returns rotation matrix as 3 x 3 np.array
    """

    if not concentric:
        t = target - np.average(target, axis=0, weights=weights)
        q = query - np.average(query, axis=0, weights=weights)
    else:
        t = target
        q = query

    c = np.dot(weights * t.T, q) if weights else np.dot(t.T, q)
    v, s, w = np.linalg.svd(c)
    d = np.identity(3)
    if np.linalg.det(c) < 0:
        d[2, 2] = -1

    return np.dot(np.dot(w.T, d), v.T)


_LARGE = 1000.  # sort of ...
_TINY = 0.001   # useful only for rmsd/rmsf calc


def rmsd(target, query=None, weights=None):
    _diff = target if query is None else query - target
    _rmsd = np.sqrt(np.average(np.sum(_diff ** 2, axis=1), axis=0, weights=weights))
    return _rmsd if _rmsd > _TINY else 0.


GAUSS_MAX_ITER = 100


def dynamic_kabsch(target, query):
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
        w = np.exp(-np.sum(_diff ** 2, axis=1) / max(_rmsd, 2.)).tolist()
    else:
        raise Exception('Dynamic Kabsch did not converge in %i steps.' % GAUSS_MAX_ITER)


def smart_flatten(l):
    """
    Function which expands and flattens a list of integers.
    m-n -> m, m+1, ..., n
    """
    fl = []
    for i in l:
        if '-' in i:
            j = i.split('-')
            if len(j) != 2:
                raise Exception('Invalid range syntax: ' + l)

            # temporary fix, nie bedzie wiecej bledow z tego powodu,
            # ale nie do konca to dziala tak jak powinno
            beg = int(re.sub(r'[^0-9]', '', j[0]))
            end = int(re.sub(r'[^0-9]', '', j[1]))

            if beg > end:
                raise Exception('The left index(%i) is greater than the right(%i)' % (beg, end))
            for k in range(beg, end + 1):
                fl.append(k)
        else:
            fl.append(int(i))
    return fl


def check_peptide_sequence(sequence):
    """
    Checks the peptide sequence for non-standard AAs.
    :param sequence: string the peptide sequence.
    :return: True is the sequence does not contain non-standard AAs. Raises error if does.
    """
    standard_one_letter_residues = AA_NAMES.keys()
    for residue in sequence:
        if residue not in standard_one_letter_residues:
            raise Exception(
                "The input peptide sequence contains a non-standard residue \"{0}\" that is currently not supported."
                "The simulation cannot be performed.".format(residue)
            )
    return True


def fix_residue(residue):
    """
    Fixes non-standard AA residues in the receptor.
    :param residue: string three-letter residue code
    :return:    if the residue is non-standard the method returns three-letter code of the appropriate substitution.
                raises exception if the residue is non-standard and there is no substitution available.
    """
    standard_three_letter_residues = AA_NAMES.values()
    known_non_standard_three_letter_residues = AA_SUB_NAMES.keys()
    if residue in standard_three_letter_residues:
        return residue
    elif residue in known_non_standard_three_letter_residues:
        modified = AA_SUB_NAMES[residue]
        warnings.warn(
            "In the current version residue \"{0}\" is not supported."
            "\"{0}\" was replaced with \"{1}\" to perform the simulation.".format(residue, modified), UserWarning
        )
        return modified
    else:
        raise Exception("The PDB file contains unknown residue \"{0}\"".format(residue))


def _chunk_lst(lst, sl_len, extend_last=None):
    """ Slices given list for slices of given len.

    Arguments:
    lst -- list to be sliced.
    sl_len -- len of one slice.
    extend_last -- value to be put in last slice in order to extend it to proper length.
    """
    slists = []
    while lst:
        slists.append(lst[:sl_len])
        lst = lst[sl_len:]
    if extend_last is not None:
        _extend_last(slists, sl_len, extend_last)
    return slists


def _extend_last(sseries, slen, token):
    try:
        sseries[-1].extend([token] * (slen - len(sseries[-1])))
    except IndexError:
        sseries.append([token] * slen)


def _fmt_res_name(atom):
    return (atom.chid + str(atom.resnum) + atom.icode).strip()


def pep2pep1(_id):
    if re.search('PEP$', _id):
        return _id + '1'
    else:
        return _id

class EprGenerator(object):
    """ class that generates EPR coefficients for peptide residues based on external output (e.g. SASA) """
    def __init__(self, file = None):
        super(EprGenerator, self).__init__()
        self.external_file = file
        self.variant_dictionary = {'easy': self.easy, 'normalized': self.normalized, 'multiplied': self.multiplied}

    def load_file(self, file = ''):
        self.residues = []
        self.external_weights = []
        self.external_parameters = []
        with open(file, 'r') as file:
            for line in file:
                residue, external_weight, external_parameter = line.split()
                self.residues.append(residue)
                self.external_weights.append(float(external_weight))
                self.external_parameters.append(float(external_parameter))

    def generate_epr(self, variant = 'easy', outfile = 'outfile.txt'):
        self.outfile = outfile
        self.variant_dictionary[variant]()

    def save_weights_file(self, outfile = '', new_weights = []):
        with open(self.outfile, 'w') as outfile:
            for residue, new_weight, external_parameter in zip(self.residues, new_weights, self.external_parameters):
                outfile.write('{}\t{}\t{}\n'.format(residue, new_weight, external_parameter))

    def easy(self):
        new_weights = self.external_weights
        self.save_weights_file(outfile=self.outfile, new_weights=new_weights)

    def normalized(self):
        new_weights = len(self.external_weights) * np.array(self.external_weights) / sum(self.external_weights)
        self.save_weights_file(outfile=self.outfile, new_weights=new_weights)

    def multiplied(self):
        new_weights = 2 * np.array(self.external_weights)
        self.save_weights_file(outfile=self.outfile, new_weights=new_weights)


def random_rotation_matrix():
    cos = np.cos
    sin = np.sin

    x = np.random.uniform(0, np.pi * 2)
    y = np.random.uniform(0, np.pi * 2)
    z = np.random.uniform(0, np.pi * 2)

    rx = np.array([[1, 0, 0], [0, cos(x), -sin(x)], [0, sin(x), cos(x)]])
    ry = np.array([[cos(y), 0, sin(y)], [0, 1, 0], [-sin(y), 0, cos(y)]])
    rz = np.array([[cos(z), -sin(z), 0], [sin(z), cos(z), 0], [0, 0, 1]])

    return (rx.dot(ry)).dot(rz)


def convert_cg_to_all(filename, work_dir='.', iter=0, reference_pdb=None,
                      renumber_flag=False):
    """
    Convert coarse-grained model to all-atom
    """
    pdb = mkstemp(prefix='.', suffix='.pdb', dir=work_dir, text=True)[1]

    atoms = []
    pattern = re.compile('ATOM.{9}CA .([A-Z]{3}) ([A-Z ])(.{5}).{27}(.{12}).*')
    with closing(filename) as f, open(pdb, 'w') as tmp:
        for line in f:
            if line.startswith('ENDMDL'):
                break
            else:
                match = re.match(pattern, line)
                if match:
                    atoms.append(match.groups())
                    tmp.write(line)
    if renumber_flag:
        sync_residues(input_pdb_path=reference_pdb,
                      output_pdb_path=tmp.name)

    output_dir = Path(work_dir) / "output_pdbs"
    input_pdb = Path(tmp.name)
    fout = f"model_{iter}.pdb"
    command = f"convert_cg2all -p {input_pdb} -o {output_dir / fout}"
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print("Standard output:")
        print(e.stdout)
        print("Standard error:")
        print(e.stderr)
        logger.critical(module_name="CG2ALL", msg=f'CG2ALL failed with exit code: {e.returncode} and error: {e.stderr}')
        raise
    except Exception as e:
        print(f"An error occurred: {e}")
        logger.warning(module_name="CG2ALL", msg=f'CG2ALL failed with error: {e}')
    finally:
        os.remove(tmp.name)


def sync_residues(input_pdb_path: Path, output_pdb_path: Path):
    input_pdb = PandasPdb().read_pdb(input_pdb_path)
    output_pdb = PandasPdb().read_pdb(output_pdb_path)
    input_atom_df = input_pdb.df['ATOM']
    hetatm_df = input_pdb.df['HETATM']
    hetatm_df['residue_name'] = hetatm_df['residue_name'].map(AA_SUB_NAMES)
    input_atom_df = pd.concat([input_atom_df,
                               hetatm_df]).dropna(subset=['residue_name'])
    output_atom_df = output_pdb.df['ATOM']
    output_atom_df.loc[output_atom_df['atom_name'] == 'CA',
                       'residue_number'] = input_atom_df.loc[
                           input_atom_df['atom_name'] == 'CA',
                           'residue_number'].values
    output_pdb.df['ATOM'] = output_atom_df
    output_pdb.to_pdb(output_pdb_path)
    return "Residues synchronized"
