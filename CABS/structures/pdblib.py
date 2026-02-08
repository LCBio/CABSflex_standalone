"""
Modernized PDB and mmCIF handler for CABS-flex.

This module provides a pure-Python interface for loading, fetching, and
analyzing protein structures. It removes dependencies on external DSSP
binaries by using MDTraj and supports the modern mmCIF format required
for large structures and AlphaFold models.
"""

import os
import gzip
import string
import requests as req
from typing import Dict, List, Optional, Tuple

# Biopython for robust structural parsing
from Bio.PDB import PDBParser, MMCIFParser

# External dependencies for secondary structure assignment
try:
    import mdtraj as md
except ImportError:
    md = None

# CABS internal imports
from CABS.config_loader import get_config_section
from CABS.constants import AA_NAMES, AA_SUB_NAMES
from CABS.io import logger
from CABS.structures.atom import Atom, Atoms
from CABS.structures.vector3d import Vector3d

_name = "PDB"


class Pdb:
    """
    Unified parser and fetcher for PDB and mmCIF structures.

    Attributes:
        atoms (Atoms): Internal CABS Atoms collection.
        source_file (str): Path to the local file used for loading.
    """

    class InvalidPdbInput(Exception):
        """Raised when a structure cannot be retrieved or parsed."""
        pass

    def __init__(
        self,
        source: str,
        selection: str = "",
        pdb_cache: str = "~",
        remove_alternative_locations: bool = True,
        fix_non_standard_aa: bool = True,
        remove_water: bool = True,
        remove_hetero: bool = True,
        verify: bool = False,
        no_exit: bool = False,
        create_from_aa: bool = False,
    ) -> None:
        """
        Loads a structure from a file path or the RCSB database.

        Args:
            source: Path to file (e.g., 'prot.cif') or PDB code (e.g., '1abc:A').
            selection: Atom selection string (e.g., 'name CA').
            pdb_cache: Directory to store downloaded structures.
            remove_alternative_locations: Only keep ' ' or 'A' locations.
            fix_non_standard_aa: Map modified residues to standard AAs.
            remove_water: Exclude HOH molecules.
            remove_hetero: Exclude HETATM records.
            verify: Unused in this version (kept for API compatibility).
            no_exit: If True, raise Exception instead of calling sys.exit.
            create_from_aa: Suppress log messages if called during reconstruction.
        """
        if not create_from_aa:
            logger.debug(_name, f"Initializing structure from: {source}")

        self.atoms = Atoms()
        self.source_file = ""

        # 1. Parse source string for ID and optional chain selection
        words = source.split(":")
        identifier = words[0]
        chains = words[1] if len(words) > 1 else ""

        # 2. Resolve structure path (Local file vs RCSB fetch)
        try:
            if os.path.exists(identifier):
                self.source_file = identifier
            else:
                self.source_file = self.fetch(identifier, pdb_cache)
        except Exception as e:
            if no_exit:
                raise Pdb.InvalidPdbInput(str(e))
            logger.exit_program(_name, f"Failed to resolve structure: {identifier}", exc=e)

        # 3. Load coordinates using Biopython
        try:
            is_cif = self.source_file.lower().endswith((".cif", ".cif.gz"))
            parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)

            if self.source_file.endswith(".gz"):
                with gzip.open(self.source_file, "rt") as handle:
                    structure = parser.get_structure("tmp", handle)
            else:
                structure = parser.get_structure("tmp", self.source_file)

            # CABS v3: Always use the first model provided in the file
            self._load_biopython_model(
                structure[0],
                remove_alt=remove_alternative_locations,
                fix_aa=fix_non_standard_aa,
                no_water=remove_water,
                no_hetero=remove_hetero
            )

            # 4. Handle Selections
            if chains:
                self.atoms = self.atoms.select(f"chain {','.join(chains)}")
            if selection:
                self.atoms = self.atoms.select(selection)

            if not len(self.atoms):
                raise Exception(f"Zero atoms were loaded from {source} after selection/filtering.")

        except Exception as e:
            if no_exit:
                raise Pdb.InvalidPdbInput(str(e))
            logger.exit_program(_name, f"Parsing error in {self.source_file}: {e}", exc=e)

    def _load_biopython_model(self, model, remove_alt, fix_aa, no_water, no_hetero):
        """
        Converts Biopython model to CABS Atoms with mmCIF chain mapping.
        Includes safety check for the 62-chain limit.
        """
        # Pool for renaming chains (e.g., 'AA' -> 'A') to maintain Fortran compatibility
        char_pool = iter(string.ascii_uppercase + string.ascii_lowercase + string.digits)
        chain_map = {}

        # Pre-identify existing single-character chains to avoid collisions
        existing_ids = {chain.id for chain in model if len(chain.id) == 1}

        for chain in model:
            chid = chain.id
            # Map multi-char mmCIF chains to 1-char for Fortran compatibility
            if len(chid) > 1 or chid == " ":
                if chid not in chain_map:
                    try:
                        new_id = next(char_pool)
                        while new_id in existing_ids:
                            new_id = next(char_pool)
                        chain_map[chid] = new_id
                        logger.debug(_name, f"Mapping mmCIF chain {chid} -> {new_id}")
                    except StopIteration:
                        logger.exit_program(_name, "Structure exceeds 62 chains (CABS limit).")
                chid = chain_map[chid]

            for residue in chain:
                resname = residue.get_resname()
                if no_water and resname == "HOH": continue
                if no_hetero and residue.id[0] != " ": continue

                if fix_aa and resname not in AA_NAMES.values():
                    resname = AA_SUB_NAMES.get(resname, resname)

                resnum, icode = residue.id[1], residue.id[2].strip()

                for atom in residue:
                    if remove_alt and atom.get_altloc() not in (" ", "A"):
                        continue

                    self.atoms.append(Atom(
                        model=model.id,
                        name=atom.get_name(),  # Use get_name()
                        resname=residue.get_resname(),
                        chid=chid,             # The mapped chain ID
                        resnum=residue.id[1],  # Residue number from Biopython residue tuple
                        icode=residue.id[2].strip(), # Insertion code from Biopython residue tuple
                        coord=Vector3d(atom.get_coord()), # <-- Correct method for coordinates
                        occ=atom.get_occupancy(), # <-- Correct method for Occupancy
                        bfac=atom.get_bfactor(),  # <-- CORRECTED to get_bfactor()
                        hetatm=(residue.id[0] != " ")
                        )
                    )

    def dssp(self, work_dir: str = "", dssp_from_aa: bool = False) -> Dict[str, str]:
        """
        Calculates secondary structure assignment using MDTraj (Pure Python).

        This replaces the external mkdssp binary call. It maps residues using
        the standard CABS format 'resnum[icode]:chid' to ensure compatibility
        with the rest of the simulation pipeline.

        Returns:
            Dict[str, str]: Map of residue IDs (e.g., "123A:A") to SS codes (H, E, C).
        """
        if md is None:
            logger.warning(_name, "MDTraj not found. Defaulting to Coil.")
            return {}

        logger.debug(_name, "Assigning secondary structure via MDTraj...")
        try:
            # MDTraj handles PDB, mmCIF, and .gz handles automatically
            traj = md.load(self.source_file)
            labels = md.compute_dssp(traj, simplified=True)[0]

            sec = {}
            for i, res in enumerate(traj.topology.residues):
                # Ensure the key format matches Atom.resid_id(): "resnum[icode]:chid"
                icode = getattr(res, 'insertion_code', '')
                key = f"{res.resSeq}{icode}:{res.chain.chain_id}"
                sec[key] = labels[i]

            logger.debug(_name, "DSSP assignment was performed with MDTraj.")

            # If log level is high enough, save the SS string to a file for the user
            if work_dir and logger.output_dssp():
                dssp_out = os.path.join(work_dir, "output_data", "DSSP_output.txt")
                logger.to_file(dssp_out, "".join(labels), f"Saved SS sequence to {dssp_out}")

            return sec
        except Exception as e:
            logger.warning(_name, f"MDTraj-DSSP assignment failed: {e}. Reverting to defaults.")
            return {}

    @staticmethod
    def fetch(identifier: str, pdb_cache: str) -> str:
        """
        Fetches structure coordinates using an API-validated mirror strategy.

        All URLs are retrieved from cabs_constants.json. Uses RCSB Data API
        to determine if a PDB format fallback is available.
        """
        # 1. Load Configuration
        try:
            remote_cfg = get_config_section("cabs_constants", "remote_services")
            mirrors = remote_cfg["mirrors"]
            api_base = remote_cfg["rcsb_api"]
        except (KeyError, FileNotFoundError):
            # Fail-safe defaults
            mirrors = ["https://files.rcsb.org/download/"]
            api_base = "https://data.rcsb.org/rest/v1/core/entry/"

        entry_id = identifier.lower()

        # 2. Setup Cache
        cache_dir = os.path.join(os.path.expanduser(pdb_cache), ".cabsPDBcache")
        os.makedirs(cache_dir, exist_ok=True)
        cif_path = os.path.join(cache_dir, f"{entry_id}.cif.gz")
        pdb_path = os.path.join(cache_dir, f"{entry_id}.pdb.gz")

        # 3. Check Local Cache
        if os.path.exists(cif_path): return cif_path
        if os.path.exists(pdb_path): return pdb_path

        # 4. API Validation
        logger.info(_name, f"Querying RCSB API for: {entry_id}")
        can_use_pdb = False
        try:
            api_resp = req.get(f"{api_base}{entry_id}", timeout=10)
            if api_resp.status_code == 200:
                metadata = api_resp.json()
                db_status = metadata.get("pdbx_database_status", {})
                can_use_pdb = (db_status.get("pdb_format_compatible", "N") == "Y")
            else:
                logger.warning(_name, f"Entry {entry_id} not indexed in API. Trying mmCIF only.")
        except req.RequestException:
            # Fallback guessing if API is unreachable (CSMs/AlphaFold don't start with digits)
            can_use_pdb = not entry_id.startswith(("af-", "ma-"))

        # 5. Multi-Mirror Download Strategy
        for base_url in mirrors:
            try:
                # Priority 1: mmCIF (Modern Standard)
                dl_url = f"{base_url}{entry_id}.cif.gz"
                r = req.get(dl_url, timeout=15)
                if r.status_code == 200:
                    with open(cif_path, "wb") as f: f.write(r.content)
                    return cif_path

                # Priority 2: PDB (Fallback if compatible)
                if can_use_pdb:
                    dl_url = f"{base_url}{entry_id}.pdb.gz"
                    r = req.get(dl_url, timeout=15)
                    if r.status_code == 200:
                        with open(pdb_path, "wb") as f: f.write(r.content)
                        return pdb_path
            except req.RequestException:
                continue # Try next mirror

        raise ConnectionError(
            f"Failed to fetch {entry_id} from all configured mirrors. "
            "Check internet connection or ID validity."
        )

    def mk_ss_header(self, dssp_from_aa: bool = False) -> str:
        """
        Generates PDB-compliant HELIX and SHEET records.

        This reconstructs the header information needed by visualization tools
        based on the secondary structure assigned by MDTraj.
        """
        dssp_data = self.dssp(dssp_from_aa=dssp_from_aa)
        if not dssp_data: return ""

        def identify_boundaries(ss_type):
            """Closure to find start and end of contiguous blocks."""
            def is_start(t):
                prev, curr, _ = t
                return curr[1] == ss_type and (prev is None or prev[1] != ss_type)
            def is_end(t):
                _, curr, next_ = t
                return curr[1] == ss_type and (next_ is None or next_[1] != ss_type)
            return is_start, is_end

        def sliding_window(seq):
            """Creates a window for boundary detection."""
            s = list(seq)
            padded = [None] + s + [None]
            return zip(padded, padded[1:], padded[2:])

        items = list(dssp_data.items())
        windows = list(sliding_window(items))
        get_id = lambda t: t[1][0]

        output = []
        ca_atoms = self.atoms.select("NAME CA")

        for ss_type, label in [("H", "HELIX"), ("E", "SHEET")]:
            is_start, is_end = identify_boundaries(ss_type)
            ranges = list(zip(map(get_id, filter(is_start, windows)),
                              map(get_id, filter(is_end, windows))))

            for i, (start, end) in enumerate(ranges, 1):
                try:
                    s_num, s_ch = start.split(":")
                    e_num, e_ch = end.split(":")
                    sa = self.atoms.select(f"RESNUM {s_num} and CHAIN {s_ch} and NAME CA")[0]
                    ea = self.atoms.select(f"RESNUM {e_num} and CHAIN {e_ch} and NAME CA")[0]

                    if label == "HELIX":
                        l = ca_atoms.atoms.index(ea) - ca_atoms.atoms.index(sa) + 1
                        # PDB HELIX Record Format
                        output.append(f"HELIX  {i:3d} {i:3d} {sa.resname:>3} {s_ch:1} {sa.resnum:4d}{sa.icode:1} {ea.resname:>3} {e_ch:1} {ea.resnum:4d}{ea.icode:1} {1:2d}                               {l:5d}\n")
                    else:
                        # PDB SHEET Record Format
                        output.append(f"SHEET  {i:3d} {i:3d} 1 {sa.resname:>3} {s_ch:1}{sa.resnum:4d}{sa.icode:1} {ea.resname:>3} {e_ch:1}{ea.resnum:4d}{ea.icode:1}  0\n")
                except (IndexError, ValueError):
                    continue

        return "".join(output)

    def save_initial_pdb(self, work_dir: str = "") -> None:
        """Saves the loaded coordinates as a standard PDB file."""
        if work_dir:
            header = self.mk_ss_header()
            path = os.path.join(work_dir, "output_pdbs", "start_all.pdb")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            self.atoms.save_to_pdb(path, header=header)

    def __repr__(self) -> str:
        return f"<PdbLoader source={self.source_file}>"
