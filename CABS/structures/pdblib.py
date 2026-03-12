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
from CABS.constants import AA_NAMES, AA_SUB_NAMES, AA_SUB_NAMES_EXTENDED
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
                no_hetero=remove_hetero,
                chains=chains
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

    def _load_biopython_model(self, model, remove_alt, fix_aa, no_water, no_hetero, chains=""):
        """
        Converts Biopython model to CABS Atoms with mmCIF chain mapping.
        Includes safety check for the 62-chain limit.
        """
        # Pool for renaming chains (e.g., 'AA' -> 'A') to maintain Fortran compatibility
        char_pool = iter(string.ascii_uppercase + string.ascii_lowercase + string.digits)
        chain_map = {}

        # Pre-identify existing single-character chains to avoid collisions
        existing_ids = {chain.id for chain in model if len(chain.id) == 1}

        # If chains selection is provided, we only load those
        selected_chains = set(chains) if chains else None

        for chain in model:
            chid = chain.id
            
            # Map multi-char mmCIF chains to 1-char for Fortran compatibility
            if len(chid) > 1 or chid == " " or chid == "":
                if chid not in chain_map:
                    try:
                        new_id = next(char_pool)
                        while new_id in existing_ids:
                            new_id = next(char_pool)
                        chain_map[chid] = new_id
                        logger.debug(_name, f"Mapping mmCIF chain {chid} -> {new_id}")
                    except StopIteration:
                        # We only fail if we actually NEEDED this chain and it's over the limit
                        # But here we don't know yet if it's selected. 
                        # We'll defer the limit check until we know it's being used.
                        pass
                chid = chain_map.get(chid, chid)

            # Check if this chain is in our selection (if any)
            if selected_chains and chid not in selected_chains and chain.id not in selected_chains:
                continue

            # If it's selected and still multi-char (and not mapped yet), we have a problem
            if len(chid) > 1 or chid == " " or chid == "":
                 logger.exit_program(_name, "Structure exceeds 62 chains (CABS limit) or has invalid chain IDs.")

            for residue in chain:
                resname = residue.get_resname()
                
                # Rescuing non-standard AAs (which are often marked as HETATMs)
                if fix_aa and resname not in AA_NAMES.values():
                    # Priority 1: Use extended mapping (maps non-standard to standard types)
                    sname = AA_SUB_NAMES_EXTENDED.get(resname)
                    if sname:
                        resname = AA_NAMES.get(sname, resname)
                    else:
                        # Priority 2: Use standard mapping (maps standard types to 1-letter codes)
                        # This handles cases where resname is already a standard 3-letter code
                        resname = AA_SUB_NAMES.get(resname, resname)

                if no_water and resname == "HOH": continue
                
                # Drop heteroatoms unless they were recognized and translated into standard AAs
                # Biopython residue.id[0] is ' ' for ATOM, 'H_RES' for HETATM, 'W' for WATER
                is_hetero = residue.id[0] != " "
                if no_hetero and is_hetero and resname not in AA_NAMES.values():
                    continue

                for atom in residue:
                    if remove_alt and atom.get_altloc() not in (" ", "A"):
                        continue

                    self.atoms.append(Atom(
                        model=model.id,
                        name=atom.get_name(),
                        resname=resname,
                        chid=chid,             # The mapped chain ID
                        resnum=residue.id[1],
                        icode=residue.id[2].strip(),
                        coord=Vector3d(atom.get_coord()),
                        occ=atom.get_occupancy(),
                        bfac=atom.get_bfactor(),
                        hetatm=(residue.id[0] != " ")
                        )
                    )

        # Final check if we have any atoms
        if not self.atoms:
            raise self.InvalidPdbInput("No atoms loaded based on selection criteria.")

        # Check if the total number of UNIQUE mapped chains exceeds 62
        unique_mapped_chains = {a.chid for a in self.atoms}
        if len(unique_mapped_chains) > 62:
             logger.exit_program(_name, "Structure exceeds 62 chains (CABS limit).")

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
            # Generate a temporary PDB for MDTraj to ensure it only sees the protein atoms
            # that CABS has already filtered (removing waters, DNA, etc.)
            temp_dir = work_dir or "."
            os.makedirs(os.path.join(temp_dir, "output_pdbs"), exist_ok=True)
            load_file = os.path.join(temp_dir, "output_pdbs", "dssp_topology.pdb")
            self.atoms.save_to_pdb(load_file)
            
            logger.info(_name, f"MDTraj loading filtered topology from: {load_file}")
            traj = md.load(load_file)
            labels = md.compute_dssp(traj, simplified=True)[0]

            sec = {}
            assigned_ss = ""
            last_chain = None
            for i, res in enumerate(traj.topology.residues):
                # Track chains to insert '+' separator for multi-chain output
                current_chain = res.chain.chain_id
                if last_chain is not None and current_chain != last_chain:
                    assigned_ss += "+"
                last_chain = current_chain

                # Ensure the key format matches Atom.resid_id(): "resnum[icode]:chid"
                icode = getattr(res, "insertion_code", "") or ""
                key = f"{res.resSeq}{icode}:{current_chain}"
                
                # MDTraj returns labels as numpy characters; convert to native str
                # Also map 'NA' or unknowns to 'C' (Coil) to avoid KeyError in CABS
                label = str(labels[i])
                if label not in ["H", "E", "T", "C"]:
                     label = "C"
                sec[key] = label
                assigned_ss += label

            logger.info(_name, "DSSP assignment was performed with MDTraj.")

            # If log level is high enough, save the SS string to a file for the user
            if work_dir and logger.output_dssp():
                dssp_out = os.path.join(work_dir, "output_data", "DSSP_output.txt")
                os.makedirs(os.path.dirname(dssp_out), exist_ok=True)
                logger.to_file(dssp_out, assigned_ss, f"Saved SS sequence to {dssp_out}")

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
