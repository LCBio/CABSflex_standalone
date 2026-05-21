import os
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import (
    Any,
    Dict,
    List,
    Optional,
    TextIO,
    Union,
)
from contextlib import closing
import numpy as np
import numpy.typing as npt

from CABS.io import logger
from CABS.utils.utils import SCModeler, CG2ALL_REPRESENTATIONS


def _format_cg_pdb_line(
    serial: int,
    atom_name: str,
    source: Dict[str, Any],
    coord: npt.NDArray[np.float64],
) -> str:
    fmt_name = f" {atom_name:<3s}"
    if len(atom_name) == 4:
        fmt_name = atom_name
    return (
        f"ATOM  {serial:5d} {fmt_name:4s}{source['alt']:1s}"
        f"{source['resname']:<4s}{source['chain']:1s}{source['resnum']:4d}"
        f"{source['icode']:1s}   {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
        f"{source['occ']:6.2f}{source['bfac']:6.2f} {source['tail']}\n"
    )


def _read_calpha_atoms(filename: Union[str, TextIO]) -> List[List[Dict[str, Any]]]:
    all_models = []
    current_model = []
    f_in = open(filename, "r") if isinstance(filename, str) else closing(filename)
    with f_in as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = []
                continue
            if line.startswith("ENDMDL"):
                if current_model:
                    all_models.append(current_model)
                current_model = []
                continue
            if not line.startswith(("ATOM", "HETATM")) or line[12:16].strip() != "CA":
                continue
            current_model.append(
                {
                    "alt": line[16],
                    "resname": line[17:21].strip(),
                    "chain": line[21],
                    "resnum": int(line[22:26]),
                    "icode": line[26],
                    "coord": np.array(
                        [
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ]
                    ),
                    "occ": float(line[54:60]),
                    "bfac": float(line[60:66]),
                    "tail": line[67:].rstrip("\n") if len(line) > 67 else "",
                }
            )
        # Handle case where there are no MODEL/ENDMDL records
        if current_model and not all_models:
            all_models.append(current_model)
    return all_models


def _write_calpha_sc_segment(
    output_file: TextIO,
    ca_atoms: List[Dict[str, Any]],
    serial: int,
) -> int:
    if len(ca_atoms) < 3:
        sc_coords = np.array([atom["coord"] for atom in ca_atoms])
    else:
        nms = [
            type("ResidueName", (), {"resname": atom["resname"]})()
            for atom in ca_atoms
        ]
        ca_coords = np.array([atom["coord"] for atom in ca_atoms])
        sc_coords = SCModeler(nms).rebuild_one(ca_coords, sc=True)
    for atom, sc_coord in zip(ca_atoms, sc_coords):
        output_file.write(_format_cg_pdb_line(serial, "CA", atom, atom["coord"]))
        serial += 1
        output_file.write(_format_cg_pdb_line(serial, "SC", atom, sc_coord))
        serial += 1
    return serial


def _write_cg2all_input_pdb(
    filename: Union[str, TextIO],
    output_file: TextIO,
    cg2all_representation: str,
) -> None:
    models = _read_calpha_atoms(filename)
    if not models:
        raise ValueError("No CA atoms found in coarse-grained model.")
    if cg2all_representation not in CG2ALL_REPRESENTATIONS:
        raise ValueError(f"Unsupported cg2all representation: {cg2all_representation}")

    # Renumber residues sequentially per-chain and strip insertion codes.
    # Chain IDs are preserved so cg2all can detect chain breaks in multichain structures.
    for ca_atoms in models:
        cur_chain = None
        r_idx = 0
        for atom in ca_atoms:
            if atom["chain"] != cur_chain:
                cur_chain = atom["chain"]
                r_idx = 0
            r_idx += 1
            atom["resnum"] = r_idx
            atom["icode"] = " "

    for i, ca_atoms in enumerate(models):
        if len(models) > 1:
            output_file.write(f"MODEL     {i+1:4d}\n")
        
        serial = 1
        if cg2all_representation == "calpha":
            for atom in ca_atoms:
                output_file.write(_format_cg_pdb_line(serial, "CA", atom, atom["coord"]))
                serial += 1
        else:
            segment = []
            for atom in ca_atoms:
                if segment and atom["chain"] != segment[-1]["chain"]:
                    serial = _write_calpha_sc_segment(output_file, segment, serial)
                    segment = []
                segment.append(atom)
            if segment:
                _write_calpha_sc_segment(output_file, segment, serial)
        
        if len(models) > 1:
            output_file.write("ENDMDL\n")


def sync_residues(input_pdb_path: Path, output_pdb_path: Path) -> str:
    """Synchronize residue numbering and chain IDs between input and output PDB files."""
    input_models: List[List[Dict[str, Any]]] = []
    current_model: List[Dict[str, Any]] = []
    saw_model_records = False

    for line in input_pdb_path.read_text().splitlines():
        if line.startswith("MODEL"):
            saw_model_records = True
            current_model = []
            continue
        if line.startswith("ENDMDL"):
            if current_model:
                input_models.append(current_model)
            current_model = []
            continue
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            current_model.append(
                {
                    "resnum": int(line[22:26]),
                    "chid": line[21],
                    "icode": line[26],
                }
            )

    if current_model:
        input_models.append(current_model)

    if not input_models:
        raise ValueError("No CA atoms found in reference PDB for residue synchronization.")

    if not saw_model_records:
        input_models = [input_models[0]]

    output_lines = output_pdb_path.read_text().splitlines(keepends=True)
    output_models: List[tuple[int, int, List[int]]] = []
    current_ca_indices: List[int] = []
    model_start = 0
    saw_output_models = False

    for idx, line in enumerate(output_lines):
        if line.startswith("MODEL"):
            saw_output_models = True
            model_start = idx
            current_ca_indices = []
            continue
        if line.startswith("ENDMDL"):
            if current_ca_indices:
                output_models.append((model_start, idx, current_ca_indices))
            current_ca_indices = []
            continue
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            current_ca_indices.append(idx)

    if current_ca_indices:
        output_models.append((model_start, len(output_lines), current_ca_indices))

    if not output_models:
        raise ValueError("No CA atoms found in reconstructed PDB for residue synchronization.")

    if not saw_output_models:
        output_models = [(0, len(output_lines), output_models[0][2])]

    if len(input_models) == 1 and len(output_models) > 1:
        input_models = [input_models[0] for _ in output_models]
    elif len(input_models) != len(output_models):
        raise ValueError(
            "Model count mismatch during residue synchronization: "
            f"reference has {len(input_models)} model(s), "
            f"output has {len(output_models)} model(s)."
        )

    logger.debug(
        module_name="sync",
        msg=f"Syncing {len(output_models)} model(s) from {input_pdb_path.name} to {output_pdb_path.name}",
    )

    for model_idx, ((model_start, model_end, output_ca_indices), input_ca_data) in enumerate(zip(output_models, input_models), start=1):
        if len(input_ca_data) < len(output_ca_indices):
            logger.debug(
                module_name="sync",
                msg=(
                    f"Reference model {model_idx} has fewer CA atoms ({len(input_ca_data)}) "
                    f"than output model ({len(output_ca_indices)}). Syncing prefix only."
                )
            )
            output_ca_indices = output_ca_indices[:len(input_ca_data)]
        elif len(input_ca_data) > len(output_ca_indices):
            raise ValueError(
                "Residue renumbering mismatch: "
                f"reference model {model_idx} has {len(input_ca_data)} CA atoms, "
                f"output model {model_idx} has {len(output_ca_indices)}."
            )

        residue_blocks: List[tuple[int, int]] = []
        block_start: Optional[int] = None
        current_resid: Optional[tuple[str, str, str, str]] = None

        for idx in range(model_start, model_end):
            line = output_lines[idx]
            if line.startswith(("ATOM", "HETATM")):
                resid = (line[21], line[22:26], line[26], line[17:21])
                if block_start is None:
                    block_start = idx
                    current_resid = resid
                elif resid != current_resid:
                    residue_blocks.append((block_start, idx))
                    block_start = idx
                    current_resid = resid
            else:
                if block_start is not None:
                    residue_blocks.append((block_start, idx))
                    block_start = None
                    current_resid = None

        if block_start is not None:
            residue_blocks.append((block_start, model_end))

        logger.debug(
            module_name="sync",
            msg=f"Model {model_idx}: {len(residue_blocks)} blocks, {len(input_ca_data)} expected CAs, {len(output_ca_indices)} output CAs"
        )
        if len(residue_blocks) < len(input_ca_data):
            ca_set = set(output_ca_indices)
            split_blocks: List[tuple[int, int]] = []
            for bs, be in residue_blocks:
                cas_in_block = [i for i in range(bs, be) if i in ca_set]
                if len(cas_in_block) <= 1:
                    split_blocks.append((bs, be))
                else:
                    # Split at midpoints between consecutive CAs
                    for k in range(len(cas_in_block)):
                        sub_start = bs if k == 0 else (cas_in_block[k - 1] + cas_in_block[k]) // 2
                        sub_end = be if k == len(cas_in_block) - 1 else (cas_in_block[k] + cas_in_block[k + 1]) // 2
                        split_blocks.append((sub_start, sub_end))
            residue_blocks = split_blocks

        if len(residue_blocks) > len(input_ca_data):
            residue_blocks = residue_blocks[:len(input_ca_data)]

        if len(residue_blocks) != len(input_ca_data):
            raise ValueError(
                "Residue block mismatch during synchronization: "
                f"reference model {model_idx} has {len(input_ca_data)} residues, "
                f"output model {model_idx} has {len(residue_blocks)} residue blocks."
            )

        for (block_start, block_end), target in zip(residue_blocks, input_ca_data):
            target_resnum = target["resnum"]
            target_chid = target["chid"]
            target_icode = target["icode"]

            for idx in range(block_start, block_end):
                line = output_lines[idx]
                if line.startswith(("ATOM", "HETATM", "TER")):
                    output_lines[idx] = (
                        f"{line[:21]}{target_chid}{target_resnum:4d}{target_icode}{line[27:]}"
                    )

    output_pdb_path.write_text("".join(output_lines))
    return "Residues and chains synchronized"


def sync_residues_with_template(
    input_pdb_path: Path,
    topology_pdb_path: Path,
    output_pdb_path: Path,
) -> str:
    """Synchronize chain IDs, residue numbers, and insertion codes using an AA topology template."""
    input_models = _read_calpha_atoms(str(input_pdb_path))
    if not input_models:
        raise ValueError("No CA atoms found in reference PDB for residue synchronization.")

    output_lines = output_pdb_path.read_text().splitlines(keepends=True)
    output_models: List[tuple[int, int, List[int]]] = []
    current_ca_indices: List[int] = []
    model_start = 0
    saw_output_models = False

    for idx, line in enumerate(output_lines):
        if line.startswith("MODEL"):
            saw_output_models = True
            model_start = idx
            current_ca_indices = []
            continue
        if line.startswith("ENDMDL"):
            if current_ca_indices:
                output_models.append((model_start, idx, current_ca_indices))
            current_ca_indices = []
            continue
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            current_ca_indices.append(idx)

    if current_ca_indices:
        output_models.append((model_start, len(output_lines), current_ca_indices))

    if not output_models:
        raise ValueError("No CA atoms found in reconstructed PDB for residue synchronization.")

    if not saw_output_models:
        output_models = [(0, len(output_lines), output_models[0][2])]

    if len(input_models) == 1 and len(output_models) > 1:
        input_models = [input_models[0] for _ in output_models]
    elif len(input_models) != len(output_models):
        raise ValueError(
            "Model count mismatch during residue synchronization: "
            f"reference has {len(input_models)} model(s), output has {len(output_models)} model(s)."
        )

    logger.debug(
        module_name="sync",
        msg=(
            f"Syncing {len(output_models)} model(s) from {input_pdb_path.name} "
            f"to {output_pdb_path.name} using template {topology_pdb_path.name}"
        ),
    )

    for model_idx, ((model_start, model_end, output_ca_indices), input_ca_data) in enumerate(zip(output_models, input_models), start=1):
        if len(input_ca_data) < len(output_ca_indices):
            logger.debug(
                module_name="sync",
                msg=(
                    f"Reference model {model_idx} has fewer CA atoms ({len(input_ca_data)}) "
                    f"than output model ({len(output_ca_indices)}). Syncing prefix only."
                )
            )
            output_ca_indices = output_ca_indices[:len(input_ca_data)]
        elif len(input_ca_data) > len(output_ca_indices):
            raise ValueError(
                "Residue count mismatch during template synchronization: "
                f"reference model {model_idx} has {len(input_ca_data)} CA atoms, "
                f"output model {model_idx} has {len(output_ca_indices)}."
            )

        residue_blocks: List[tuple[int, int]] = []
        block_start = None
        current_resid = None

        for idx in range(model_start, model_end):
            line = output_lines[idx]
            if line.startswith(("ATOM", "HETATM")):
                resid = (line[21], line[22:26], line[26], line[17:21])
                if block_start is None:
                    block_start = idx
                    current_resid = resid
                elif resid != current_resid:
                    residue_blocks.append((block_start, idx))
                    block_start = idx
                    current_resid = resid
            else:
                if block_start is not None:
                    residue_blocks.append((block_start, idx))
                    block_start = None
                    current_resid = None

        if block_start is not None:
            residue_blocks.append((block_start, model_end))

        if len(residue_blocks) < len(input_ca_data):
            ca_set = set(output_ca_indices)
            split_blocks: List[tuple[int, int]] = []
            for bs, be in residue_blocks:
                cas_in_block = [i for i in range(bs, be) if i in ca_set]
                if len(cas_in_block) <= 1:
                    split_blocks.append((bs, be))
                else:
                    for k in range(len(cas_in_block)):
                        sub_start = bs if k == 0 else (cas_in_block[k - 1] + cas_in_block[k]) // 2
                        sub_end = be if k == len(cas_in_block) - 1 else (cas_in_block[k] + cas_in_block[k + 1]) // 2
                        split_blocks.append((sub_start, sub_end))
            residue_blocks = split_blocks

        if len(residue_blocks) > len(input_ca_data):
            residue_blocks = residue_blocks[:len(input_ca_data)]

        if len(residue_blocks) != len(input_ca_data):
            raise ValueError(
                "Residue block mismatch during template synchronization: "
                f"reference model {model_idx} has {len(input_ca_data)} residues, "
                f"output model {model_idx} has {len(residue_blocks)} residue blocks."
            )

        for (block_start, block_end), target in zip(residue_blocks, input_ca_data):
            target_resnum = target["resnum"]
            target_chid = target["chain"]
            target_icode = target["icode"]

            for idx in range(block_start, block_end):
                line = output_lines[idx]
                if line.startswith(("ATOM", "HETATM", "TER")):
                    output_lines[idx] = (
                        f"{line[:21]}{target_chid}{target_resnum:4d}{target_icode}{line[27:]}"
                    )

    output_pdb_path.write_text("".join(output_lines))
    return "Residues synchronized using topology template"


def minimize_pdb_energy(pdb_path: Path) -> None:
    """
    Perform a quick vacuum energy minimization on the given PDB file using OpenMM
    to resolve steric clashes and improve clashscores. Supports both single-model
    and multi-model trajectory structures.
    """
    try:
        from openmm.app import PDBFile, ForceField, Simulation, Modeller, NoCutoff
        from openmm import LangevinIntegrator, Platform
        from openmm.unit import nanometer, picosecond, kelvin
    except ImportError:
        logger.warning(
            module_name="CG2ALL",
            msg="openmm package not found. Skipping energy minimization for reconstructed structure.",
        )
        return

    try:
        # Load PDB file
        pdb = PDBFile(str(pdb_path))
        num_frames = pdb.getNumFrames()
        
        forcefield = ForceField('amber19-all.xml')
        
        # Force the CPU platform for robustness, avoiding GPU memory/driver limits under concurrent workers
        try:
            platform = Platform.getPlatformByName('CPU')
        except Exception:
            platform = None
        
        if num_frames > 1:
            logger.debug(
                module_name="CG2ALL",
                msg=f"Performing vacuum energy minimization frame-by-frame on multi-model PDB: {pdb_path} ({num_frames} models)",
            )
            all_minimized_positions = []
            minimized_topology = None
            
            for i in range(num_frames):
                modeller = Modeller(pdb.topology, pdb.getPositions(frame=i))
                modeller.addHydrogens(forcefield)
                
                system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None)
                integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picosecond)
                
                if platform:
                    simulation = Simulation(modeller.topology, system, integrator, platform)
                else:
                    simulation = Simulation(modeller.topology, system, integrator)
                    
                simulation.context.setPositions(modeller.positions)
                
                simulation.minimizeEnergy(maxIterations=500)
                
                pos = simulation.context.getState(getPositions=True).getPositions()
                all_minimized_positions.append(pos)
                if minimized_topology is None:
                    minimized_topology = modeller.topology
            
            # Save the multi-model structure back to the same path
            with open(str(pdb_path), 'w') as f:
                PDBFile.writeHeader(minimized_topology, f)
                for idx, pos in enumerate(all_minimized_positions):
                    PDBFile.writeModel(minimized_topology, pos, f, modelIndex=idx+1)
                PDBFile.writeFooter(minimized_topology, f)
                
            logger.debug(module_name="CG2ALL", msg=f"Successfully minimized energy and resolved clashes for multi-model PDB: {pdb_path}")
        else:
            logger.debug(module_name="CG2ALL", msg=f"Performing vacuum energy minimization on: {pdb_path}")
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens(forcefield)
            
            system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None)
            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picosecond)
            
            if platform:
                simulation = Simulation(modeller.topology, system, integrator, platform)
            else:
                simulation = Simulation(modeller.topology, system, integrator)
                
            simulation.context.setPositions(modeller.positions)
            
            simulation.minimizeEnergy(maxIterations=500)
            
            with open(str(pdb_path), 'w') as f:
                PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True).getPositions(), f)
                
            logger.debug(module_name="CG2ALL", msg=f"Successfully minimized energy and resolved clashes for: {pdb_path}")
    except Exception as e:
        logger.warning(
            module_name="CG2ALL",
            msg=f"OpenMM energy minimization failed with error: {e}. Reconstructed structure was kept as-is.",
        )


def convert_cg_to_all(
    filename: Union[str, TextIO],
    work_dir: str = ".",
    iter: int = 0,
    reference_pdb: Optional[str] = None,
    renumber_flag: bool = False,
    env_prefix: Optional[str] = None,
    output_filename: Optional[str] = None,
    cg2all_representation: str = "calpha",
    minimize_flag: bool = True,
) -> str:
    """
    Convert coarse-grained model to all-atom
    """
    with NamedTemporaryFile(
        prefix=".", suffix=".pdb", dir=work_dir, mode="w", delete=False
    ) as tmp_file:
        pdb = tmp_file.name
        _write_cg2all_input_pdb(filename, tmp_file, cg2all_representation)

    if renumber_flag:
        reference_path = None
        start_all_path = Path(work_dir) / "output_pdbs" / "start_all.pdb"
        start_path = Path(work_dir) / "output_pdbs" / "start.pdb"
        if start_all_path.exists():
            reference_path = start_all_path
        elif start_path.exists():
            reference_path = start_path
        elif reference_pdb:
            candidate = Path(str(reference_pdb).split(":")[0])
            if candidate.exists():
                reference_path = candidate

        if reference_path is None:
            raise FileNotFoundError(
                "Could not resolve a PDB file for residue renumbering."
            )

        # Note: synchronizing the input file is often ignored by the external tool,
        # so we primarily rely on output synchronization after the run.
        pass

    output_dir = Path(work_dir) / "output_pdbs"
    input_pdb = Path(pdb)
    fout = output_filename or f"model_{iter}.pdb"
    cg_model = CG2ALL_REPRESENTATIONS[cg2all_representation]
    # Modify the subprocess call to use micromamba run if an environment prefix is provided.
    # This ensures all environment variables and dependencies (like torch, dgl) are correctly set up.
    if env_prefix:
        # Try to find micromamba executable
        mamba_exe = os.environ.get("MAMBA_EXE") or os.path.expanduser("~/.local/bin/micromamba")
        if not os.path.exists(mamba_exe):
            mamba_exe = "micromamba"  # Fallback to PATH

        command_parts = [
            mamba_exe, "run", "-p", env_prefix,
            "convert_cg2all",
            "-p", str(input_pdb),
            "-o", str(output_dir / fout),
            "--cg", cg_model,
            "--device", "cpu"
        ]
    else:
        # Fallback to the default path if no specific env is passed
        command_parts = [
            "convert_cg2all",
            "-p", str(input_pdb),
            "-o", str(output_dir / fout),
            "--cg", cg_model,
            "--device", "cpu"
        ]

    # Prepare a clean environment for micromamba run
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    env.pop("PYTHONUSERBASE", None)

    try:
        result = subprocess.run(
            command_parts,
            shell=False,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        
        # After successful conversion, synchronize residues in the OUTPUT file.
        # We always sync from the input CG file first to ensure chains are preserved.
        output_file_path = output_dir / fout
        if output_file_path.exists():
            # Pass 1: Sync to CABS assignments (using original filename as reference)
            ref_path = Path(pdb)
            try:
                sync_residues(input_pdb_path=ref_path, output_pdb_path=output_file_path)
            except Exception as e:
                logger.warning(
                    module_name="CG2ALL",
                    msg=f"Initial CABS assignment residue synchronization warning: {e}. Keeping default assignments."
                )
            
            # Pass 2: Sync to original PDB numbering (using reference_path and start_all.pdb as template)
            if renumber_flag and reference_path:
                start_all_path = Path(work_dir) / "output_pdbs" / "start_all.pdb"
                if start_all_path.exists():
                    try:
                        sync_residues_with_template(
                            input_pdb_path=reference_path,
                            topology_pdb_path=start_all_path,
                            output_pdb_path=output_file_path
                        )
                    except Exception as e:
                        logger.warning(
                            module_name="CG2ALL",
                            msg=f"Template-based residue synchronization warning: {e}. Falling back to standard synchronization."
                        )
                        try:
                            sync_residues(input_pdb_path=reference_path, output_pdb_path=output_file_path)
                        except Exception as e2:
                            logger.warning(
                                module_name="CG2ALL",
                                msg=f"Fallback synchronization also failed: {e2}. Structure was kept with default numbering."
                            )
                else:
                    # Fallback to sync_residues if template is missing
                    try:
                        sync_residues(input_pdb_path=reference_path, output_pdb_path=output_file_path)
                    except Exception as e:
                        logger.warning(
                            module_name="CG2ALL",
                            msg=f"Standard synchronization failed: {e}. Structure was kept with default numbering."
                        )
            
            if minimize_flag:
                minimize_pdb_energy(output_file_path)
            
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
