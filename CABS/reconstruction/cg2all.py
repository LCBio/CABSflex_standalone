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
        if len(input_ca_data) != len(output_ca_indices):
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

    topology_lines = topology_pdb_path.read_text().splitlines()
    residue_atom_counts: List[int] = []
    current_residue = None
    current_count = 0

    for line in topology_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resid = (line[21], line[22:26], line[26], line[17:21])
        if current_residue is None:
            current_residue = resid
            current_count = 1
        elif resid == current_residue:
            current_count += 1
        else:
            residue_atom_counts.append(current_count)
            current_residue = resid
            current_count = 1
    if current_residue is not None:
        residue_atom_counts.append(current_count)

    if not residue_atom_counts:
        raise ValueError("No atom residues found in topology helper PDB for residue synchronization.")

    output_lines = output_pdb_path.read_text().splitlines(keepends=True)
    output_models: List[List[int]] = []
    current_atom_indices: List[int] = []
    saw_output_models = False

    for idx, line in enumerate(output_lines):
        if line.startswith("MODEL"):
            saw_output_models = True
            current_atom_indices = []
            continue
        if line.startswith("ENDMDL"):
            if current_atom_indices:
                output_models.append(current_atom_indices)
            current_atom_indices = []
            continue
        if line.startswith(("ATOM", "HETATM")):
            current_atom_indices.append(idx)

    if current_atom_indices:
        output_models.append(current_atom_indices)

    if not output_models:
        raise ValueError("No atom records found in reconstructed PDB for residue synchronization.")

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
            f"to {output_pdb_path.name} using topology {topology_pdb_path.name}"
        ),
    )

    atoms_per_model = sum(residue_atom_counts)
    for model_idx, (input_ca_data, atom_indices) in enumerate(zip(input_models, output_models), start=1):
        if len(input_ca_data) != len(residue_atom_counts):
            raise ValueError(
                "Residue count mismatch during synchronization: "
                f"reference model {model_idx} has {len(input_ca_data)} residues, "
                f"topology has {len(residue_atom_counts)}."
            )
        if len(atom_indices) != atoms_per_model:
            raise ValueError(
                "Atom count mismatch during synchronization: "
                f"output model {model_idx} has {len(atom_indices)} atoms, "
                f"topology expects {atoms_per_model}."
            )

        cursor = 0
        for residue_info, atom_count in zip(input_ca_data, residue_atom_counts):
            target_resnum = residue_info["resnum"]
            target_chid = residue_info["chain"]
            target_icode = residue_info["icode"]
            for atom_idx in atom_indices[cursor:cursor + atom_count]:
                line = output_lines[atom_idx]
                output_lines[atom_idx] = (
                    f"{line[:21]}{target_chid}{target_resnum:4d}{target_icode}{line[27:]}"
                )
            cursor += atom_count

    output_pdb_path.write_text("".join(output_lines))
    return "Residues synchronized using topology template"


def convert_cg_to_all(
    filename: Union[str, TextIO],
    work_dir: str = ".",
    iter: int = 0,
    reference_pdb: Optional[str] = None,
    renumber_flag: bool = False,
    env_prefix: Optional[str] = None,
    output_filename: Optional[str] = None,
    cg2all_representation: str = "calpha",
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
        start_path = Path(work_dir) / "output_pdbs" / "start.pdb"
        start_all_path = Path(work_dir) / "output_pdbs" / "start_all.pdb"
        if start_path.exists():
            reference_path = start_path
        elif start_all_path.exists():
            reference_path = start_all_path
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
        
        # After successful conversion, synchronize residues in the OUTPUT file
        output_file_path = output_dir / fout
        if output_file_path.exists() and renumber_flag:
            sync_residues(input_pdb_path=reference_path, output_pdb_path=output_file_path)
            
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
