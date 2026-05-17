import os
import subprocess
from pathlib import Path
from typing import Optional, List
import mdtraj as md
from concurrent.futures import ProcessPoolExecutor, as_completed

from CABS.io import logger
from CABS.utils.utils import convert_cg_to_all

from CABS.utils.utils import CG2ALL_REPRESENTATIONS
from CABS.reconstruction.cg2all import _read_calpha_atoms, _write_calpha_sc_segment, _format_cg_pdb_line
from CABS.reconstruction.cg2all import sync_residues_with_template, sync_residues

def _task_wrapper(kwargs):
    """Simple wrapper for parallel execution of convert_cg_to_all."""
    return convert_cg_to_all(**kwargs)


def _write_single_model_cg2all_topology(
    filename: str,
    output_file,
    cg_model: str,
) -> None:
    """Write a single-frame coarse-grained topology for batch DCD reconstruction."""
    models = _read_calpha_atoms(filename)
    if not models:
        raise ValueError("No CA atoms found in coarse-grained model.")

    ca_atoms = models[0]
    serial = 1
    if cg_model == "CalphaBasedModel":
        for atom in ca_atoms:
            output_file.write(_format_cg_pdb_line(serial, "CA", atom, atom["coord"]))
            serial += 1
        return

    segment = []
    for atom in ca_atoms:
        if segment and atom["chain"] != segment[-1]["chain"]:
            serial = _write_calpha_sc_segment(output_file, segment, serial)
            segment = []
        segment.append(atom)
    if segment:
        _write_calpha_sc_segment(output_file, segment, serial)

def reconstruct_trajectory(
    topology_pdb: str,
    trajectory_file: str,
    output_dcd: str,
    output_pdb: Optional[str] = None,
    cg_model: str = "CalphaSCModel",
    env_prefix: Optional[str] = None,
    batch_size: Optional[int] = None,
    n_proc: Optional[int] = None,
    device: str = "cpu",
    renumber_flag: bool = False,
) -> None:
    """
    Reconstructs an entire trajectory in one optimized DCD batch call.
    Automatically handles PDB files by converting them to DCD first if needed.
    """
    temp_dcd = None
    temp_output_pdb = None
    input_traj = trajectory_file
    
    # If input is a PDB, convert to DCD and create matching topology
    if trajectory_file.endswith(".pdb"):
        logger.debug(module_name="CG2ALL", msg=f"Converting multi-model PDB {trajectory_file} to temporary DCD...")
        traj = md.load(trajectory_file)
        if traj.n_frames > 1:
            temp_dcd = trajectory_file.replace(".pdb", "_temp.dcd")
            traj.save_dcd(temp_dcd)
            input_traj = temp_dcd
            
            # Create matching single-frame topology PDB using CABS-native icode-aware writer
            temp_top = trajectory_file.replace(".pdb", "_temp_top.pdb")
            with open(temp_top, "w") as f_top:
                _write_single_model_cg2all_topology(trajectory_file, f_top, cg_model)
            topology_pdb = temp_top
            logger.debug(module_name="CG2ALL", msg=f"Created matching topology PDB: {temp_top}")
    
    if env_prefix:
        mamba_exe = os.environ.get("MAMBA_EXE") or os.path.expanduser("~/.local/bin/micromamba")
        if not os.path.exists(mamba_exe):
            mamba_exe = "micromamba"
        command_parts = [
            mamba_exe, "run", "-p", env_prefix,
            "convert_cg2all"
        ]
    else:
        command_parts = ["convert_cg2all"]

    command_parts.extend([
        "-p", topology_pdb,
        "-d", input_traj,
        "-o", output_dcd,
    ])

    if output_pdb:
        temp_output_pdb = output_pdb.replace(".pdb", "_last_frame.pdb")
        command_parts.extend(["-opdb", temp_output_pdb])

    # Determine dynamic defaults for batch and proc based on hardware
    if n_proc is None or n_proc <= 0:
        # Cap data loading to 4 threads; more than that is usually overkill for I/O
        n_proc = min(4, os.cpu_count() or 1)
    if batch_size is None or batch_size <= 0:
        # Scale compute batching dynamic with hardware
        batch_size = max(1, (os.cpu_count() or 1) // 2)

    command_parts.extend([
        "--cg", cg_model,
        "--batch", str(batch_size),
        "--proc", str(n_proc),
        "--device", "cpu",
        "--fix"
    ])

    logger.info(module_name="CG2ALL", msg=f"Running batch trajectory reconstruction: {' '.join(command_parts)}")
    
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    
    try:
        subprocess.run(
            command_parts,
            shell=False,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        
        # Convert the authoritative reconstructed DCD into the final multi-model
        # PDB using the last-frame all-atom PDB only as topology.
        if output_pdb and temp_output_pdb and os.path.exists(output_dcd) and os.path.exists(temp_output_pdb):
            logger.debug(module_name="CG2ALL", msg=f"Converting reconstructed DCD to final PDB: {output_pdb}")
            aa_traj = md.load(output_dcd, top=temp_output_pdb)
            aa_traj.save_pdb(output_pdb)
        # Always synchronize residues and chains from the input coarse-grained trajectory.
        # This ensures that CABS-assigned metadata is preserved in the all-atom reconstruction.
        if output_pdb and os.path.exists(output_pdb):
            logger.debug(module_name="CG2ALL", msg=f"Synchronizing chains and residues for: {output_pdb}")
            sync_residues(
                input_pdb_path=Path(trajectory_file),
                output_pdb_path=Path(output_pdb),
            )
            
            # If renumbering to original PDB is requested, perform a second pass with the topology template.
            if renumber_flag:
                # Use topology_pdb as the AA template for precise renumbering
                sync_residues_with_template(
                    input_pdb_path=Path(trajectory_file),
                    topology_pdb_path=Path(topology_pdb),
                    output_pdb_path=Path(output_pdb),
                )
            if os.path.exists(output_dcd):
                os.remove(output_dcd)
            if os.path.exists(temp_output_pdb):
                os.remove(temp_output_pdb)
            
    except subprocess.CalledProcessError as e:
        logger.critical(module_name="CG2ALL", msg=f"Batch reconstruction failed: {e.stderr}")
        raise Exception("cg2all batch reconstruction failed")
    finally:
        if temp_dcd and os.path.exists(temp_dcd):
            os.remove(temp_dcd)
        if 'temp_top' in locals() and temp_top and os.path.exists(temp_top):
            os.remove(temp_top)
        if temp_output_pdb and os.path.exists(temp_output_pdb):
            os.remove(temp_output_pdb)

def reconstruct_parallel(
    input_files: List[str],
    output_dir: str,
    workers: int = 4,
    env_prefix: Optional[str] = None,
    cg_model: str = "CalphaBasedModel",
    renumber: bool = False,
    reference_pdb: Optional[str] = None,
    work_dir: str = "."
) -> None:
    """
    Parallel batch processor for individual PDB files (clusters or replicas).
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    tasks = []
    for i, input_file in enumerate(input_files):
        out_name = Path(input_file).name
        tasks.append({
            "filename": input_file,
            "work_dir": work_dir,
            "iter": i,
            "reference_pdb": reference_pdb,
            "renumber_flag": renumber,
            "env_prefix": env_prefix,
            "output_filename": out_name,
            "cg2all_representation": cg_model
        })

    logger.info(module_name="CG2ALL", msg=f"Launching {len(tasks)} parallel reconstruction tasks with {workers} workers.")
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_task_wrapper, task): task for task in tasks}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                task = futures[future]
                logger.warning(module_name="CG2ALL", msg=f"Task {task['filename']} failed: {e}")

def reconstruct_job_outputs(job, output_folder: str) -> None:
    """
    Orchestrates parallel all-atom reconstruction for Clusters and Trajectories.
    This function scans the output folder for saved CA models and triggers reconstruction.
    """
    if job.aa_method != "cg2all":
        return
    
    import glob
    mode = str(job.aa_rebuild).upper()

    # 1. Parallel Cluster Reconstruction (Mode C or A)
    # Map internal CABS representation names to cg2all ones
    # Resolve representation to cg2all-native name
    cg_model = CG2ALL_REPRESENTATIONS.get(job.cg2all_representation, "CalphaSCModel")

    if "C" in mode or "A" in mode:
        # Scan for cluster files saved by job.save_models()
        cluster_files = sorted(
            path
            for path in glob.glob(os.path.join(output_folder, "cluster_*.pdb"))
            if "_all_atom" not in Path(path).stem
        )
        
        if cluster_files:
            logger.info(module_name="CG2ALL", msg=f"Launching batch reconstruction for {len(cluster_files)} clusters...")
            for cluster_file in cluster_files:
                cluster_path = Path(cluster_file)
                batch_name = cluster_path.stem
                out_dcd = os.path.join(output_folder, f"{batch_name}_all_atom.dcd")
                out_pdb = os.path.join(output_folder, f"{batch_name}_all_atom.pdb")
                
                # Use the cluster file itself as its own topology for the batch run
                reconstruct_trajectory(
                    topology_pdb=cluster_file,
                    trajectory_file=cluster_file,
                    output_dcd=out_dcd,
                    output_pdb=out_pdb,
                    cg_model=cg_model,
                    env_prefix=job.cg2all_env_prefix,
                    batch_size=getattr(job, "batch_size", None),
                    n_proc=getattr(job, "aa_rebuild_workers", None),
                    renumber_flag=job.renumber,
                )

    # 2. Optimized Trajectory Reconstruction (DCD Native Path) (Mode T or A)
    if "T" in mode or "A" in mode:
        # Scan for replica files saved by job.save_models()
        replica_files = sorted(
            path
            for path in glob.glob(os.path.join(output_folder, "replica_*.pdb"))
            if "_all_atom" not in Path(path).stem
        )
        # Also check for single-replica naming
        if not replica_files:
            if os.path.exists(os.path.join(output_folder, "replica.pdb")):
                replica_files = [os.path.join(output_folder, "replica.pdb")]
        
        if replica_files:
            logger.info(module_name="CG2ALL", msg=f"Launching optimized batch reconstruction for {len(replica_files)} trajectories...")
            for replica_file in replica_files:
                replica_path = Path(replica_file)
                batch_name = replica_path.stem
                out_dcd = os.path.join(output_folder, f"{batch_name}_all_atom.dcd")
                out_pdb = os.path.join(output_folder, f"{batch_name}_all_atom.pdb")
                
                # Use Medoid 0 as a stable single-model topology for the DCD
                topology_pdb = os.path.join(output_folder, "model_0.pdb")
                if not os.path.exists(topology_pdb):
                     # Fallback to start.pdb if medoids were not built
                     topology_pdb = os.path.join(output_folder, "start.pdb")
                
                reconstruct_trajectory(
                    topology_pdb=topology_pdb,
                    trajectory_file=str(replica_path),
                    output_dcd=out_dcd,
                    output_pdb=out_pdb,
                    cg_model=cg_model,
                    env_prefix=job.cg2all_env_prefix,
                    batch_size=getattr(job, "batch_size", None),
                    n_proc=getattr(job, "aa_rebuild_workers", None),
                    renumber_flag=job.renumber,
                )
