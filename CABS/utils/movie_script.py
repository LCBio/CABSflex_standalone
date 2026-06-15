"""Generate a standalone make_movies.py script for CABSdock runs."""
from __future__ import annotations

import os
from pathlib import Path
from typing import List


_SCRIPT_HEADER = """\
#!/usr/bin/env python3
import shutil
import subprocess
import sys
from pathlib import Path

WORK_DIR = Path(__file__).parent
N_MODELS = {n_models}
FRAMES_PER_REPLICA = {frames_per_replica}
MODEL0_REPLICA_ID = {model0_replica_id}
PEPTIDE_CHAINS = {peptide_chains!r}
RECEPTOR_CHAINS = {receptor_chains!r}
REFERENCE_PDB = {reference_pdb!r}
HAS_REFERENCE = {has_reference}
"""

_SCRIPT_BODY = '''\

def _require(tool):
    path = shutil.which(tool)
    if not path:
        print(f"{tool} not found in PATH. Install it and re-run.")
        sys.exit(1)
    return path


def _pdb_chains(pdb_path):
    chains = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 22:
                c = line[21].strip()
                if c and c not in chains:
                    chains.append(c)
    return chains


def _count_models(pdb_path):
    return sum(1 for line in open(pdb_path) if line.startswith("MODEL"))


def _chain_sizes(pdb_path):
    sizes = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 22:
                c = line[21].strip()
                if c:
                    sizes[c] = sizes.get(c, 0) + 1
    return sizes


def _resolve_chains(pdb_path, min_atoms=10):
    """Return (rec_chains, pep_chains) from the PDB using RECEPTOR_CHAINS as anchor."""
    sizes = _chain_sizes(pdb_path)
    rec = [c for c in sizes if c in set(RECEPTOR_CHAINS) and sizes[c] >= min_atoms]
    pep = [c for c in sizes if c not in set(RECEPTOR_CHAINS) and sizes[c] >= min_atoms]
    if not rec and sizes:
        ordered = list(sizes.keys())
        pep = [ordered[-1]]
        rec = ordered[:-1]
    return rec, pep


_FACING_HELPER_SCRIPT = """\
# check_facing.py — called by ChimeraX via runscript after view orient.
# Args: model_id  pep_chains(comma-sep)  rec_chains(comma-sep)
# Rotates 180 degrees around Y if the receptor is in front of the peptide
# (i.e. peptide is on the far side from the viewer in the current camera view).
import sys
import numpy as np
from chimerax.core.commands import run

model_id = int(sys.argv[1])
pep_chains = set(sys.argv[2].split(","))
rec_chains = set(sys.argv[3].split(","))

target = next(
    (m for m in session.models if hasattr(m, "atoms") and m.id == (model_id,)),
    None,
)
if target is None:
    raise SystemExit(0)

pep_coords = np.array([a.coord for a in target.atoms if a.residue.chain_id in pep_chains])
rec_coords = np.array([a.coord for a in target.atoms if a.residue.chain_id in rec_chains])
if len(pep_coords) == 0 or len(rec_coords) == 0:
    raise SystemExit(0)

pep_c = pep_coords.mean(axis=0)
rec_c = rec_coords.mean(axis=0)

# camera.position maps camera coords to scene coords.
# Column 2 of the rotation is the camera +Z axis in scene space (pointing toward viewer).
# dot(rec_c - pep_c, cam_z) > 0  →  receptor is closer to viewer than peptide
# → peptide is on the far side → rotate 180 degrees around Y to bring it forward.
cam_z = session.main_view.camera.position.matrix[:3, 2]
if float(np.dot(rec_c - pep_c, cam_z)) > 0:
    run(session, "turn y 180")
"""


def _write_tight_camera(fh, focus_sel, facing_cmd=None):
    fh.write(f"select {focus_sel}\\n")
    fh.write("view orient\\n")
    fh.write("~select\\n")
    if facing_cmd:
        fh.write(facing_cmd + "\\n")
    fh.write("zoom 0.6\\n")


def _write_standard_views(fh, out_dir, prefix):
    for i in range(1, 11):
        fh.write(f\'save "{out_dir / f\'{prefix}snapshot_{i}.png\'}" width 800 height 600\\n\')
        if i < 10:
            fh.write("turn y 36\\n")


def _write_receptor_surface_explicit(fh, model_id, chains):
    sel = f"#{model_id}/" + ",".join(chains)
    fh.write(f"surface {sel}\\n")
    fh.write(f"color {sel} #5b84b1 target a\\n")
    fh.write(f"color {sel} #b3b3b3 target s\\n")
    fh.write(f"transparency {sel} 60 target s\\n")


def _read_csv(path):
    rmsds, energies = [], []
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    rmsds.append(float(parts[0]))
                    energies.append(float(parts[1]))
                except ValueError:
                    continue
    return rmsds, energies


def _find_best_replica():
    csv_files = sorted(WORK_DIR.glob("plots/E_RMSD_*_total.csv"))
    rmsds, _ = _read_csv(csv_files[0])
    best_idx = rmsds.index(min(rmsds))
    return best_idx // FRAMES_PER_REPLICA


def _write_rotation_cxc(model_pdbs, script_path, frames_dir, snapshots_dir, helper_path):
    frames_dir.mkdir(parents=True, exist_ok=True)
    snapshots_dir.mkdir(parents=True, exist_ok=True)
    rec, pep = _resolve_chains(model_pdbs[0])
    check_model_id = 2 if (HAS_REFERENCE and REFERENCE_PDB) else 1
    facing_cmd = f\'runscript "{helper_path}" {check_model_id} {",".join(pep)} {",".join(rec)}\'

    with open(script_path, "w") as fh:
        fh.write("set bgColor white\\n")
        if HAS_REFERENCE and REFERENCE_PDB:
            fh.write(f\'open "{REFERENCE_PDB}" id 1\\n\')
            for i, p in enumerate(model_pdbs):
                fh.write(f\'open "{p}" id {i + 2}\\n\')
            fh.write("dssp\\n")
            rec_spec_ref = "#1/" + ",".join(rec)
            for i in range(len(model_pdbs)):
                match_spec = f"#{i + 2}/" + ",".join(rec)
                fh.write(f"matchmaker {match_spec} to {rec_spec_ref} pairing ss\\n")
            fh.write("hide atoms\\nhide cartoons\\n")
            fh.write("transparency 0 target c\\n")
            fh.write(f"show {rec_spec_ref} cartoon\\n")
            fh.write(f"color {rec_spec_ref} blue\\n")
            _write_receptor_surface_explicit(fh, 1, rec)
            pep_ref_spec = "#1/" + ",".join(pep)
            fh.write(f"show {pep_ref_spec} cartoon\\n")
            fh.write(f"color {pep_ref_spec} green\\n")
            fh.write(f"transparency {pep_ref_spec} 50 target c\\n")
            for i in range(len(model_pdbs)):
                pep_spec = f"#{i + 2}/" + ",".join(pep)
                fh.write(f"show {pep_spec} cartoon\\n")
                fh.write(f"color {pep_spec} red\\n")
            orient_sel = pep_ref_spec
        else:
            for i, p in enumerate(model_pdbs):
                fh.write(f\'open "{p}" id {i + 1}\\n\')
            fh.write("dssp\\n")
            rec_spec_0 = "#1/" + ",".join(rec)
            for i in range(1, len(model_pdbs)):
                match_spec = f"#{i + 1}/" + ",".join(rec)
                fh.write(f"matchmaker {match_spec} to {rec_spec_0} pairing ss\\n")
            fh.write("hide atoms\\nhide cartoons\\n")
            fh.write("transparency 0 target c\\n")
            fh.write(f"show {rec_spec_0} cartoon\\n")
            fh.write(f"color {rec_spec_0} blue\\n")
            _write_receptor_surface_explicit(fh, 1, rec)
            for i in range(len(model_pdbs)):
                pep_spec = f"#{i + 1}/" + ",".join(pep)
                fh.write(f"show {pep_spec} cartoon\\n")
                fh.write(f"color {pep_spec} red\\n")
            orient_sel = "#1/" + ",".join(pep)

        _write_tight_camera(fh, orient_sel, facing_cmd)
        _write_standard_views(fh, snapshots_dir, "top10_")
        _write_tight_camera(fh, orient_sel, facing_cmd)
        for fr in range(1, 721):
            fh.write(f\'turn y 0.5\\nsave "{frames_dir / f\'frame_{fr:04d}.png\'}" width 800 height 600\\n\')


def _write_trajectory_cxc(replica_path, script_path, frames_dir, helper_path):
    frames_dir.mkdir(parents=True, exist_ok=True)
    rec, pep = _resolve_chains(replica_path)
    n_frames = _count_models(replica_path)
    check_model_id = 2 if (HAS_REFERENCE and REFERENCE_PDB) else 1
    facing_cmd = f\'runscript "{helper_path}" {check_model_id} {",".join(pep)} {",".join(rec)}\'

    with open(script_path, "w") as fh:
        fh.write("set bgColor white\\n")
        if HAS_REFERENCE and REFERENCE_PDB:
            fh.write(f\'open "{REFERENCE_PDB}" id 1\\n\')
            fh.write(f\'open "{replica_path}" id 2 coordsets true\\n\')
            fh.write("dssp #1\\n")
            rec_spec_ref = "#1/" + ",".join(rec)
            traj_rec_spec = "#2/" + ",".join(rec)
            fh.write(f"matchmaker {traj_rec_spec} to {rec_spec_ref} pairing ss\\n")
            fh.write("hide atoms\\nhide cartoons\\n")
            fh.write("transparency 0 target c\\n")
            fh.write(f"show {rec_spec_ref} cartoon\\n")
            fh.write(f"color {rec_spec_ref} blue\\n")
            _write_receptor_surface_explicit(fh, 1, rec)
            pep_ref_spec = "#1/" + ",".join(pep)
            fh.write(f"show {pep_ref_spec} cartoon\\n")
            fh.write(f"color {pep_ref_spec} green\\n")
            fh.write(f"transparency {pep_ref_spec} 50 target c\\n")
            traj_pep_spec = "#2/" + ",".join(pep)
            fh.write(f"show {traj_pep_spec} cartoon\\n")
            fh.write(f"color {traj_pep_spec} red\\n")
            orient_sel = pep_ref_spec
        else:
            fh.write(f\'open "{replica_path}" id 1 coordsets true\\n\')
            fh.write("hide atoms\\nhide cartoons\\n")
            fh.write("transparency 0 target c\\n")
            rec_spec = "#1/" + ",".join(rec)
            pep_spec = "#1/" + ",".join(pep)
            fh.write(f"show {rec_spec} cartoon\\n")
            fh.write(f"color {rec_spec} blue\\n")
            _write_receptor_surface_explicit(fh, 1, rec)
            fh.write(f"show {pep_spec} cartoon\\n")
            fh.write(f"color {pep_spec} red\\n")
            orient_sel = pep_spec

        _write_tight_camera(fh, orient_sel, facing_cmd)
        traj_pep = "#2/" + ",".join(pep) if (HAS_REFERENCE and REFERENCE_PDB) else "#1/" + ",".join(pep)
        for fr in range(1, n_frames + 1):
            fh.write(f\'coordset {traj_pep.split("/")[0]} {fr}\\nsave "{frames_dir / f\'frame_{fr:04d}.png\'}" width 800 height 600\\n\')


def _stitch_frames(ffmpeg, frames_dir, output_path, framerate=30):
    subprocess.run(
        [ffmpeg, "-y", "-framerate", str(framerate), "-i", str(frames_dir / "frame_%04d.png"),
         "-c:v", "libx264", "-pix_fmt", "yuv420p", str(output_path)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    )
    for png in frames_dir.glob("*.png"):
        png.unlink()
    try:
        frames_dir.rmdir()
    except OSError:
        pass


def _combine_movies(ffmpeg, left, right, output):
    subprocess.run(
        [ffmpeg, "-y", "-i", str(left), "-i", str(right),
         "-filter_complex", "[0:v]scale=-1:1080[v0];[1:v]scale=-1:1080[v1];[v0][v1]hstack=shortest=1",
         "-c:v", "libx264", "-pix_fmt", "yuv420p", str(output)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    )


def _make_energy_rmsd_animation(csv_path, output_path, ylabel):
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    import numpy as np

    rmsds, energies = _read_csv(csv_path)
    if not rmsds:
        return

    n = len(rmsds)
    orig_idx = np.linspace(0, n - 1, n)
    interp_idx = np.linspace(0, n - 1, n * 2)
    x_interp = np.interp(interp_idx, orig_idx, rmsds)
    y_interp = np.interp(interp_idx, orig_idx, energies)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(rmsds, energies, s=12, alpha=0.4, color="#E83A5D", edgecolors="none")
    dot, = ax.plot([], [], "o", color="red", markersize=14, zorder=10)
    ax.set_xlabel("RMSD (\\u00c5)", fontsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.grid(True, linestyle=":", alpha=0.5)
    plt.tight_layout()

    def update(frame):
        dot.set_data([x_interp[frame]], [y_interp[frame]])
        return [dot]

    anim = animation.FuncAnimation(fig, update, frames=len(x_interp), interval=1000 / 60, blit=True)
    anim.save(str(output_path), writer="ffmpeg", fps=60, dpi=100)
    plt.close(fig)


def main():
    chimerax = _require("chimerax")
    ffmpeg = _require("ffmpeg")
    movies_dir = WORK_DIR / "movies"
    movies_dir.mkdir(exist_ok=True)

    helper_path = movies_dir / "check_facing.py"
    helper_path.write_text(_FACING_HELPER_SCRIPT)

    model_pdbs = sorted(WORK_DIR.glob("output_pdbs/model_*.pdb"))[:N_MODELS]
    if model_pdbs:
        print("Rendering top-10 model rotation...")
        script = movies_dir / "top10_rotation.cxc"
        frames = movies_dir / "top10_frames"
        snapshots = movies_dir / "snapshots"
        _write_rotation_cxc(model_pdbs, script, frames, snapshots, helper_path)
        subprocess.run([chimerax, "--offscreen", "--exit", str(script)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        _stitch_frames(ffmpeg, frames, movies_dir / "Top_10_Models.mp4", framerate=60)

    replica_id = _find_best_replica() if HAS_REFERENCE else MODEL0_REPLICA_ID
    replica_path = None
    for _suffix in ("_all_atom.pdb", "_aa.pdb", ".pdb"):
        _candidate = WORK_DIR / f"output_pdbs/replica_{replica_id}{_suffix}"
        if _candidate.exists():
            replica_path = _candidate
            break
    traj_movie = None
    if replica_path is not None:
        print(f"Rendering trajectory animation (replica {replica_id})...")
        script = movies_dir / f"trajectory_{replica_id}.cxc"
        frames = movies_dir / f"traj_frames_{replica_id}"
        _write_trajectory_cxc(replica_path, script, frames, helper_path)
        subprocess.run([chimerax, "--offscreen", "--exit", str(script)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        traj_movie = movies_dir / f"Trajectory_{replica_id}.mp4"
        _stitch_frames(ffmpeg, frames, traj_movie, framerate=30)

    if HAS_REFERENCE:
        for energy_type, plot_name, combined_name, ylabel in (
            ("total",       "Total_Energy_vs_RMSD",       f"Combined_Trajectory_{replica_id}_Total_Energy_vs_RMSD",       "Total energy"),
            ("interaction", "Interaction_Energy_vs_RMSD", f"Combined_Trajectory_{replica_id}_Interaction_Energy_vs_RMSD", "Interaction energy"),
        ):
            csv_files = sorted(WORK_DIR.glob(f"plots/E_RMSD_*_{energy_type}.csv"))
            if not csv_files:
                continue
            print(f"Generating E vs RMSD animation ({energy_type})...")
            plot_movie = movies_dir / f"{plot_name}.mp4"
            _make_energy_rmsd_animation(csv_files[0], plot_movie, ylabel)
            if traj_movie and traj_movie.exists() and plot_movie.exists():
                print(f"Combining trajectory + {energy_type} energy plot...")
                _combine_movies(ffmpeg, traj_movie, plot_movie, movies_dir / f"{combined_name}.mp4")

    print(f"Done. Movies written to {movies_dir}")


if __name__ == "__main__":
    main()
'''


def write_make_movies_script(
    work_dir: str | Path,
    n_models: int,
    frames_per_replica: int,
    model0_replica_id: int,
    peptide_chains: List[str],
    receptor_chains: List[str],
    has_reference: bool,
    reference_pdb: str = "",
) -> Path:
    work_dir = Path(work_dir)
    header = _SCRIPT_HEADER.format(
        n_models=n_models,
        frames_per_replica=frames_per_replica,
        model0_replica_id=model0_replica_id,
        peptide_chains=peptide_chains,
        receptor_chains=receptor_chains,
        reference_pdb=reference_pdb,
        has_reference=has_reference,
    )
    script_path = work_dir / "make_movies.py"
    script_path.write_text(header + _SCRIPT_BODY, encoding="utf-8")
    return script_path
