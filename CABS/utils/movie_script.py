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
# The reference PDB may use different chain letters than the docked output
# (e.g. original PDB numbering vs. CABS-assigned chains), so its receptor and
# peptide chains are tracked separately and used for all "#1/..." selections.
REFERENCE_RECEPTOR_CHAINS = {reference_receptor_chains!r}
REFERENCE_PEPTIDE_CHAINS = {reference_peptide_chains!r}
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
            rec_spec_ref = "#1/" + ",".join(REFERENCE_RECEPTOR_CHAINS)
            for i in range(len(model_pdbs)):
                match_spec = f"#{i + 2}/" + ",".join(rec)
                fh.write(f"matchmaker {match_spec} to {rec_spec_ref} pairing ss\\n")
            fh.write("hide atoms\\nhide cartoons\\n")
            fh.write("transparency 0 target c\\n")
            fh.write(f"show {rec_spec_ref} cartoon\\n")
            fh.write(f"color {rec_spec_ref} blue\\n")
            _write_receptor_surface_explicit(fh, 1, REFERENCE_RECEPTOR_CHAINS)
            pep_ref_spec = "#1/" + ",".join(REFERENCE_PEPTIDE_CHAINS)
            fh.write(f"show {pep_ref_spec} cartoon\\n")
            fh.write(f"color {pep_ref_spec} green\\n")
            for i in range(len(model_pdbs)):
                pep_spec = f"#{i + 2}/" + ",".join(pep)
                fh.write(f"show {pep_spec} cartoon\\n")
                fh.write(f"color {pep_spec} red\\n")
            # Orient on the whole reference complex (receptor + peptide), not the
            # peptide alone: the peptide's own principal axes are short/noisy and
            # don't reflect the receptor's channel/cleft geometry the peptide sits in.
            orient_sel = "#1"
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
            # Orient on the whole first model (receptor + peptide) for the same reason.
            orient_sel = "#1"

        _write_tight_camera(fh, orient_sel, facing_cmd)
        _write_standard_views(fh, snapshots_dir, "top10_")
        _write_tight_camera(fh, orient_sel, facing_cmd)
        for fr in range(1, 721):
            fh.write(f\'turn y 0.5\\nsave "{frames_dir / f\'frame_{fr:04d}.png\'}" width 800 height 600\\n\')


def _write_trajectory_cxc(replica_path, script_path, frames_dir, helper_path, flexible_receptor=False):
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
            rec_spec_ref = "#1/" + ",".join(REFERENCE_RECEPTOR_CHAINS)
            traj_rec_spec = "#2/" + ",".join(rec)
            fh.write(f"matchmaker {traj_rec_spec} to {rec_spec_ref} pairing ss\\n")
            fh.write("hide atoms\\nhide cartoons\\n")
            fh.write("transparency 0 target c\\n")
            pep_ref_spec = "#1/" + ",".join(REFERENCE_PEPTIDE_CHAINS)
            traj_pep_spec = "#2/" + ",".join(pep)
            if flexible_receptor:
                # Receptor cartoon/surface from the trajectory itself (model #2), so it
                # visibly flexes frame-to-frame instead of staying pinned to the reference.
                fh.write(f"show {traj_rec_spec} cartoon\\n")
                fh.write(f"color {traj_rec_spec} blue\\n")
                _write_receptor_surface_explicit(fh, 2, rec)
                # Orient on the whole trajectory complex (receptor + peptide), not the
                # peptide alone, so multi-chain receptor clefts/channels are captured.
                orient_sel = "#2"
            else:
                fh.write(f"show {rec_spec_ref} cartoon\\n")
                fh.write(f"color {rec_spec_ref} blue\\n")
                _write_receptor_surface_explicit(fh, 1, REFERENCE_RECEPTOR_CHAINS)
                orient_sel = "#1"
            fh.write(f"show {pep_ref_spec} cartoon\\n")
            fh.write(f"color {pep_ref_spec} green\\n")
            fh.write(f"show {traj_pep_spec} cartoon\\n")
            fh.write(f"color {traj_pep_spec} red\\n")
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
            # Orient on the whole complex (receptor + peptide) for the same reason.
            orient_sel = "#1"

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


def _combine_movies(ffmpeg, left, right, output, fps=60):
    # left (trajectory, 30fps) and right (energy plot, 60fps) have matching total
    # duration but different frame rates. Normalize to 60 (not 30): the fps filter
    # only duplicates/drops frames to hit the target rate without changing duration,
    # so 60 keeps the already-smooth energy animation untouched and just duplicates
    # the trajectory's frames, whereas 30 would needlessly drop half the energy frames.
    subprocess.run(
        [ffmpeg, "-y", "-i", str(left), "-i", str(right),
         "-filter_complex",
         f"[0:v]scale=-1:1080,fps={fps}[v0];[1:v]scale=-1:1080,fps={fps}[v1];[v0][v1]hstack=shortest=1",
         "-c:v", "libx264", "-pix_fmt", "yuv420p", str(output)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    )


def _make_energy_rmsd_animation(
    csv_path, output_path, ylabel, replica_id=None, frames_per_replica=None,
    show_other_replicas=True, other_replica_color=None, rmsd_log_scale=False,
):
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import FixedLocator, ScalarFormatter

    rmsds, energies = _read_csv(csv_path)
    if not rmsds:
        return
    rmsds = np.asarray(rmsds)
    energies = np.asarray(energies)

    # The CSV holds all replicas concatenated, replica-major (all frames of
    # replica 0, then replica 1, ...). The moving dot must only trace the
    # selected replica's slice (matching the single-replica trajectory video's
    # duration); the rest of the multi-replica dataset is revealed progressively
    # as background scatter, synchronized by simulation step across all replicas
    # (they were sampled in parallel, so step i of every replica lines up in time).
    n_replicas = None
    if replica_id is not None and frames_per_replica and len(rmsds) % frames_per_replica == 0:
        n_replicas = len(rmsds) // frames_per_replica
        start = replica_id * frames_per_replica
        end = start + frames_per_replica
        traj_rmsds = rmsds[start:end]
        traj_energies = energies[start:end]
        rmsds_by_step = rmsds.reshape(n_replicas, frames_per_replica)
        energies_by_step = energies.reshape(n_replicas, frames_per_replica)
    else:
        traj_rmsds, traj_energies = rmsds, energies
    if len(traj_rmsds) == 0:
        return

    n = len(traj_rmsds)
    orig_idx = np.linspace(0, n - 1, n)
    interp_idx = np.linspace(0, n - 1, n * 2)
    x_interp = np.interp(interp_idx, orig_idx, traj_rmsds)
    y_interp = np.interp(interp_idx, orig_idx, traj_energies)

    fig, ax = plt.subplots(figsize=(8, 6))
    if rmsd_log_scale:
        ax.set_xscale("log")
        # Plain, explicit tick positions -- no default log-scale decade formatting.
        ax.xaxis.set_major_locator(FixedLocator([1, 2.5, 5, 10, 20, 50]))
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_locator(FixedLocator([]))
        x_min = float(rmsds[rmsds > 0].min()) if np.any(rmsds > 0) else 0.01
        ax.set_xlim(max(x_min * 0.8, 0.01), float(rmsds.max()) * 1.2)
        xlabel = "RMSD (\\u00c5, log scale)"
    else:
        # Default matplotlib linear-axis ticks/formatting.
        x_pad = 0.05 * (float(rmsds.max()) - float(rmsds.min()) or 1.0)
        ax.set_xlim(float(rmsds.min()) - x_pad, float(rmsds.max()) + x_pad)
        xlabel = "RMSD (\\u00c5)"
    y_pad = 0.05 * (float(energies.max()) - float(energies.min()) or 1.0)
    ax.set_ylim(float(energies.min()) - y_pad, float(energies.max()) + y_pad)

    own_color = "#E83A5D"
    other_color = other_replica_color or own_color
    # Explicit zorder: other-replica points stay behind the selected replica's
    # own points (and the moving dot), regardless of draw order, so the
    # highlighted replica is never visually buried under the background.
    other_points = ax.scatter([], [], s=12, alpha=0.4, color=other_color, edgecolors="none", zorder=1)
    own_points = ax.scatter([], [], s=12, alpha=0.4, color=own_color, edgecolors="none", zorder=2)
    dot, = ax.plot([], [], "o", color="red", markersize=14, zorder=10)
    ax.set_xlabel(xlabel, fontsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.grid(True, linestyle=":", alpha=0.5)
    plt.tight_layout()

    def update(frame):
        dot.set_data([x_interp[frame]], [y_interp[frame]])
        if n_replicas is not None:
            cur_step = min(int(interp_idx[frame]) + 1, frames_per_replica)
            own_points.set_offsets(np.column_stack([
                rmsds_by_step[replica_id, :cur_step], energies_by_step[replica_id, :cur_step],
            ]))
            if show_other_replicas:
                other_mask = np.ones(n_replicas, dtype=bool)
                other_mask[replica_id] = False
                other_points.set_offsets(np.column_stack([
                    rmsds_by_step[other_mask, :cur_step].ravel(),
                    energies_by_step[other_mask, :cur_step].ravel(),
                ]))
            else:
                other_points.set_offsets(np.empty((0, 2)))
        else:
            own_points.set_offsets(np.column_stack([rmsds, energies]))
            other_points.set_offsets(np.empty((0, 2)))
        return [dot, own_points, other_points]

    anim = animation.FuncAnimation(fig, update, frames=len(x_interp), interval=1000 / 60, blit=True)
    anim.save(str(output_path), writer="ffmpeg", fps=60, dpi=100)
    plt.close(fig)


_THINGS = [
    "top10", "traj-flexible", "traj-rigid", "rmsd-vs-tot-eng", "rmsd-vs-int-eng",
    "combined-tot-eng", "combined-int-eng",
]


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--hide-other-replicas", action="store_true",
                         help="Energy vs RMSD animations: don't show points from replicas other than the selected one")
    parser.add_argument("--other-replica-color", default=None,
                         help="Color for other-replica points (default: same color as the selected replica)")
    parser.add_argument("--rmsd-log-scale", action="store_true",
                         help="Use a log scale (ticks at 1, 2.5, 5, 10, 20, 50) for the RMSD axis in energy vs "
                              "RMSD animations. Default: linear scale with standard matplotlib ticks.")
    parser.add_argument("--generate-only", nargs="+", choices=_THINGS, default=None,
                         help="Generate only these specific outputs instead of everything. "
                              f"Choices: {', '.join(_THINGS)}.")
    args = parser.parse_args()
    only = args.generate_only

    def want(name):
        return only is None or name in only

    chimerax = _require("chimerax")
    ffmpeg = _require("ffmpeg")
    movies_dir = WORK_DIR / "movies"
    movies_dir.mkdir(exist_ok=True)

    helper_path = movies_dir / "check_facing.py"
    helper_path.write_text(_FACING_HELPER_SCRIPT)

    if want("top10"):
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

    # The flexible variant (receptor taken from the trajectory itself) is the
    # important one — it actually shows receptor motion during sampling.
    # The static/rigid variant (receptor pinned to the reference) is only
    # meaningful when a reference structure is available for comparison.
    traj_movies = {}
    if replica_path is not None:
        variants = [("Flexible", True, "traj-flexible")]
        if HAS_REFERENCE:
            variants.append(("", False, "traj-rigid"))
        for suffix, flexible, thing in variants:
            label = f"_{suffix}" if suffix else ""
            movie = movies_dir / f"Trajectory_{replica_id}{label}.mp4"
            if want(thing):
                print(f"Rendering trajectory animation (replica {replica_id}{label})...")
                script = movies_dir / f"trajectory_{replica_id}{label.lower()}.cxc"
                frames = movies_dir / f"traj_frames_{replica_id}{label.lower()}"
                _write_trajectory_cxc(replica_path, script, frames, helper_path, flexible_receptor=flexible)
                subprocess.run([chimerax, "--offscreen", "--exit", str(script)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
                _stitch_frames(ffmpeg, frames, movie, framerate=30)
                traj_movies[suffix] = movie
            elif movie.exists():
                # Not regenerating this run, but reuse it for the combine step below.
                traj_movies[suffix] = movie

    if HAS_REFERENCE:
        # Tag the filename when other-replica context points are hidden, so this
        # doesn't silently overwrite the default (with-others) output under the
        # same name -- matches the "_SoloReplica" convention used for comparisons.
        solo_suffix = "_SoloReplica" if args.hide_other_replicas else ""
        for energy_type, plot_name, combined_prefix, ylabel, thing, combined_thing in (
            ("total",       "Total_Energy_vs_RMSD",       f"Combined_Trajectory_{replica_id}_Total_Energy_vs_RMSD",       "Total energy",       "rmsd-vs-tot-eng", "combined-tot-eng"),
            ("interaction", "Interaction_Energy_vs_RMSD", f"Combined_Trajectory_{replica_id}_Interaction_Energy_vs_RMSD", "Interaction energy", "rmsd-vs-int-eng", "combined-int-eng"),
        ):
            plot_movie = movies_dir / f"{plot_name}{solo_suffix}.mp4"
            # A requested combined output needs the plot too, even if the plot
            # itself wasn't explicitly requested -- generate it if missing.
            if want(thing) or (want(combined_thing) and not plot_movie.exists()):
                csv_files = sorted(WORK_DIR.glob(f"plots/E_RMSD_*_{energy_type}.csv"))
                if csv_files:
                    print(f"Generating E vs RMSD animation ({energy_type})...")
                    _make_energy_rmsd_animation(
                        csv_files[0], plot_movie, ylabel,
                        replica_id=replica_id, frames_per_replica=FRAMES_PER_REPLICA,
                        show_other_replicas=not args.hide_other_replicas,
                        other_replica_color=args.other_replica_color,
                        rmsd_log_scale=args.rmsd_log_scale,
                    )
            if not want(combined_thing) or not plot_movie.exists():
                continue
            for suffix, traj_movie in traj_movies.items():
                if not traj_movie.exists():
                    continue
                label = f"_{suffix}" if suffix else ""
                print(f"Combining trajectory{label} + {energy_type} energy plot...")
                _combine_movies(ffmpeg, traj_movie, plot_movie, movies_dir / f"{combined_prefix}{solo_suffix}{label}.mp4")

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
    reference_receptor_chains: List[str] | None = None,
    reference_peptide_chains: List[str] | None = None,
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
        reference_receptor_chains=reference_receptor_chains or receptor_chains,
        reference_peptide_chains=reference_peptide_chains or peptide_chains,
    )
    script_path = work_dir / "make_movies.py"
    script_path.write_text(header + _SCRIPT_BODY, encoding="utf-8")
    return script_path
