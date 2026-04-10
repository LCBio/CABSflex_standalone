"""Generate ChimeraX (.cxc) scripts for CABS simulations.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, List, Optional


TEMPLATE_COLOR_RMSF = """# Color by RMSF
open model_*.pdb
set bgColor white
color byattribute a:bfactor
lighting intensity 0.7
graphics silhouettes true
"""

TEMPLATE_COLOR_BY_CHAIN = """# Color by chain
open model_*.pdb
set bgColor white
color bychain
lighting intensity 0.7
graphics silhouettes true
"""

TEMPLATE_COLOR_SS = """# Color by secondary structure
open model_*.pdb
set bgColor white
select helix
color sel #6b1a56ff
select strand
color sel dark goldenrod
select coil
color sel #cfcfcfff
select clear
lighting intensity 0.6
graphics silhouettes true
"""

TEMPLATE_RMSF_WORM = """# Color by RMSF and display as worm
open model_*.pdb
set bgColor white
color byattribute a:bfactor
Worm bfactor
lighting intensity 0.7
graphics silhouettes true
"""

TEMPLATE_RECORD_MOVIE = """# Record animation
movie record size 3840,2160
turn y 2 180; wait 180
movie stop
movie encode output my_movie.mp4 quality highest
"""


DEFAULT_PRESETS = {
    "color_rmsf": TEMPLATE_COLOR_RMSF,
    "color_chain": TEMPLATE_COLOR_BY_CHAIN,
    "color_ss": TEMPLATE_COLOR_SS,
    "rmsf_worm": TEMPLATE_RMSF_WORM,
    "record_movie": TEMPLATE_RECORD_MOVIE,
}


def write_cxc(out_dir: Path, name: str, content: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.cxc"
    with open(path, "w") as fh:
        fh.write(content)
    return path


def _generate_restraints_script(scripts_dir: Path, restraints_file: str, start_pdb_path: Optional[str] = None) -> Optional[Path]:
    # Simple restraints translator: create ChimeraX distance commands between CA atoms
    restraints_lines: List[str] = []
    try:
        with open(restraints_file, "r") as fh:
            for idx, line in enumerate(fh):
                parts = line.split()
                if len(parts) >= 2:
                    rc1 = parts[0].split(":")
                    rc2 = parts[1].split(":")
                    if len(rc1) == 2 and len(rc2) == 2:
                        res1, chain1 = rc1
                        res2, chain2 = rc2
                        # ChimeraX distance command between CA atoms
                        cmd = f"distance rest_{idx} #0:{chain1}:{res1}@CA #0:{chain2}:{res2}@CA"
                        restraints_lines.append(cmd)
    except Exception:
        return None

    if not restraints_lines:
        return None

    content_lines = []
    if start_pdb_path:
        rel = os.path.relpath(start_pdb_path, scripts_dir)
        content_lines.append(f"open {rel}")
    content_lines.extend(restraints_lines)
    content = "\n".join(content_lines) + "\n"
    return write_cxc(scripts_dir, "load_restraints", content)


def generate_chimerax_scripts(work_dir: str | Path, start_pdb_path: Optional[str], models_pdbs: List[str], restraints_file: Optional[str] = None, presets: Optional[Iterable[str]] = None) -> List[Path]:
    """Generate ChimeraX `.cxc` scripts in `work_dir`.
    """
    work_dir = Path(work_dir)
    if presets is None:
        presets = list(DEFAULT_PRESETS.keys())
    model_pattern = "model_*.pdb"
    if models_pdbs:
        model_pattern = os.path.relpath(str(Path(models_pdbs[0]).parent / "model_*.pdb"), work_dir)

    generated: List[Path] = []
    for preset in presets:
        tpl = DEFAULT_PRESETS.get(preset)
        if not tpl:
            continue
        generated.append(write_cxc(work_dir, preset, tpl.replace("model_*.pdb", model_pattern)))

    if restraints_file and os.path.exists(restraints_file):
        r = _generate_restraints_script(work_dir, restraints_file, start_pdb_path)
        if r:
            generated.append(r)

    return generated


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="chimerax_generator")
    p.add_argument("--sim-dir", required=True, help="Simulation directory to scan for PDBs (not used for wildcard scripts but recommended)")
    p.add_argument("--out-dir", default="./chimerax_scripts", help="Output directory for .cxc files")
    p.add_argument("--presets", help="Comma-separated presets to generate (color_rmsf,color_chain,color_ss,rmsf_worm,record_movie)")
    p.add_argument("--restraints", help="Path to restraints.txt to generate distances")
    return p


def main(argv: Optional[List[str]] = None) -> None:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    presets = None
    if args.presets:
        presets = [s.strip() for s in args.presets.split(",") if s.strip()]
    generated = generate_chimerax_scripts(args.out_dir, None, [], args.restraints, presets)
    for p in generated:
        print(p)


if __name__ == "__main__":
    main()
