from pathlib import Path
from types import SimpleNamespace

from CABS.utils import utils


def _pdb_line(serial, name, resname, chain, resnum, x, y, z):
    fmt_name = f" {name:<3s}"
    return (
        f"ATOM  {serial:5d} {fmt_name:4s} {resname:<4s}{chain:1s}{resnum:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}\n"
    )


def _write_calpha_pdb(path: Path):
    path.write_text(
        "".join(
            [
                _pdb_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
                _pdb_line(2, "CA", "VAL", "A", 2, 3.8, 0.0, 0.0),
                _pdb_line(3, "CA", "GLY", "A", 3, 7.6, 0.0, 0.0),
            ]
        )
    )


def _capture_cg2all(monkeypatch):
    captured = {}

    def fake_run(command_parts, **kwargs):
        input_path = Path(command_parts[command_parts.index("-p") + 1])
        captured["command"] = command_parts
        captured["input_pdb"] = input_path.read_text()
        return SimpleNamespace(stdout="ok")

    monkeypatch.setattr(utils.subprocess, "run", fake_run)
    return captured


def test_convert_cg_to_all_uses_calpha_model(tmp_path, monkeypatch):
    input_pdb = tmp_path / "model.pdb"
    _write_calpha_pdb(input_pdb)
    (tmp_path / "output_pdbs").mkdir()
    captured = _capture_cg2all(monkeypatch)

    utils.convert_cg_to_all(str(input_pdb), work_dir=str(tmp_path))

    assert captured["command"][captured["command"].index("--cg") + 1] == "CalphaBasedModel"
    assert captured["input_pdb"].count(" CA ") == 3
    assert " SC " not in captured["input_pdb"]


def test_convert_cg_to_all_uses_calpha_sc_model(tmp_path, monkeypatch):
    input_pdb = tmp_path / "model.pdb"
    _write_calpha_pdb(input_pdb)
    (tmp_path / "output_pdbs").mkdir()
    captured = _capture_cg2all(monkeypatch)

    utils.convert_cg_to_all(
        str(input_pdb),
        work_dir=str(tmp_path),
        cg2all_representation="calpha-sc",
    )

    assert captured["command"][captured["command"].index("--cg") + 1] == "CalphaSCModel"
    assert captured["input_pdb"].count(" CA ") == 3
    assert captured["input_pdb"].count(" SC ") == 3
