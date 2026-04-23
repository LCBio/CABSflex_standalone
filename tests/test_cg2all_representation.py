from pathlib import Path
from types import SimpleNamespace

from CABS.utils import utils
from CABS.reconstruction.cg2all import convert_cg_to_all, sync_residues, sync_residues_with_template
from CABS.reconstruction.cg2all_trajectory import (
    _write_single_model_cg2all_topology,
    reconstruct_job_outputs,
    reconstruct_trajectory,
)


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


def _write_multimodel_calpha_pdb(path: Path):
    path.write_text(
        "".join(
            [
                "MODEL        1\n",
                _pdb_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
                _pdb_line(2, "CA", "VAL", "A", 2, 3.8, 0.0, 0.0),
                "ENDMDL\n",
                "MODEL        2\n",
                _pdb_line(1, "CA", "ALA", "A", 1, 1.0, 1.0, 1.0),
                _pdb_line(2, "CA", "VAL", "A", 2, 4.8, 1.0, 1.0),
                "ENDMDL\n",
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

    convert_cg_to_all(str(input_pdb), work_dir=str(tmp_path))

    assert captured["command"][captured["command"].index("--cg") + 1] == "CalphaBasedModel"
    assert captured["input_pdb"].count(" CA ") == 3
    assert " SC " not in captured["input_pdb"]


def test_convert_cg_to_all_uses_calpha_sc_model(tmp_path, monkeypatch):
    input_pdb = tmp_path / "model.pdb"
    _write_calpha_pdb(input_pdb)
    (tmp_path / "output_pdbs").mkdir()
    captured = _capture_cg2all(monkeypatch)

    convert_cg_to_all(
        str(input_pdb),
        work_dir=str(tmp_path),
        cg2all_representation="calpha-sc",
    )

    assert captured["command"][captured["command"].index("--cg") + 1] == "CalphaSCModel"
    assert captured["input_pdb"].count(" CA ") == 3
    assert captured["input_pdb"].count(" SC ") == 3


def test_utils_writer_preserves_sc_for_cg2all_native_model_names(tmp_path):
    input_pdb = tmp_path / "cluster_0.pdb"
    output_pdb = tmp_path / "cluster_0_input_for_cg2all.pdb"
    _write_calpha_pdb(input_pdb)

    with output_pdb.open("w") as handle:
        utils._write_cg2all_input_pdb(str(input_pdb), handle, "CalphaSCModel")

    written = output_pdb.read_text()

    assert written.count(" CA ") == 3
    assert written.count(" SC ") == 3


def test_batch_topology_writer_uses_only_first_model(tmp_path):
    input_pdb = tmp_path / "cluster_0.pdb"
    output_pdb = tmp_path / "cluster_0_temp_top.pdb"
    _write_multimodel_calpha_pdb(input_pdb)

    with output_pdb.open("w") as handle:
        _write_single_model_cg2all_topology(str(input_pdb), handle, "CalphaSCModel")

    written = output_pdb.read_text()

    assert "MODEL" not in written
    assert "ENDMDL" not in written
    assert written.count(" CA ") == 2
    assert written.count(" SC ") == 2


def test_sync_residues_preserves_insertion_codes_for_multimodel_outputs(tmp_path):
    reference = tmp_path / "replica.pdb"
    reconstructed = tmp_path / "replica_all_atom.pdb"

    reference.write_text(
        "".join(
            [
                "MODEL        1\n",
                _pdb_line(1, "CA", "LYS", "H", 60, 0.0, 0.0, 0.0)[:26] + "A" + _pdb_line(1, "CA", "LYS", "H", 60, 0.0, 0.0, 0.0)[27:],
                _pdb_line(2, "CA", "ILE", "H", 60, 1.0, 0.0, 0.0)[:26] + "B" + _pdb_line(2, "CA", "ILE", "H", 60, 1.0, 0.0, 0.0)[27:],
                "ENDMDL\n",
                "MODEL        2\n",
                _pdb_line(1, "CA", "LYS", "H", 60, 0.5, 0.5, 0.5)[:26] + "A" + _pdb_line(1, "CA", "LYS", "H", 60, 0.5, 0.5, 0.5)[27:],
                _pdb_line(2, "CA", "ILE", "H", 60, 1.5, 0.5, 0.5)[:26] + "B" + _pdb_line(2, "CA", "ILE", "H", 60, 1.5, 0.5, 0.5)[27:],
                "ENDMDL\n",
            ]
        )
    )

    reconstructed.write_text(
        "".join(
            [
                "MODEL        1\n",
                _pdb_line(1, "N", "LYS", "A", 60, 0.0, 0.0, 0.0),
                _pdb_line(2, "CA", "LYS", "A", 60, 0.1, 0.0, 0.0),
                _pdb_line(3, "N", "ILE", "A", 60, 1.0, 0.0, 0.0),
                _pdb_line(4, "CA", "ILE", "A", 60, 1.1, 0.0, 0.0),
                "ENDMDL\n",
                "MODEL        2\n",
                _pdb_line(1, "N", "LYS", "A", 60, 0.5, 0.5, 0.5),
                _pdb_line(2, "CA", "LYS", "A", 60, 0.6, 0.5, 0.5),
                _pdb_line(3, "N", "ILE", "A", 60, 1.5, 0.5, 0.5),
                _pdb_line(4, "CA", "ILE", "A", 60, 1.6, 0.5, 0.5),
                "ENDMDL\n",
            ]
        )
    )

    sync_residues(reference, reconstructed)
    synced = reconstructed.read_text()

    assert "LYS H  60A" in synced
    assert "ILE H  60B" in synced
    assert " N   ILE H  60A" not in synced


def test_reconstruct_job_outputs_ignores_existing_all_atom_pdbs(tmp_path, monkeypatch):
    output_dir = tmp_path / "output_pdbs"
    output_dir.mkdir()

    for name in [
        "cluster_0.pdb",
        "cluster_0_all_atom.pdb",
        "cluster_1.pdb",
        "replica_0.pdb",
        "replica_0_all_atom.pdb",
        "model_0.pdb",
    ]:
        (output_dir / name).write_text("")

    captured = []

    def fake_reconstruct_trajectory(**kwargs):
        captured.append(Path(kwargs["trajectory_file"]).name)

    monkeypatch.setattr(
        "CABS.reconstruction.cg2all_trajectory.reconstruct_trajectory",
        fake_reconstruct_trajectory,
    )

    job = SimpleNamespace(
        aa_method="cg2all",
        aa_rebuild="A",
        cg2all_representation="calpha-sc",
        cg2all_env_prefix=None,
        aa_rebuild_workers=1,
        batch_size=1,
        renumber=False,
    )

    reconstruct_job_outputs(job, str(output_dir))

    assert captured == ["cluster_0.pdb", "cluster_1.pdb", "replica_0.pdb"]


def test_reconstruct_trajectory_uses_last_frame_topology_helper_and_cleans_up(tmp_path, monkeypatch):
    trajectory_file = tmp_path / "replica.pdb"
    trajectory_file.write_text(
        "".join(
            [
                "MODEL        1\n",
                _pdb_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
                "ENDMDL\n",
                "MODEL        2\n",
                _pdb_line(1, "CA", "ALA", "A", 1, 1.0, 1.0, 1.0),
                "ENDMDL\n",
            ]
        )
    )

    output_dcd = tmp_path / "replica_all_atom.dcd"
    output_pdb = tmp_path / "replica_all_atom.pdb"
    captured = {}

    class FakeInputTraj:
        n_frames = 2

        def save_dcd(self, path):
            Path(path).write_text("temp dcd")

    class FakeOutputTraj:
        def save_pdb(self, path):
            Path(path).write_text(
                "MODEL        1\n"
                + _pdb_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0)
                + _pdb_line(2, "CA", "ALA", "A", 1, 0.1, 0.0, 0.0)
                + "ENDMDL\n"
            )

    def fake_md_load(path, top=None):
        if top is None:
            return FakeInputTraj()
        captured["mdtraj_topology"] = top
        return FakeOutputTraj()

    def fake_run(command_parts, **kwargs):
        captured["command"] = command_parts
        Path(command_parts[command_parts.index("-o") + 1]).write_text("aa dcd")
        helper_pdb = Path(command_parts[command_parts.index("-opdb") + 1])
        helper_pdb.write_text(_pdb_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0))
        return SimpleNamespace(stdout="ok")

    monkeypatch.setattr("CABS.reconstruction.cg2all_trajectory.md.load", fake_md_load)
    monkeypatch.setattr("CABS.reconstruction.cg2all_trajectory.subprocess.run", fake_run)

    reconstruct_trajectory(
        topology_pdb=str(trajectory_file),
        trajectory_file=str(trajectory_file),
        output_dcd=str(output_dcd),
        output_pdb=str(output_pdb),
        cg_model="CalphaSCModel",
    )

    assert captured["command"][captured["command"].index("-o") + 1] == str(output_dcd)
    assert captured["command"][captured["command"].index("-opdb") + 1].endswith("_last_frame.pdb")
    assert captured["mdtraj_topology"].endswith("_last_frame.pdb")
    assert output_pdb.exists()
    assert not output_dcd.exists()
    assert not Path(captured["mdtraj_topology"]).exists()


def test_sync_residues_with_template_handles_flattened_output_blocks(tmp_path):
    reference = tmp_path / "replica.pdb"
    topology = tmp_path / "replica_all_atom_last_frame.pdb"
    output = tmp_path / "replica_all_atom.pdb"

    reference.write_text(
        "".join(
            [
                "MODEL        1\n",
                _pdb_line(1, "CA", "GLN", "H", 170, 0.0, 0.0, 0.0),
                _pdb_line(2, "CA", "SER", "H", 171, 1.0, 0.0, 0.0),
                "ENDMDL\n",
                "MODEL        2\n",
                _pdb_line(1, "CA", "GLN", "H", 170, 0.5, 0.5, 0.5),
                _pdb_line(2, "CA", "SER", "H", 171, 1.5, 0.5, 0.5),
                "ENDMDL\n",
            ]
        )
    )

    topology.write_text(
        "".join(
            [
                _pdb_line(1, "N", "GLN", "A", 1, 0.0, 0.0, 0.0),
                _pdb_line(2, "CA", "GLN", "A", 1, 0.1, 0.0, 0.0),
                _pdb_line(3, "C", "GLN", "A", 1, 0.2, 0.0, 0.0),
                _pdb_line(4, "N", "SER", "A", 2, 1.0, 0.0, 0.0),
                _pdb_line(5, "CA", "SER", "A", 2, 1.1, 0.0, 0.0),
                _pdb_line(6, "C", "SER", "A", 2, 1.2, 0.0, 0.0),
            ]
        )
    )

    output.write_text(
        "".join(
            [
                "MODEL        1\n",
                _pdb_line(1, "N", "GLN", "A", 170, 0.0, 0.0, 0.0),
                _pdb_line(2, "CA", "GLN", "A", 170, 0.1, 0.0, 0.0),
                _pdb_line(3, "C", "GLN", "A", 170, 0.2, 0.0, 0.0),
                _pdb_line(4, "N", "GLN", "A", 170, 1.0, 0.0, 0.0),
                _pdb_line(5, "CA", "SER", "A", 170, 1.1, 0.0, 0.0),
                _pdb_line(6, "C", "SER", "A", 170, 1.2, 0.0, 0.0),
                "ENDMDL\n",
                "MODEL        2\n",
                _pdb_line(1, "N", "GLN", "A", 170, 0.5, 0.5, 0.5),
                _pdb_line(2, "CA", "GLN", "A", 170, 0.6, 0.5, 0.5),
                _pdb_line(3, "C", "GLN", "A", 170, 0.7, 0.5, 0.5),
                _pdb_line(4, "N", "GLN", "A", 170, 1.5, 0.5, 0.5),
                _pdb_line(5, "CA", "SER", "A", 170, 1.6, 0.5, 0.5),
                _pdb_line(6, "C", "SER", "A", 170, 1.7, 0.5, 0.5),
                "ENDMDL\n",
            ]
        )
    )

    sync_residues_with_template(reference, topology, output)
    synced = output.read_text()

    assert "N   GLN H 170" in synced
    assert "CA  SER H 171" in synced
    assert synced.count(" H 170") == 6
    assert synced.count(" H 171") == 6
