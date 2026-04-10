import pytest

from CABS.core.job import DockTask
from CABS.structures.atom import Atom, Atoms
from CABS.structures.vector3d import Vector3d


AA3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}

NATIVE_1D4T_PEPTIDE = "KSLTIYAQVQK"


def make_chain(sequence, chid):
    atoms = []
    for index, aa in enumerate(sequence, start=1):
        atoms.append(
            Atom(
                model=0,
                name="CA",
                resname=AA3[aa],
                chid=chid,
                resnum=index,
                coord=Vector3d(float(index), 0.0, 0.0),
                occ=1.0,
                bfac=0.0,
                hetatm=False,
            )
        )
    return Atoms(atoms)


class FailingStrictTrajectory:
    def __init__(self, template):
        self.template = template

    def align_to(self, *args, **kwargs):
        raise ValueError("No sequential similarity between input and reference.")


class ExactStrictTrajectory:
    def __init__(self, template):
        self.template = template

    def align_to(self, ref_stc, ref_chain, model_chain, align_mth="SW", kwargs=None):
        ref_atoms = ref_stc.select(f"name CA and chain {ref_chain}")
        model_atoms = self.template.select(f"name CA and chain {model_chain}")
        pairs = tuple(zip(ref_atoms, model_atoms))
        return Atoms([a for a, _ in pairs]), Atoms([b for _, b in pairs]), pairs


def make_task(reference_atoms, model_atoms, strict_trajectory):
    task = DockTask.__new__(DockTask)
    task.reference = (reference_atoms, "A", "B")
    task.trajectory = strict_trajectory(model_atoms)
    task.align = "SW"
    task.align_peptide_options = {}
    task.peptide_rmsd_mode = "auto"
    task.peptide_min_aligned_length = 5
    task.peptide_min_aligned_fraction = 0.5
    task.peptide_min_identity = 0.5
    task.peptide_max_gap_fraction = 0.4
    task.peptide_min_contiguous_block = 3
    return task


def test_strict_peptide_alignment_selected_for_identical_sequences():
    reference_atoms = make_chain("GPPPAMPARPT", "B")
    model_atoms = make_chain("GPPPAMPARPT", "B")
    task = make_task(reference_atoms, model_atoms, ExactStrictTrajectory)

    selection = task._select_peptide_alignment("B", "B")

    assert selection["mode"] == "strict"
    assert selection["metrics"]["aligned_length"] == 11
    assert selection["metrics"]["mismatch_count"] == 0


def test_overlap_peptide_alignment_selected_for_flanking_overhang():
    reference_atoms = make_chain("GPPPAMPARPT", "B")
    model_atoms = make_chain("EGPPPAMPARPT", "B")
    task = make_task(reference_atoms, model_atoms, FailingStrictTrajectory)

    selection = task._select_peptide_alignment("B", "B")

    assert selection["mode"] == "overlap"
    assert selection["metrics"]["aligned_length"] == 11
    assert selection["metrics"]["excluded_model_residues"] == 1
    assert selection["results"]["strict"]["status"] == "failed"


def test_mutational_peptide_alignment_selected_when_exact_overlap_is_too_fragmented():
    reference_atoms = make_chain("GPGPGP", "B")
    model_atoms = make_chain("GAGAGA", "B")
    task = make_task(reference_atoms, model_atoms, FailingStrictTrajectory)

    selection = task._select_peptide_alignment("B", "B")

    assert selection["mode"] == "mutational"
    assert selection["metrics"]["aligned_length"] >= 5
    assert selection["metrics"]["identity"] >= 0.5
    assert selection["results"]["overlap"]["status"] == "rejected"


@pytest.mark.parametrize(
    ("variant_sequence", "strict_trajectory", "expected_mode"),
    [
        (NATIVE_1D4T_PEPTIDE, ExactStrictTrajectory, "strict"),
        ("Q" + NATIVE_1D4T_PEPTIDE, FailingStrictTrajectory, "overlap"),
        ("KALTVYVQIQR", FailingStrictTrajectory, "mutational"),
        ("RALTVFVIVQR", FailingStrictTrajectory, None),
    ],
)
def test_1d4t_peptide_panel_covers_strict_overlap_mutational_and_rejection(
    variant_sequence, strict_trajectory, expected_mode
):
    """
    Sequence-only regression panel inspired by 1D4T (SAP chain A with the SLAM peptide).

    Native peptide resolved in the crystal structure:
    KSLTIYAQVQK

    Variants cover:
    - strict: exact native peptide
    - overlap: flanking overhang that should align to the resolved core
    - mutational: screening-like variant with several substitutions
    - rejection: variant too dissimilar to pass current thresholds
    """
    reference_atoms = make_chain(NATIVE_1D4T_PEPTIDE, "B")
    model_atoms = make_chain(variant_sequence, "B")
    task = make_task(reference_atoms, model_atoms, strict_trajectory)

    selection = task._select_peptide_alignment("B", "B")

    assert selection["mode"] == expected_mode
    if expected_mode == "strict":
        assert selection["metrics"]["aligned_length"] == len(NATIVE_1D4T_PEPTIDE)
    elif expected_mode == "overlap":
        assert selection["metrics"]["excluded_model_residues"] > 0
        assert selection["metrics"]["identity"] == 1.0
    elif expected_mode == "mutational":
        assert selection["metrics"]["identity"] >= task.peptide_min_identity
        assert selection["results"]["overlap"]["status"] == "rejected"
    else:
        assert selection["metrics"] is None
