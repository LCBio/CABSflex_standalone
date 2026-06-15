"""
Regression test for bug: CABSflex failing on proteins containing MSE (selenomethionine)
residues with the error:
  "No sequential similarity between input and reference according to used alignment
   method (trivial)."

Root cause: MSE residues appear as HETATM records in PDB files. The trivial alignment
calls select("name CA and not HETERO"), which excluded MSE CA atoms from the reference
structure but not from the CABS trajectory template (which converts all residues to
standard AAs). The length mismatch caused trivial alignment to raise AlignError → ValueError.

Fix: MSE was added to extended_amino_acids.json, causing it to be re-labelled as MET and
marked hetatm=False during PDB loading (commit 61f7660).

Affected PDB entries: 2gzr, 2nvp, 2phz, 3r5s (all contain MSE HETATM records).
"""

import textwrap
import tempfile
import os
import pytest

from CABS.structures.pdblib import Pdb
from CABS.utils.align import TrivialAlign, AlignError, align_to


# Minimal PDB with a 5-residue chain A where residues 2 and 4 are MSE (HETATM).
# Only CA atoms are included to keep the fixture small.
_MSE_PDB_CONTENT = textwrap.dedent("""\
    ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  5.00           C
    HETATM    2  CA  MSE A   2       2.000   0.000   0.000  1.00  5.00           C
    ATOM      3  CA  GLY A   3       3.000   0.000   0.000  1.00  5.00           C
    HETATM    4  CA  MSE A   4       4.000   0.000   0.000  1.00  5.00           C
    ATOM      5  CA  LEU A   5       5.000   0.000   0.000  1.00  5.00           C
    END
""")


@pytest.fixture
def mse_pdb_path(tmp_path):
    pdb_file = tmp_path / "mse_test.pdb"
    pdb_file.write_text(_MSE_PDB_CONTENT)
    return str(pdb_file)


class TestMseLoading:
    """MSE residues must be treated as standard amino acids (MET) after loading."""

    def test_mse_not_marked_as_hetatm(self, mse_pdb_path):
        """After loading, no CA atom in the chain should be flagged as hetatm."""
        pdb = Pdb(mse_pdb_path)
        ca_atoms = pdb.atoms.select("name CA")
        hetatm_ca = [a for a in ca_atoms.atoms if a.hetatm]
        assert hetatm_ca == [], (
            f"MSE CA atoms still marked as HETATM: "
            f"{[(a.resname, a.resnum) for a in hetatm_ca]}"
        )

    def test_mse_converted_to_met(self, mse_pdb_path):
        """MSE residues must be re-labelled as MET."""
        pdb = Pdb(mse_pdb_path)
        for atom in pdb.atoms.select("name CA").atoms:
            assert atom.resname != "MSE", (
                f"Residue {atom.resnum} still named MSE; expected MET."
            )

    def test_ca_count_unaffected_by_hetero_filter(self, mse_pdb_path):
        """
        select('name CA') and select('name CA and not HETERO') must return
        the same count — this is the direct precondition for trivial alignment.
        """
        pdb = Pdb(mse_pdb_path)
        ca_total = pdb.atoms.select("name CA")
        ca_no_hetero = pdb.atoms.select("name CA and not HETERO")
        assert len(ca_total) == len(ca_no_hetero), (
            f"HETERO filter removes CA atoms: total={len(ca_total)}, "
            f"after filter={len(ca_no_hetero)}. "
            f"MSE residues are still being excluded."
        )

    def test_mse_ca_count_equals_five(self, mse_pdb_path):
        """All 5 CA atoms (including ex-MSE positions) must be present."""
        pdb = Pdb(mse_pdb_path)
        ca_atoms = pdb.atoms.select("name CA")
        assert len(ca_atoms) == 5, (
            f"Expected 5 CA atoms, got {len(ca_atoms)}. "
            f"MSE residues may have been dropped."
        )


class TestTrivialAlignmentWithMse:
    """
    Trivial alignment of a structure against itself must succeed even when
    the original structure contains MSE residues.
    This is the exact failure mode reported in the bug.
    """

    def test_trivial_align_self_succeeds(self, mse_pdb_path):
        """Aligning the MSE structure against itself (trivial) must not raise."""
        pdb = Pdb(mse_pdb_path)
        chains = list(pdb.atoms.list_chains().keys())
        # This call mirrors what job.py does at the trajectory-loading step.
        # Before the fix it raised:
        #   ValueError: No sequential similarity between input and reference
        #               according to used alignment method (trivial).
        ref_sstc, tmp_sstc, aln = align_to(
            pdb.atoms, chains, pdb.atoms, chains, align_mth="trivial"
        )
        assert len(aln) == 5, (
            f"Trivial alignment returned {len(aln)} pairs; expected 5."
        )

    def test_trivial_align_raises_on_size_mismatch(self, tmp_path):
        """
        Sanity-check: trivial alignment must still fail when structures genuinely
        differ in length (guard against over-permissive fixes).
        """
        short_pdb_content = textwrap.dedent("""\
            ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  5.00           C
            ATOM      2  CA  GLY A   2       2.000   0.000   0.000  1.00  5.00           C
            END
        """)
        short_pdb_file = tmp_path / "short.pdb"
        short_pdb_file.write_text(short_pdb_content)

        long_pdb_content = textwrap.dedent("""\
            ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00  5.00           C
            ATOM      2  CA  GLY A   2       2.000   0.000   0.000  1.00  5.00           C
            ATOM      3  CA  LEU A   3       3.000   0.000   0.000  1.00  5.00           C
            END
        """)
        long_pdb_file = tmp_path / "long.pdb"
        long_pdb_file.write_text(long_pdb_content)

        short_pdb = Pdb(str(short_pdb_file))
        long_pdb = Pdb(str(long_pdb_file))

        short_chains = list(short_pdb.atoms.list_chains().keys())
        long_chains = list(long_pdb.atoms.list_chains().keys())

        with pytest.raises(ValueError, match="No sequential similarity"):
            align_to(
                short_pdb.atoms, short_chains,
                long_pdb.atoms, long_chains,
                align_mth="trivial",
            )
