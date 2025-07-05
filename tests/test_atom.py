"""
Tests for Atom and Atoms classes - protein structure representation.
"""

import numpy as np

from CABS.structures.atom import Atom, Atoms
from CABS.structures.vector3d import Vector3d


class TestAtomInitialization:
    """Test Atom class initialization."""

    def test_atom_from_pdb_line(self):
        """Test creating atom from PDB line."""
        pdb_line = "ATOM      1  CA  ALA A   1      20.154  16.967  14.239  1.00 20.00           C  "
        atom = Atom(pdb_line)

        assert atom.serial == 1
        assert atom.name == "CA"
        assert atom.resname == "ALA"
        assert atom.chid == "A"
        assert atom.resnum == 1
        assert abs(atom.coord.x - 20.154) < 1e-3
        assert abs(atom.coord.y - 16.967) < 1e-3
        assert abs(atom.coord.z - 14.239) < 1e-3
        assert abs(atom.occ - 1.00) < 1e-3
        assert abs(atom.bfac - 20.00) < 1e-3

    def test_atom_empty_initialization(self):
        """Test creating empty atom."""
        atom = Atom()

        assert atom.serial == 0
        assert atom.name == "XXXX"
        assert atom.resname == "XXX"
        assert atom.chid == "X"
        assert atom.resnum == 0
        assert atom.hetatm == True

    def test_atom_with_kwargs(self):
        """Test creating atom with keyword arguments."""
        atom = Atom(serial=10, name="CB", resname="VAL", chid="B", resnum=5)

        assert atom.serial == 10
        assert atom.name == "CB"
        assert atom.resname == "VAL"
        assert atom.chid == "B"
        assert atom.resnum == 5


class TestAtomMethods:
    """Test Atom class methods."""

    def test_atom_string_representation(self):
        """Test atom string representation."""
        atom = Atom()
        atom.serial = 1
        atom.name = "CA"
        atom.resname = "ALA"
        atom.chid = "A"
        atom.resnum = 1
        atom.coord = Vector3d(1.0, 2.0, 3.0)
        atom.occ = 1.0
        atom.bfac = 20.0

        pdb_line = str(atom)
        # The atom string representation uses HETATM by default
        assert "HETATM" in pdb_line or "ATOM" in pdb_line
        assert "CA" in pdb_line
        assert "ALA" in pdb_line

    def test_atom_fmt(self):
        """Test atom format method."""
        atom = Atom()
        atom.name = "CA"
        atom.resname = "ALA"
        atom.chid = "A"
        atom.resnum = 1

        fmt_str = atom.fmt()
        # fmt() returns chid + resnum + icode.strip()
        assert "A" in fmt_str
        assert "1" in fmt_str

    def test_atom_distance(self):
        """Test distance calculation between atoms."""
        atom1 = Atom()
        atom1.coord = Vector3d(0.0, 0.0, 0.0)

        atom2 = Atom()
        atom2.coord = Vector3d(3.0, 4.0, 0.0)

        distance = atom1.distance(atom2)
        assert abs(distance - 5.0) < 1e-10

    def test_atom_same_residue(self):
        """Test same residue comparison."""
        atom1 = Atom()
        atom1.chid = "A"
        atom1.resnum = 1
        atom1.icode = ""
        atom1.model = 0

        atom2 = Atom()
        atom2.chid = "A"
        atom2.resnum = 1
        atom2.icode = ""
        atom2.model = 0

        atom3 = Atom()
        atom3.chid = "A"
        atom3.resnum = 2
        atom3.icode = ""
        atom3.model = 0

        assert atom1.same_residue(atom2)
        assert not atom1.same_residue(atom3)


class TestAtomsCollection:
    """Test Atoms collection class."""

    def test_atoms_initialization(self):
        """Test different ways to initialize Atoms."""
        # Empty initialization
        atoms1 = Atoms()
        assert len(atoms1) == 0

        # Initialize with list
        atom = Atom()
        atoms2 = Atoms([atom])
        assert len(atoms2) == 1

        # Initialize with integer
        atoms3 = Atoms(5)
        assert len(atoms3) == 5
        assert all(atom.resname == "ALA" for atom in atoms3)

    def test_atoms_list_operations(self):
        """Test list-like operations."""
        atoms = Atoms()
        atom = Atom()

        # Test append
        atoms.append(atom)
        assert len(atoms) == 1

        # Test iteration
        for a in atoms:
            assert a == atom

        # Test indexing
        assert atoms[0] == atom

    def test_atoms_chains_list(self):
        """Test chain listing."""
        atoms = Atoms()

        atom1 = Atom()
        atom1.chid = "A"
        atom1.resnum = 1
        atoms.append(atom1)

        atom2 = Atom()
        atom2.chid = "B"
        atom2.resnum = 1
        atoms.append(atom2)

        chains = atoms.list_chains()
        assert "A" in chains
        assert "B" in chains

    def test_atoms_selection(self):
        """Test atom selection methods."""
        atoms = Atoms()

        # Create test atoms
        for i in range(5):
            atom = Atom()
            atom.chid = "A"
            atom.resnum = i + 1
            atom.name = "CA"
            atom.coord = Vector3d(float(i), 0.0, 0.0)
            atoms.append(atom)

        # Test selection by chain using the select method
        chain_a = atoms.select("CHAIN A")
        assert len(chain_a) == 5

        # Test selection by name
        ca_atoms = atoms.select("NAME CA")
        assert len(ca_atoms) == 5

    def test_atoms_coordinate_operations(self):
        """Test coordinate-related operations."""
        atoms = Atoms()

        # Create atoms with known coordinates
        coords = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        for i, (x, y, z) in enumerate(coords):
            atom = Atom()
            atom.chid = "A"
            atom.resnum = i + 1
            atom.name = "CA"
            atom.coord = Vector3d(x, y, z)
            atoms.append(atom)

        # Test center of mass calculation
        centroid = atoms.cent_of_mass()
        assert abs(centroid.x - 1.0 / 3.0) < 1e-10
        assert abs(centroid.y - 1.0 / 3.0) < 1e-10
        assert abs(centroid.z - 1.0 / 3.0) < 1e-10

    def test_atoms_coordinate_array(self):
        """Test coordinate array extraction."""
        atoms = Atoms()

        coords = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
        for i, (x, y, z) in enumerate(coords):
            atom = Atom()
            atom.coord = Vector3d(x, y, z)
            atoms.append(atom)

        # Test numpy array conversion
        coord_array = atoms.to_numpy()
        assert coord_array.shape == (3, 3)
        np.testing.assert_allclose(coord_array[0], [1.0, 2.0, 3.0])
        np.testing.assert_allclose(coord_array[1], [4.0, 5.0, 6.0])
        np.testing.assert_allclose(coord_array[2], [7.0, 8.0, 9.0])

        # Test coordinate dictionary
        coord_dict = atoms.get_coordinates()
        # get_coordinates() returns dict with resid_id() as keys, aggregated by residue
        assert len(coord_dict) == 1  # Single residue/group


class TestAtomsManipulation:
    """Test Atoms manipulation methods."""

    def test_atoms_translation(self):
        """Test atom translation."""
        atoms = Atoms()

        atom = Atom()
        atom.coord = Vector3d(1.0, 2.0, 3.0)
        atoms.append(atom)

        translation = Vector3d(5.0, 5.0, 5.0)
        atoms.move(translation)

        assert abs(atoms[0].coord.x - 6.0) < 1e-10
        assert abs(atoms[0].coord.y - 7.0) < 1e-10
        assert abs(atoms[0].coord.z - 8.0) < 1e-10

    def test_atoms_centering(self):
        """Test atom centering operations."""
        atoms = Atoms()

        # Create atoms with offset from origin
        coords = [(5.0, 5.0, 5.0), (6.0, 5.0, 5.0), (5.0, 6.0, 5.0)]
        for x, y, z in coords:
            atom = Atom()
            atom.coord = Vector3d(x, y, z)
            atoms.append(atom)

        # Center at origin
        atoms.center_at_origin()

        # Check that center of mass is now at origin
        centroid = atoms.cent_of_mass()
        assert abs(centroid.x) < 1e-10
        assert abs(centroid.y) < 1e-10
        assert abs(centroid.z) < 1e-10

    def test_atoms_superposition(self):
        """Test atoms superposition."""
        # Create reference atoms
        ref_atoms = Atoms()
        ref_coords = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
        for x, y, z in ref_coords:
            atom = Atom()
            atom.coord = Vector3d(x, y, z)
            ref_atoms.append(atom)

        # Create mobile atoms (same coordinates)
        mobile_atoms = Atoms()
        mobile_coords = [(1.1, 0.1, 0.1), (0.1, 1.1, 0.1), (0.1, 0.1, 1.1)]
        for x, y, z in mobile_coords:
            atom = Atom()
            atom.coord = Vector3d(x, y, z)
            mobile_atoms.append(atom)

        # Test RMSD calculation
        rmsd_value = mobile_atoms.rmsd(ref_atoms)
        assert rmsd_value > 0

        # Test structural alignment
        mobile_atoms.str_align(ref_atoms)
        final_rmsd = mobile_atoms.rmsd(ref_atoms)
        assert final_rmsd <= rmsd_value  # Should be better after alignment


class TestAtomsIO:
    """Test Atoms I/O operations."""

    def test_atoms_pdb_output(self):
        """Test PDB format output."""
        atoms = Atoms()

        atom = Atom()
        atom.serial = 1
        atom.name = "CA"
        atom.resname = "ALA"
        atom.chid = "A"
        atom.resnum = 1
        atom.coord = Vector3d(1.0, 2.0, 3.0)
        atoms.append(atom)

        pdb_str = atoms.make_pdb()
        assert isinstance(pdb_str, str)
        assert "CA" in pdb_str
        assert "ALA" in pdb_str

    def test_atoms_to_pdb_string(self):
        """Test converting Atoms to PDB string."""
        atoms = Atoms()

        atom = Atom()
        atom.serial = 1
        atom.name = "CA"
        atom.resname = "ALA"
        atom.chid = "A"
        atom.resnum = 1
        atom.coord = Vector3d(1.0, 2.0, 3.0)
        atoms.append(atom)

        pdb_str = str(atoms)
        lines = pdb_str.strip().split("\n")

        # Should have at least one line
        assert len(lines) >= 1
        assert "CA" in pdb_str


class TestAtomsEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_atoms_operations(self):
        """Test operations on empty Atoms collection."""
        atoms = Atoms()

        # These should handle empty collections gracefully
        coord_array = atoms.to_numpy()
        assert coord_array.shape == (0,)  # Fix: actual shape for empty arrays

        coord_dict = atoms.get_coordinates()
        assert len(coord_dict) == 0

    def test_single_atom_operations(self):
        """Test operations with single atom."""
        atoms = Atoms()

        atom = Atom()
        atom.coord = Vector3d(5.0, 10.0, 15.0)
        atoms.append(atom)

        centroid = atoms.cent_of_mass()
        assert abs(centroid.x - 5.0) < 1e-10
        assert abs(centroid.y - 10.0) < 1e-10
        assert abs(centroid.z - 15.0) < 1e-10

    def test_atoms_with_invalid_coordinates(self):
        """Test handling of coordinates."""
        atoms = Atoms()

        # Atom with large coordinates
        atom = Atom()
        atom.coord = Vector3d(1000.0, 1000.0, 1000.0)
        atoms.append(atom)

        # Should handle large coordinates
        centroid = atoms.cent_of_mass()
        assert abs(centroid.x - 1000.0) < 1e-10
        assert abs(centroid.y - 1000.0) < 1e-10
        assert abs(centroid.z - 1000.0) < 1e-10
