#!/usr/bin/env python3
"""
Simple, working test suite for CABS that matches the actual implementation.
"""

import numpy as np

from CABS.structures.atom import Atom, Atoms
from CABS.structures.vector3d import Vector3d
from CABS.utils.utils import rmsd


def test_vector3d_basic_operations():
    """Test basic Vector3d operations."""
    v1 = Vector3d(1, 2, 3)
    v2 = Vector3d(4, 5, 6)

    # Test addition
    result = v1 + v2
    assert result.x == 5
    assert result.y == 7
    assert result.z == 9

    # Test length
    v = Vector3d(3, 4, 0)
    assert abs(v.length() - 5.0) < 1e-10


def test_atom_creation():
    """Test Atom creation and basic properties."""
    atom = Atom()
    atom.serial = 1
    atom.name = "CA"
    atom.resname = "ALA"
    atom.coord = Vector3d(1.0, 2.0, 3.0)

    assert atom.serial == 1
    assert atom.name == "CA"
    assert atom.resname == "ALA"
    assert atom.coord.x == 1.0


def test_atoms_collection():
    """Test Atoms collection basic functionality."""
    atoms = Atoms()

    # Test empty collection
    assert len(atoms) == 0

    # Add atoms
    atom1 = Atom()
    atom1.coord = Vector3d(1.0, 0.0, 0.0)
    atoms.append(atom1)

    atom2 = Atom()
    atom2.coord = Vector3d(0.0, 1.0, 0.0)
    atoms.append(atom2)

    assert len(atoms) == 2


def test_atoms_get_coordinates():
    """Test the actual get_coordinates method."""
    atoms = Atoms()

    atom = Atom()
    atom.serial = 1
    atom.resnum = 1
    atom.chid = "A"
    atom.coord = Vector3d(1.0, 2.0, 3.0)
    atoms.append(atom)

    coords = atoms.get_coordinates()
    assert isinstance(coords, dict)
    # The key format is resnum:chid
    key = list(coords.keys())[0]
    assert isinstance(coords[key], Vector3d)


def test_rmsd_calculation():
    """Test RMSD calculation with numpy arrays."""
    # Create simple test data
    coords1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    coords2 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    result = rmsd(coords1, coords2)
    assert result == 0.0


def test_atoms_string_representation():
    """Test string representation of atoms."""
    atoms = Atoms()
    atom = Atom()
    atom.serial = 1
    atom.name = "CA"
    atoms.append(atom)

    str_repr = str(atoms)
    assert isinstance(str_repr, str)
    assert len(str_repr) > 0


def test_atoms_sequence_creation():
    """Test creating Atoms from sequence string."""
    atoms = Atoms("ALA")
    assert len(atoms) == 3
    assert atoms[0].resname == "ALA"
    assert atoms[1].resname == "LEU"
    assert atoms[2].resname == "ALA"


def test_vector3d_edge_cases():
    """Test Vector3d with small numbers."""
    tiny = Vector3d(1e-15, 1e-15, 1e-15)
    length = tiny.length()
    assert length >= 0
    assert not np.isnan(length)
    assert not np.isinf(length)


def test_numerical_stability():
    """Test numerical operations don't produce NaN/Inf."""
    # Test with very small coordinates
    coords1 = np.array([[1e-12, 1e-12, 1e-12], [2e-12, 2e-12, 2e-12]])
    coords2 = np.array([[1.1e-12, 1.1e-12, 1.1e-12], [2.1e-12, 2.1e-12, 2.1e-12]])

    result = rmsd(coords1, coords2)
    assert not np.isnan(result)
    assert not np.isinf(result)
    assert result >= 0


if __name__ == "__main__":
    # Run tests if script is executed directly
    test_vector3d_basic_operations()
    test_atom_creation()
    test_atoms_collection()
    test_atoms_get_coordinates()
    test_rmsd_calculation()
    test_atoms_string_representation()
    test_atoms_sequence_creation()
    test_vector3d_edge_cases()
    test_numerical_stability()
    print("All tests passed!")
