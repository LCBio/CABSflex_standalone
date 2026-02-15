"""
Tests for contact map analysis module.
"""

import sys
from unittest.mock import MagicMock

import types

# Create a mock package for matplotlib
mock_mpl = MagicMock()
sys.modules["matplotlib"] = mock_mpl

# Mock matplotlib.pyplot
sys.modules["matplotlib.pyplot"] = MagicMock()

# Mock matplotlib.ticker
sys.modules["matplotlib.ticker"] = MagicMock()

# Mock matplotlib.axes (Crucial fix)
mock_axes = MagicMock()
sys.modules["matplotlib.axes"] = mock_axes
# Also ensure it can be imported as from matplotlib.axes import Axes
mock_mpl.axes = mock_axes

# Mock matplotlib.colors
sys.modules["matplotlib.colors"] = MagicMock()

import numpy as np
import pytest

from CABS.analysis.cmap import ContactMapFactory, ContactMap
from CABS.structures.atom import Atom, Atoms
from CABS.structures.vector3d import Vector3d


class TestContactMapFactory:
    """Test ContactMapFactory class."""
    
    def setup_method(self):
        """Setup test fixtures."""
        self.atoms = Atoms()
        # Create Chain A (residues 1-2)
        for i in range(2):
            atom = Atom()
            atom.chid = "A"
            atom.resnum = i + 1
            atom.name = "CA"
            atom.coord = Vector3d(float(i * 10.0), 0.0, 0.0)
            self.atoms.append(atom)
            
        # Create Chain B (residues 1-2)
        for i in range(2):
            atom = Atom()
            atom.chid = "B"
            atom.resnum = i + 1
            atom.name = "CA"
            atom.coord = Vector3d(float(i * 10.0) + 3.0, 0.0, 0.0) # Dist 3.0 from A
            self.atoms.append(atom)

        # Initialize Factory
        self.factory = ContactMapFactory(chains1="A", chains2="B", temp=self.atoms)

    def test_dimensions(self):
        """Test correct dimension calculation."""
        # 2 atoms in A, 2 atoms in B -> (2, 2)
        assert self.factory.dims == (2, 2)
        
    def test_mk_dmtx(self):
        """Test distance matrix calculation logic."""
        # Create dummy coordinates matching template order: A1, A2, B1, B2
        coords = np.array([
            [0.0, 0.0, 0.0],    # A1
            [10.0, 0.0, 0.0],   # A2
            [0.0, 3.0, 0.0],    # B1 (Dist 3.0 to A1)
            [10.0, 4.0, 0.0]    # B2 (Dist 4.0 to A2)
        ])
        
        # Calculate distance matrix using factory method
        dmtx = self.factory.mk_dmtx(coords)
        
        assert dmtx.shape == (2, 2)
        assert abs(dmtx[0, 0] - 3.0) < 1e-5 # A1-B1
        assert abs(dmtx[1, 1] - 4.0) < 1e-5 # A2-B2

    def test_mk_cmtx(self):
        """Test contact boolean matrix creation."""
        # Distances:
        # A1-B1: 3.0 (Contact < 3.5)
        # A1-B2: 10.0 (No contact)
        # A2-B1: 10.0 (No contact)
        # A2-B2: 4.0 (No contact > 3.5)
        
        dmtx = np.array([[3.0, 10.0], [10.0, 4.0]])
        cmtx = self.factory.mk_cmtx(dmtx, thr=3.5)
        
        # Check logic
        assert cmtx[0, 0] == 1.0
        assert cmtx[1, 1] == 0.0


class TestContactMap:
    """Test ContactMap object behavior."""
    
    def test_addition(self):
        """Test adding contact maps."""
        # Map 1: 10 frames
        cm1 = ContactMap(np.array([[1, 0], [0, 0]]), ["A1", "A2"], ["B1", "B2"], n=10)
        # Map 2: 5 frames
        cm2 = ContactMap(np.array([[0, 1], [0, 0]]), ["A1", "A2"], ["B1", "B2"], n=5)
        
        cm_sum = cm1 + cm2
        
        assert cm_sum.n == 15
        np.testing.assert_array_equal(cm_sum.cmtx, np.array([[1, 1], [0, 0]]))

    def test_zero_diagonal(self):
        """Test zeroing diagonal."""
        cm = ContactMap(np.eye(2), ["A", "B"], ["A", "B"], 10)
        cm.zero_diagonal()
        np.testing.assert_array_equal(cm.cmtx, np.zeros((2, 2)))
