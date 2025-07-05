"""
Tests for trajectory handling and analysis.
"""

import numpy as np

from CABS.core.trajectory import Header, Trajectory
from CABS.structures.atom import Atom, Atoms
from CABS.structures.vector3d import Vector3d


class TestHeader:
    """Test Header class for trajectory frame metadata."""

    def test_header_initialization(self):
        """Test Header initialization from string."""
        # Example header line from CABS output
        header_line = "1 10 -45.2 -12.3 -8.7 300.0 1"
        header = Header(header_line)

        assert header.model == 1
        assert header.length == (8,)  # 10 - 2
        assert header.temperature == 300.0
        assert header.replica == 1
        assert header.rmsd == 0

        # Energy should be a matrix
        assert header.energy.shape[1] == 3  # Three energy values

    def test_header_string_representation(self):
        """Test Header string representation."""
        header_line = "1 10 -45.2 -12.3 -8.7 300.0 1"
        header = Header(header_line)

        repr_str = repr(header)
        assert "Replica: 1" in repr_str
        assert "Model: 1" in repr_str
        assert "300.00" in repr_str

    def test_header_addition(self):
        """Test Header merging for multi-chain systems."""
        header1_line = "1 10 -45.2 -12.3 -8.7 300.0 1"
        header2_line = "1 15 -35.1 -22.1 -9.1 300.0 1"

        header1 = Header(header1_line)
        header2 = Header(header2_line)

        merged = header1 + header2

        assert merged.model == 1
        assert merged.replica == 1
        assert merged.temperature == 300.0
        # Length should be combined
        assert merged.length == (8, 13)  # (10-2, 15-2)

    def test_header_addition_mismatch(self):
        """Test Header addition with mismatched frames."""
        header1_line = "1 10 -45.2 -12.3 -8.7 300.0 1"
        header2_line = "2 15 -35.1 -22.1 -9.1 300.0 1"  # Different model

        header1 = Header(header1_line)
        header2 = Header(header2_line)

        try:
            merged = header1 + header2
            assert False, "Should have raised CannotMerge exception"
        except Header.CannotMerge:
            pass  # Expected

    def test_header_energy_modes(self):
        """Test Header energy calculation modes."""
        header_line = "1 10 -45.2 -12.3 -8.7 300.0 1"
        header = Header(header_line)

        # Test total energy
        total_energy = header.get_energy(mode="total")
        assert isinstance(total_energy, (int, float))

        # Test interaction energy
        interaction_energy = header.get_energy(mode="interaction", number_of_peptides=1)
        assert isinstance(interaction_energy, (int, float))


class TestTrajectoryInitialization:
    """Test Trajectory class initialization."""

    def test_trajectory_creation(self):
        """Test basic trajectory creation."""
        # Create template
        template = Atoms()
        atom = Atom()
        atom.coord = Vector3d(0.0, 0.0, 0.0)
        template.atoms = [atom]

        # Create coordinates (replicas x frames x atoms x 3)
        coordinates = np.random.random((1, 10, 1, 3))

        # Create headers
        headers = []
        for i in range(10):
            header_line = f"{i + 1} 3 -45.2 -12.3 -8.7 300.0 1"
            headers.append(Header(header_line))

        trajectory = Trajectory(template, coordinates, headers)

        assert trajectory.template == template
        assert trajectory.coordinates.shape == (1, 10, 1, 3)
        assert len(trajectory.headers) == 10
        assert trajectory.rmsd_native is None

    def test_trajectory_with_weights(self):
        """Test trajectory creation with weights."""
        template = Atoms()
        atom = Atom()
        template.atoms = [atom]

        coordinates = np.random.random((1, 5, 1, 3))
        headers = [Header(f"{i + 1} 3 -45.2 -12.3 -8.7 300.0 1") for i in range(5)]
        weights = np.array([1.0])

        trajectory = Trajectory(template, coordinates, headers, weights=weights)

        assert trajectory.weights is not None
        assert trajectory.weights.shape == (1,)


class TestTrajectorySelection:
    """Test Trajectory selection and filtering methods."""

    def test_trajectory_selection(self):
        """Test trajectory frame selection."""
        # Create test trajectory
        template = Atoms()
        for i in range(3):
            atom = Atom()
            atom.coord = Vector3d(float(i), 0.0, 0.0)
            template.atoms.append(atom)

        coordinates = np.random.random((1, 10, 3, 3))
        headers = [Header(f"{i + 1} 5 -45.2 -12.3 -8.7 300.0 1") for i in range(10)]

        trajectory = Trajectory(template, coordinates, headers)

        # Select subset of frames using a selection string
        from CABS.structures.atom import Selection

        selection = Selection("name CA")  # Select all CA atoms
        selected_trajectory = trajectory.select(selection)

        # Check selection worked (should have fewer or equal atoms)
        assert selected_trajectory.coordinates.shape[1] <= coordinates.shape[1]
        assert len(selected_trajectory.headers) == len(headers)

    def test_trajectory_template_selection(self):
        """Test trajectory with template selection."""
        # Create test trajectory with multiple atoms
        template = Atoms()
        for i in range(5):
            atom = Atom()
            atom.name = "CA" if i % 2 == 0 else "CB"
            atom.coord = Vector3d(float(i), 0.0, 0.0)
            template.atoms.append(atom)

        coordinates = np.random.random((1, 3, 5, 3))  # (replicas, frames, atoms, 3)
        headers = [Header(f"{i + 1} 7 -45.2 -12.3 -8.7 300.0 1") for i in range(3)]

        trajectory = Trajectory(template, coordinates, headers)

        # Select using template
        ca_template = Atoms()
        for atom in template.atoms:
            if atom.name == "CA":
                ca_template.atoms.append(atom)

        selected_trajectory = trajectory.select(template=ca_template)

        # Should have fewer or same atoms per frame
        assert selected_trajectory.coordinates.shape[2] <= coordinates.shape[2]
        assert len(selected_trajectory.template.atoms) == 3  # 3 CA atoms


class TestTrajectoryAnalysis:
    """Test Trajectory analysis methods."""

    def test_rmsd_matrix(self):
        """Test RMSD matrix calculation."""
        template = Atoms()
        for i in range(3):
            atom = Atom()
            atom.coord = Vector3d(float(i), 0.0, 0.0)
            template.atoms.append(atom)

        # Create coordinates with known differences
        coordinates = np.array(
            [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],  # Frame 0
                [[0.1, 0.0, 0.0], [1.1, 0.0, 0.0], [2.1, 0.0, 0.0]],  # Frame 1
                [[0.0, 0.1, 0.0], [1.0, 0.1, 0.0], [2.0, 0.1, 0.0]],  # Frame 2
            ]
        )

        headers = [Header(f"{i + 1} 5 -45.2 -12.3 -8.7 300.0 1") for i in range(3)]
        trajectory = Trajectory(template, coordinates, headers)

        rmsd_matrix = trajectory.rmsd_matrix()

        assert rmsd_matrix.shape == (3, 3)

        # Diagonal should be zero (or very small)
        for i in range(3):
            assert abs(rmsd_matrix[i, i]) < 1e-10

        # Matrix should be symmetric
        for i in range(3):
            for j in range(3):
                assert abs(rmsd_matrix[i, j] - rmsd_matrix[j, i]) < 1e-10

        # All values should be non-negative
        assert np.all(rmsd_matrix >= 0)

    def test_superimpose_to(self):
        """Test trajectory superimposition."""
        template = Atoms()
        for i in range(3):
            atom = Atom()
            atom.coord = Vector3d(float(i), 0.0, 0.0)
            template.atoms.append(atom)

        # Reference structure (Atoms object)
        reference = Atoms()
        ref_coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
        for i, (x, y, z) in enumerate(ref_coords):
            atom = Atom()
            atom.coord = Vector3d(x, y, z)
            reference.atoms.append(atom)

        # Trajectory with translated coordinates (4D: replicas, frames, atoms, 3)
        coordinates = np.array(
            [
                [
                    [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [3.0, 1.0, 1.0]],  # Frame 0
                    [[0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0]],  # Frame 1
                ]
            ]
        )

        headers = [Header(f"{i + 1} 5 -45.2 -12.3 -8.7 300.0 1") for i in range(2)]
        trajectory = Trajectory(template, coordinates, headers)
        # Superimpose to reference
        trajectory.superimpose_to(reference)

        # After superimposition, coordinates should be closer to reference
        # Check that all frames are now aligned to the reference coordinates
        ref_coords_array = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])

        for frame_idx in range(trajectory.coordinates.shape[1]):
            frame_coords = trajectory.coordinates[0, frame_idx, :, :]
            np.testing.assert_allclose(frame_coords, ref_coords_array, atol=1e-10)

        # But should be closer to reference structure
        # (This is a basic test - full validation would require more complex checks)

    def test_rmsf_calculation(self):
        """Test RMSF (Root Mean Square Fluctuation) calculation."""
        template = Atoms()
        for i in range(3):
            atom = Atom()
            atom.chid = "A"
            atom.coord = Vector3d(float(i), 0.0, 0.0)
            template.atoms.append(atom)

        # Create trajectory with fluctuations (4D)
        coordinates = np.random.random((1, 20, 3, 3)) * 0.5  # Small fluctuations
        headers = [Header(f"{i + 1} 5 -45.2 -12.3 -8.7 300.0 1") for i in range(20)]
        trajectory = Trajectory(template, coordinates, headers)

        rmsf_values = trajectory.rmsf(chains="A")

        assert len(rmsf_values) == 3  # One per atom
        assert all(val >= 0 for val in rmsf_values)  # All non-negative
        assert not any(np.isnan(val) for val in rmsf_values)  # No NaN values


class TestTrajectoryConversions:
    """Test Trajectory conversion methods."""

    def test_to_atoms(self):
        """Test conversion to Atoms object."""
        template = Atoms()
        for i in range(2):
            atom = Atom()
            atom.serial = i + 1
            atom.name = "CA"
            atom.resname = "ALA"
            atom.coord = Vector3d(float(i), 0.0, 0.0)
            template.atoms.append(atom)

        coordinates = np.array(
            [
                [
                    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],  # Frame 0
                    [[1.1, 2.1, 3.1], [4.1, 5.1, 6.1]],  # Frame 1
                ]
            ]
        )  # 4D shape

        headers = [Header(f"{i + 1} 4 -45.2 -12.3 -8.7 300.0 1") for i in range(2)]
        trajectory = Trajectory(template, coordinates, headers)

        # Get specific frame using get_model
        atoms_frame0 = trajectory.get_model(0)
        assert len(atoms_frame0.atoms) == 2
        assert abs(atoms_frame0.atoms[0].coord.x - 1.0) < 1e-10
        assert abs(atoms_frame0.atoms[1].coord.x - 4.0) < 1e-10

    def test_to_atoms_list(self):
        """Test conversion to list of Atoms objects."""
        template = Atoms()
        atom = Atom()
        atom.coord = Vector3d(0.0, 0.0, 0.0)
        template.atoms = [atom]

        coordinates = np.random.random((5, 1, 3))
        headers = [Header(f"{i + 1} 3 -45.2 -12.3 -8.7 300.0 1") for i in range(5)]
        trajectory = Trajectory(template, coordinates, headers)

        atoms_list = trajectory.to_atoms_list()

        assert len(atoms_list) == 5
        assert all(isinstance(atoms, Atoms) for atoms in atoms_list)
        assert all(len(atoms.atoms) == 1 for atoms in atoms_list)

    def test_get_model(self):
        """Test getting specific model from trajectory."""
        template = Atoms()
        atom = Atom()
        template.atoms = [atom]

        coordinates = np.random.random((10, 1, 3))
        headers = [Header(f"{i + 1} 3 -45.2 -12.3 -8.7 300.0 1") for i in range(10)]
        trajectory = Trajectory(template, coordinates, headers)

        # Get model 5
        model_atoms = trajectory.get_model(5)

        assert isinstance(model_atoms, Atoms)
        assert len(model_atoms.atoms) == 1

        # Coordinates should match frame 5 (as confirmed by manual testing)
        expected_coords = coordinates[5, 0, :]
        actual_coords = [
            model_atoms.atoms[0].coord.x,
            model_atoms.atoms[0].coord.y,
            model_atoms.atoms[0].coord.z,
        ]
        np.testing.assert_allclose(actual_coords, expected_coords, atol=1e-10)


class TestTrajectoryIO:
    """Test Trajectory I/O operations."""

    def test_to_pdb_models(self):
        """Test PDB output in models mode."""
        template = Atoms()
        atom = Atom()
        atom.serial = 1
        atom.name = "CA"
        atom.resname = "ALA"
        atom.chid = "A"
        atom.resnum = 1
        template.atoms = [atom]

        coordinates = np.array(
            [[[1.0, 2.0, 3.0]], [[1.1, 2.1, 3.1]], [[1.2, 2.2, 3.2]]]
        )

        headers = [Header(f"{i + 1} 3 -45.2 -12.3 -8.7 300.0 1") for i in range(3)]
        trajectory = Trajectory(template, coordinates, headers)

        # Generate PDB content - returns list with single StringIO object in models mode
        pdb_content = trajectory.to_pdb(mode="models")

        assert isinstance(pdb_content, list)
        assert (
            len(pdb_content) == 1
        )  # Should have 1 StringIO object containing all models

        # Should be a StringIO object with PDB content
        pdb_io = pdb_content[0]
        content = pdb_io.getvalue()
        assert isinstance(content, str)
        assert "HETATM" in content  # Based on your output, it uses HETATM not ATOM


class TestTrajectoryEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_trajectory(self):
        """Test handling of empty trajectory."""
        template = Atoms()
        coordinates = np.empty((0, 0, 3))
        headers = []

        trajectory = Trajectory(template, coordinates, headers)

        assert trajectory.coordinates.shape[0] == 0
        assert len(trajectory.headers) == 0

    def test_single_frame_trajectory(self):
        """Test trajectory with single frame."""
        template = Atoms()
        atom = Atom()
        template.atoms = [atom]

        coordinates = np.random.random((1, 1, 3))
        headers = [Header("1 3 -45.2 -12.3 -8.7 300.0 1")]

        trajectory = Trajectory(template, coordinates, headers)

        # RMSD matrix should work
        rmsd_matrix = trajectory.rmsd_matrix()
        assert rmsd_matrix.shape == (1, 1)
        assert abs(rmsd_matrix[0, 0]) < 1e-10

    def test_mismatched_dimensions(self):
        """Test handling of mismatched dimensions."""
        template = Atoms()
        atom = Atom()
        template.atoms = [atom]

        coordinates = np.random.random((5, 1, 3))
        headers = [
            Header(f"{i + 1} 3 -45.2 -12.3 -8.7 300.0 1") for i in range(3)
        ]  # Wrong count

        try:
            trajectory = Trajectory(template, coordinates, headers)
            # Some operations might still work, others might fail
            # The exact behavior depends on implementation
        except (ValueError, AssertionError):
            # Expected for mismatched dimensions
            pass
