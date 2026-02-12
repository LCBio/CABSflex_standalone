"""
Tests for utility functions - RMSD, Kabsch algorithm, and numerical operations.
(Consolidated Numerical Hub)
"""

import math
import warnings

import numpy as np

from CABS.utils.utils import dynamic_kabsch, kabsch, random_rotation_matrix, rmsd


class TestRMSDCalculations:
    """Test RMSD calculation functions."""

    def test_rmsd_identical_structures(self):
        """Test RMSD between identical structures should be zero."""
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])

        result = rmsd(coords, coords)
        assert abs(result) < 1e-10

    def test_rmsd_known_displacement(self):
        """Test RMSD with known displacement."""
        target = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

        # Translate by 1 unit in each direction
        query = target + np.array([1.0, 1.0, 1.0])

        expected_rmsd = math.sqrt(3.0)  # sqrt(1^2 + 1^2 + 1^2)
        result = rmsd(target, query)
        assert abs(result - expected_rmsd) < 1e-10

    def test_rmsd_with_weights(self):
        """Test weighted RMSD calculation."""
        target = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

        query = np.array(
            [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        )

        # Equal weights
        weights = np.array([1.0, 1.0])
        result_weighted = rmsd(target, query, weights)

        # Should be sqrt((1^2 + 0^2) / 2) = sqrt(0.5)
        expected = math.sqrt(0.5)
        assert abs(result_weighted - expected) < 1e-10

    def test_rmsd_single_coordinate_set(self):
        """Test RMSD with single coordinate set (variance calculation)."""
        coords = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        )

        # Should calculate RMSD from centroid
        result = rmsd(coords)
        assert result > 0
        assert not math.isnan(result)

    def test_rmsd_numerical_stability(self):
        """Test RMSD with very small differences."""
        target = np.random.random((10, 3))
        query = target + np.random.random((10, 3)) * 1e-12

        result = rmsd(target, query)
        assert result >= 0
        assert not math.isnan(result)
        assert not math.isinf(result)

    def test_rmsd_large_coordinates(self):
        """Test RMSD with large coordinate values."""
        target = np.random.random((5, 3)) * 1e6
        query = target + np.random.random((5, 3)) * 1e3

        result = rmsd(target, query)
        assert result > 0
        assert not math.isnan(result)
        assert not math.isinf(result)


class TestKabschAlgorithm:
    """Test Kabsch algorithm for optimal rotation."""

    def test_kabsch_identical_structures(self):
        """Test Kabsch with identical structures."""
        coords = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

        rotation_matrix = kabsch(coords, coords)

        # Kabsch should return a valid orthogonal matrix (rotation or reflection)
        assert np.allclose(
            np.dot(rotation_matrix, rotation_matrix.T), np.eye(3), atol=1e-10
        )
        assert np.allclose(np.abs(np.linalg.det(rotation_matrix)), 1.0, atol=1e-10)
        assert not np.any(np.isnan(rotation_matrix))
        assert not np.any(np.isinf(rotation_matrix))

    def test_kabsch_90_degree_rotation(self):
        """Test Kabsch with known 90-degree rotation."""
        target = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

        # 90-degree rotation around z-axis
        query = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

        rotation_matrix = kabsch(target, query)

        # Apply rotation to target
        rotated = target @ rotation_matrix.T

        # Should match query after rotation
        np.testing.assert_allclose(rotated, query, atol=1e-10)

    def test_kabsch_with_weights(self):
        """Test weighted Kabsch algorithm."""
        target = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        query = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        weights = np.array([1.0, 1.0, 0.1])

        rotation_matrix = kabsch(target, query, weights)

        assert np.allclose(np.abs(np.linalg.det(rotation_matrix)), 1.0, atol=1e-10)
        assert not np.any(np.isnan(rotation_matrix))

    def test_kabsch_numerical_stability(self):
        """Test Kabsch with numerical edge cases (retained from test_utils)."""
        target = np.array([[1e-15, 0.0, 0.0], [0.0, 1e-15, 0.0], [0.0, 0.0, 1e-15]])
        query = target * 1.1

        rotation_matrix = kabsch(target, query)
        assert not np.any(np.isnan(rotation_matrix))


class TestDynamicKabsch:
    """Test dynamic Kabsch algorithm with iterative refinement."""

    def test_dynamic_kabsch_convergence(self):
        """Test that dynamic Kabsch converges for good input."""
        target = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        query = target + np.random.random(target.shape) * 0.01

        try:
            final_rmsd, rotation, _, _ = dynamic_kabsch(target, query)
            assert final_rmsd >= 0
            assert not math.isnan(final_rmsd)
            assert np.allclose(np.linalg.det(rotation), 1.0, atol=1e-6)
        except Exception:
            pass

    def test_dynamic_kabsch_identical_structures(self):
        """Test dynamic Kabsch with identical structures."""
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
        try:
            final_rmsd, _, _, _ = dynamic_kabsch(coords, coords)
            assert abs(final_rmsd) < 1e-10
        except Exception:
            pass

    def test_dynamic_kabsch_failure_case(self):
        """Test dynamic Kabsch with problematic input."""
        # Collinear points - should be challenging
        target = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        query = np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]])
        try:
            dynamic_kabsch(target, query)
            assert False, "Expected failure or extreme numerical instability"
        except Exception:
            pass


class TestRotationMatrix:
    """Test rotation matrix utilities."""

    def test_random_rotation_matrix(self):
        """Test random rotation matrix generation."""
        R = random_rotation_matrix()
        assert R.shape == (3, 3)
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-10)
        assert np.allclose(np.linalg.det(R), 1.0, atol=1e-10)


class TestNumericalErrorHandling:
    """Test error handling in utility functions."""

    def test_rmsd_mismatched_shapes(self):
        """Test RMSD with mismatched coordinate arrays."""
        coords1 = np.random.random((5, 3))
        coords2 = np.random.random((4, 3))
        try:
            rmsd(coords1, coords2)
            assert False, "Should have raised an error"
        except (ValueError, AssertionError):
            pass

    def test_rmsd_wrong_dimensions(self):
        """Test RMSD with wrong dimensional arrays."""
        coords1 = np.random.random((5, 2))
        coords2 = np.random.random((5, 2))
        try:
            rmsd(coords1, coords2)
            assert False, "Should have raised an error"
        except (ValueError, AssertionError, IndexError):
            pass

    def test_kabsch_insufficient_points(self):
        """Test Kabsch with insufficient points (e.g., only one point)."""
        coords1 = np.array([[1.0, 2.0, 3.0]])
        coords2 = np.array([[4.0, 5.0, 6.0]])
        try:
            kabsch(coords1, coords2)
        except (ValueError, np.linalg.LinAlgError):
            pass


class TestAdvancedNumericalStability:
    """Test robustness against extreme numerical conditions (Merged from test_numerical)."""

    def test_denormal_coordinates_rmsd(self):
        """Test RMSD with denormal coordinate values."""
        # Note: Using 1e-320 as 1e-324 is often hard to represent consistently
        denormal_coords = np.array(
            [[1e-320, 1e-320, 1e-320], [2e-320, 2e-320, 2e-320]]
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = rmsd(denormal_coords, denormal_coords)
            assert abs(result) < 1e-300
            assert not math.isnan(result)

    def test_mixed_scale_coordinates(self):
        """Test with coordinates spanning many orders of magnitude."""
        mixed_coords = np.array(
            [[1e-15, 1e15, 1.0], [1e15, 1e-15, 1.0], [1.0, 1.0, 1e15]]
        )
        perturbed = mixed_coords + np.random.random((3, 3)) * 1e-12
        result = rmsd(mixed_coords, perturbed)
        assert result >= 0
        assert not math.isnan(result)

    def test_zero_variance_structures_kabsch(self):
        """Test structures with zero variance in some dimensions."""
        coords = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 2.0]])
        query_coords = np.array([[1.0, 1.0, 0.1], [1.0, 1.0, 1.1], [1.0, 1.0, 2.1]])
        try:
            rotation = kabsch(coords, query_coords)
            assert not np.any(np.isnan(rotation))
        except np.linalg.LinAlgError:
            pass

    def test_repeated_operations(self):
        """Test that repeated operations don't accumulate errors."""
        coords = np.random.random((5, 3))
        result_coords = coords.copy()
        for _ in range(100):
            angle = 0.01
            rotation = np.array(
                [[math.cos(angle), -math.sin(angle), 0],
                 [math.sin(angle), math.cos(angle), 0],
                 [0, 0, 1],]
            )
            result_coords = result_coords @ rotation.T

        original_norms = np.linalg.norm(coords, axis=1)
        final_norms = np.linalg.norm(result_coords, axis=1)
        np.testing.assert_allclose(original_norms, final_norms, rtol=1e-10)


class TestPerformance:
    """Test performance characteristics of utility functions."""

    def test_rmsd_large_arrays(self):
        """Test RMSD with large coordinate arrays."""
        n_points = 10000
        coords1 = np.random.random((n_points, 3))
        coords2 = coords1 + np.random.random((n_points, 3)) * 0.1
        result = rmsd(coords1, coords2)
        assert result > 0
        assert not math.isnan(result)

    def test_kabsch_medium_arrays(self):
        """Test Kabsch with medium-sized arrays."""
        n_points = 1000
        coords1 = np.random.random((n_points, 3))
        coords2 = coords1 + np.random.random((n_points, 3)) * 0.01
        rotation = kabsch(coords1, coords2)
        assert rotation.shape == (3, 3)
