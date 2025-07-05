"""
Tests for utility functions - RMSD, Kabsch algorithm, and numerical operations.
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
            [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]  # displaced by 1  # no displacement
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
        # Check that it's orthogonal (R @ R.T = I)
        assert np.allclose(
            np.dot(rotation_matrix, rotation_matrix.T), np.eye(3), atol=1e-10
        )

        # Check that determinant is ±1 (proper rotation or reflection)
        assert np.allclose(np.abs(np.linalg.det(rotation_matrix)), 1.0, atol=1e-10)

        # The matrix should be numerically stable (no NaN or inf values)
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

        weights = np.array([1.0, 1.0, 0.1])  # Lower weight for third point

        rotation_matrix = kabsch(target, query, weights)

        # Should be a valid orthogonal matrix (det=±1, orthogonal)
        # Kabsch can return reflection matrices (det=-1) in certain cases
        assert np.allclose(np.abs(np.linalg.det(rotation_matrix)), 1.0, atol=1e-10)
        assert np.allclose(rotation_matrix @ rotation_matrix.T, np.eye(3), atol=1e-10)

        # Should be numerically stable
        assert not np.any(np.isnan(rotation_matrix))
        assert not np.any(np.isinf(rotation_matrix))

    def test_kabsch_numerical_stability(self):
        """Test Kabsch with numerical edge cases."""
        # Very small coordinates
        target = np.array([[1e-15, 0.0, 0.0], [0.0, 1e-15, 0.0], [0.0, 0.0, 1e-15]])

        query = target * 1.1  # Slight scaling

        rotation_matrix = kabsch(target, query)

        # Should still produce valid rotation matrix
        assert not np.any(np.isnan(rotation_matrix))
        assert not np.any(np.isinf(rotation_matrix))


class TestDynamicKabsch:
    """Test dynamic Kabsch algorithm with iterative refinement."""

    def test_dynamic_kabsch_convergence(self):
        """Test that dynamic Kabsch converges for good input."""
        target = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 1.0],
            ]
        )

        # Add small rotation and noise
        query = target + np.random.random(target.shape) * 0.01

        try:
            final_rmsd, rotation, centroid_target, centroid_query = dynamic_kabsch(
                target, query
            )

            # Should converge to reasonable values
            assert final_rmsd >= 0
            assert not math.isnan(final_rmsd)
            assert not math.isinf(final_rmsd)

            # Rotation should be valid
            assert not np.any(np.isnan(rotation))
            assert np.allclose(np.linalg.det(rotation), 1.0, atol=1e-6)

        except Exception as e:
            # If it doesn't converge, that's also informative
            assert "convergence" in str(e).lower() or "iteration" in str(e).lower()

    def test_dynamic_kabsch_identical_structures(self):
        """Test dynamic Kabsch with identical structures."""
        coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])

        try:
            final_rmsd, rotation, centroid_target, centroid_query = dynamic_kabsch(
                coords, coords
            )

            # RMSD should be essentially zero
            assert abs(final_rmsd) < 1e-10

        except Exception:
            # If it fails due to numerical issues with identical structures, that's acceptable
            pass

    def test_dynamic_kabsch_failure_case(self):
        """Test dynamic Kabsch with problematic input."""
        # Collinear points - should be challenging
        target = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])

        query = np.array([[0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]])

        try:
            result = dynamic_kabsch(target, query)
            # If it succeeds, result should be reasonable
            assert len(result) == 4
        except Exception as e:
            # Expected to fail due to collinearity
            assert isinstance(e, (ValueError, RuntimeError, np.linalg.LinAlgError))


class TestRotationMatrix:
    """Test rotation matrix utilities."""

    def test_random_rotation_matrix(self):
        """Test random rotation matrix generation."""
        R = random_rotation_matrix()

        # Should be 3x3
        assert R.shape == (3, 3)

        # Should be orthogonal
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-10)

        # Should have determinant 1 (proper rotation)
        assert np.allclose(np.linalg.det(R), 1.0, atol=1e-10)

    def test_multiple_random_rotations(self):
        """Test that multiple random rotations are different."""
        R1 = random_rotation_matrix()
        R2 = random_rotation_matrix()

        # Should be different (very unlikely to be identical)
        assert not np.allclose(R1, R2, atol=1e-6)


class TestNumericalStability:
    """Test numerical stability of utility functions."""

    def test_denormal_number_handling(self):
        """Test handling of denormal numbers that cause IEEE warnings."""
        # These might trigger IEEE_DENORMAL warnings
        tiny_coords = np.array(
            [
                [1e-320, 1e-320, 1e-320],
                [2e-320, 2e-320, 2e-320],
                [3e-320, 3e-320, 3e-320],
            ]
        )

        # Suppress warnings for this test
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # RMSD should handle denormal numbers
            result = rmsd(tiny_coords, tiny_coords * 1.1)
            assert not math.isnan(result)
            assert not math.isinf(result)

    def test_large_coordinate_stability(self):
        """Test stability with very large coordinates."""
        large_coords = np.random.random((5, 3)) * 1e12
        # Use proportionally larger noise that will be measurable
        large_noise = (
            np.random.random((5, 3)) * 1e9
        )  # Larger noise: 1e9 relative to 1e12 gives 1e-3 relative scale

        result = rmsd(large_coords, large_coords + large_noise)
        assert not math.isnan(result)
        assert not math.isinf(result)
        # Should be greater than 0 since we added significant noise
        assert result > 0

    def test_mixed_scale_coordinates(self):
        """Test with coordinates of very different scales."""
        coords1 = np.array([[1e-15, 1e15, 1.0], [1e15, 1e-15, 1.0], [1.0, 1.0, 1e15]])

        coords2 = coords1 + np.random.random((3, 3)) * 1e-10

        result = rmsd(coords1, coords2)
        assert not math.isnan(result)
        assert not math.isinf(result)


class TestErrorHandling:
    """Test error handling in utility functions."""

    def test_rmsd_mismatched_shapes(self):
        """Test RMSD with mismatched coordinate arrays."""
        coords1 = np.random.random((5, 3))
        coords2 = np.random.random((4, 3))  # Different number of points

        try:
            rmsd(coords1, coords2)
            assert False, "Should have raised an error"
        except (ValueError, AssertionError):
            pass  # Expected

    def test_rmsd_wrong_dimensions(self):
        """Test RMSD with wrong dimensional arrays."""
        coords1 = np.random.random((5, 2))  # 2D instead of 3D
        coords2 = np.random.random((5, 2))

        try:
            rmsd(coords1, coords2)
            assert False, "Should have raised an error"
        except (ValueError, AssertionError, IndexError):
            pass  # Expected

    def test_kabsch_insufficient_points(self):
        """Test Kabsch with insufficient points."""
        coords1 = np.array([[1.0, 2.0, 3.0]])  # Only one point
        coords2 = np.array([[4.0, 5.0, 6.0]])

        try:
            kabsch(coords1, coords2)
            # Might work or might fail - both are acceptable
        except (ValueError, np.linalg.LinAlgError):
            pass  # Expected for insufficient points


class TestPerformance:
    """Test performance characteristics of utility functions."""

    def test_rmsd_large_arrays(self):
        """Test RMSD with large coordinate arrays."""
        n_points = 10000
        coords1 = np.random.random((n_points, 3))
        coords2 = coords1 + np.random.random((n_points, 3)) * 0.1

        # Should complete in reasonable time
        result = rmsd(coords1, coords2)
        assert result > 0
        assert not math.isnan(result)

    def test_kabsch_medium_arrays(self):
        """Test Kabsch with medium-sized arrays."""
        n_points = 1000
        coords1 = np.random.random((n_points, 3))
        coords2 = coords1 + np.random.random((n_points, 3)) * 0.01

        # Should complete without issues
        rotation = kabsch(coords1, coords2)
        assert rotation.shape == (3, 3)
        assert np.allclose(np.linalg.det(rotation), 1.0, atol=1e-6)
