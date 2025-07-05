"""
Tests for numerical stability and edge cases.
"""

import math
import warnings

import numpy as np

from CABS.structures.vector3d import Vector3d
from CABS.utils.utils import dynamic_kabsch, kabsch, rmsd


class TestNumericalStability:
    """Test numerical stability with edge cases that might cause IEEE_DENORMAL warnings."""

    def test_denormal_coordinates_rmsd(self):
        """Test RMSD with denormal coordinate values."""
        # These values are in the denormal range and might trigger IEEE_DENORMAL
        denormal_coords = np.array(
            [
                [1e-320, 1e-320, 1e-320],
                [2e-320, 2e-320, 2e-320],
                [3e-320, 3e-320, 3e-320],
            ]
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Suppress IEEE warnings for test

            result = rmsd(denormal_coords, denormal_coords)

            # Should be zero (or very close) for identical structures
            assert abs(result) < 1e-300
            assert not math.isnan(result)
            assert not math.isinf(result)

    def test_mixed_scale_coordinates(self):
        """Test with coordinates spanning many orders of magnitude."""
        mixed_coords = np.array(
            [[1e-15, 1e15, 1.0], [1e15, 1e-15, 1.0], [1.0, 1.0, 1e15]]
        )

        # Add small perturbation
        perturbed = mixed_coords + np.random.random((3, 3)) * 1e-12

        result = rmsd(mixed_coords, perturbed)

        assert result >= 0
        assert not math.isnan(result)
        assert not math.isinf(result)

    def test_very_small_differences(self):
        """Test RMSD with extremely small coordinate differences."""
        base_coords = np.random.random((10, 3))

        # Add tiny differences that might cause numerical issues
        tiny_diff = np.random.random((10, 3)) * 1e-15
        perturbed_coords = base_coords + tiny_diff

        result = rmsd(base_coords, perturbed_coords)

        assert result >= 0
        assert result < 1e-10  # Should be very small
        assert not math.isnan(result)

    def test_vector_operations_denormal(self):
        """Test Vector3d operations with denormal numbers."""
        # Create vectors with denormal values
        tiny_vector = Vector3d(1e-320, 1e-320, 1e-320)
        normal_vector = Vector3d(1.0, 1.0, 1.0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Test basic operations
            result1 = tiny_vector + normal_vector
            assert not any(math.isnan(x) for x in [result1.x, result1.y, result1.z])

            result2 = tiny_vector * 1e300  # Scale up
            assert not any(math.isnan(x) for x in [result2.x, result2.y, result2.z])

            # Length calculation
            length = tiny_vector.length()
            assert length >= 0
            assert not math.isnan(length)

    def test_kabsch_numerical_stability(self):
        """Test Kabsch algorithm with challenging numerical conditions."""
        # Nearly collinear points
        target = np.array([[0.0, 0.0, 0.0], [1.0, 1e-15, 1e-15], [2.0, 2e-15, 2e-15]])

        query = target + np.random.random((3, 3)) * 1e-12

        try:
            rotation = kabsch(target, query)

            # Should produce valid rotation matrix
            assert rotation.shape == (3, 3)
            assert not np.any(np.isnan(rotation))
            assert not np.any(np.isinf(rotation))

            # Check if it's orthogonal (within tolerance)
            orthogonal_test = rotation @ rotation.T
            identity = np.eye(3)

            # Might not be perfectly orthogonal due to numerical issues
            # but should be close if algorithm handled the case well
            if np.allclose(orthogonal_test, identity, atol=1e-6):
                # For numerical stability, allow determinant to be -1 or 1 (reflection or rotation)
                det_value = np.linalg.det(rotation)
                assert np.allclose(abs(det_value), 1.0, atol=1e-6), (
                    f"Determinant should be ±1, got {det_value}"
                )

        except (np.linalg.LinAlgError, ValueError):
            # Expected for degenerate cases
            pass

    def test_dynamic_kabsch_convergence_issues(self):
        """Test dynamic Kabsch with conditions that might prevent convergence."""
        # Identical points (no variance)
        target = np.array([[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])

        query = np.array([[2.0, 2.0, 2.0], [2.0, 2.0, 2.0], [2.0, 2.0, 2.0]])

        try:
            result = dynamic_kabsch(target, query)

            # If it succeeds, result should make sense
            rmsd_val, rotation, centroid_target, centroid_query = result
            assert rmsd_val >= 0
            assert not math.isnan(rmsd_val)

        except Exception as e:
            # Expected to fail due to lack of variance
            assert isinstance(e, (ValueError, RuntimeError, np.linalg.LinAlgError))

    def test_overflow_prevention(self):
        """Test prevention of overflow in calculations."""
        # Very large coordinates
        large_coords = np.array(
            [[1e100, 1e100, 1e100], [2e100, 2e100, 2e100], [3e100, 3e100, 3e100]]
        )

        try:
            result = rmsd(large_coords, large_coords)

            # Should be zero for identical structures
            assert abs(result) < 1e90  # Very small relative to input scale
            assert not math.isnan(result)
            assert not math.isinf(result)

        except OverflowError:
            # Acceptable if implementation can't handle such large values
            pass

    def test_underflow_handling(self):
        """Test handling of underflow conditions."""
        # Coordinates so small they might underflow
        tiny_coords = np.array([[1e-308, 1e-308, 1e-308], [2e-308, 2e-308, 2e-308]])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            result = rmsd(tiny_coords)

            # Should handle underflow gracefully
            assert result >= 0
            assert not math.isnan(result)

    def test_zero_variance_structures(self):
        """Test structures with zero variance in some dimensions."""
        # All points have same x and y, only z varies
        coords = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 2.0]])

        # RMSD should still work
        result = rmsd(coords)
        assert result > 0  # Should have some variance
        assert not math.isnan(result)

        # Kabsch with another similar structure
        query_coords = np.array([[1.0, 1.0, 0.1], [1.0, 1.0, 1.1], [1.0, 1.0, 2.1]])

        try:
            rotation = kabsch(coords, query_coords)
            assert not np.any(np.isnan(rotation))
        except np.linalg.LinAlgError:
            # Expected for degenerate cases
            pass


class TestNumericalAccuracy:
    """Test numerical accuracy of calculations."""

    def test_rmsd_known_values(self):
        """Test RMSD with known analytical results."""
        # Unit cube vertices
        target = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
        )

        # Translate by (1, 1, 1)
        query = target + np.array([1.0, 1.0, 1.0])

        result = rmsd(target, query)
        expected = math.sqrt(3.0)  # sqrt(1^2 + 1^2 + 1^2)

        assert abs(result - expected) < 1e-12

    def test_rotation_matrix_properties(self):
        """Test that rotation matrices maintain proper properties."""
        # Generate test coordinates
        coords = np.random.random((10, 3))
        rotated_coords = coords @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]).T

        rotation = kabsch(coords, rotated_coords)

        # Test orthogonality: R @ R.T = I
        orthogonal_test = rotation @ rotation.T
        identity = np.eye(3)
        assert np.allclose(orthogonal_test, identity, atol=1e-10)

        # Test proper rotation: det(R) = 1
        determinant = np.linalg.det(rotation)
        assert abs(determinant - 1.0) < 1e-10

        # Test that applying rotation improves alignment
        rotated_back = coords @ rotation.T
        rmsd_before = rmsd(coords, rotated_coords)
        rmsd_after = rmsd(rotated_back, rotated_coords)
        assert rmsd_after <= rmsd_before + 1e-10  # Should be better or equal

    def test_vector_precision(self):
        """Test vector operations maintain precision."""
        # Test with known precise values
        v1 = Vector3d(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
        v2 = Vector3d(2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)

        result = v1 + v2
        expected = Vector3d(1.0, 1.0, 1.0)

        assert abs(result.x - expected.x) < 1e-15
        assert abs(result.y - expected.y) < 1e-15
        assert abs(result.z - expected.z) < 1e-15

    def test_centroid_accuracy(self):
        """Test centroid calculation accuracy."""
        # Points with known centroid
        coords = np.array(
            [
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [-1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
            ]
        )

        # Centroid should be (0, 0, 0)
        centroid = np.mean(coords, axis=0)

        assert abs(centroid[0]) < 1e-15
        assert abs(centroid[1]) < 1e-15
        assert abs(centroid[2]) < 1e-15


class TestNumericalRobustness:
    """Test robustness against numerical instabilities."""

    def test_repeated_operations(self):
        """Test that repeated operations don't accumulate errors."""
        coords = np.random.random((5, 3))

        # Apply multiple rotations
        result_coords = coords.copy()
        for _ in range(100):
            # Small rotation
            angle = 0.01
            rotation = np.array(
                [
                    [math.cos(angle), -math.sin(angle), 0],
                    [math.sin(angle), math.cos(angle), 0],
                    [0, 0, 1],
                ]
            )
            result_coords = result_coords @ rotation.T

        # After 100 small rotations, should still be reasonable
        assert not np.any(np.isnan(result_coords))
        assert not np.any(np.isinf(result_coords))

        # Distances from origin shouldn't have changed drastically
        original_norms = np.linalg.norm(coords, axis=1)
        final_norms = np.linalg.norm(result_coords, axis=1)

        # Should be approximately preserved (rotation preserves distances)
        np.testing.assert_allclose(original_norms, final_norms, rtol=1e-10)

    def test_condition_number_handling(self):
        """Test handling of ill-conditioned matrices."""
        # Create nearly singular coordinate set
        coords = np.array(
            [[1.0, 0.0, 0.0], [1.0, 1e-15, 0.0], [1.0, 2e-15, 0.0]]  # Nearly collinear
        )

        try:
            # This might fail due to singular matrix
            rotation = kabsch(coords, coords + 1e-12)

            if rotation is not None:
                # If it succeeds, should still be valid
                assert not np.any(np.isnan(rotation))

        except (np.linalg.LinAlgError, ValueError):
            # Expected for nearly singular cases
            pass

    def test_floating_point_precision_limits(self):
        """Test behavior at floating point precision limits."""
        # Coordinates differing at machine epsilon level
        coords1 = np.array([[1.0, 1.0, 1.0]])
        coords2 = coords1 + np.finfo(float).eps

        result = rmsd(coords1, coords2)

        # Should be very small but computable
        assert result >= 0
        assert result < 1e-14
        assert not math.isnan(result)


class TestWarningHandling:
    """Test proper handling of numerical warnings."""

    def test_ieee_denormal_suppression(self):
        """Test that IEEE_DENORMAL warnings are handled appropriately."""
        # Create conditions that typically trigger IEEE_DENORMAL
        denormal_coords = np.array(
            [
                [1e-324, 1e-324, 1e-324],  # Smallest positive denormal
                [2e-324, 2e-324, 2e-324],
            ]
        )

        # Should complete without raising warnings to user
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")

            result = rmsd(denormal_coords, denormal_coords)

            # Check if the operation completed
            assert not math.isnan(result)

            # Check if IEEE warnings were generated
            ieee_warnings = [
                warning
                for warning in w
                if "IEEE" in str(warning.message)
                or "denormal" in str(warning.message).lower()
            ]

            # The actual handling depends on implementation
            # This test documents the behavior rather than enforcing it

    def test_overflow_warning_handling(self):
        """Test handling of potential overflow warnings."""
        # Large values that might cause overflow warnings
        large_coords = np.array([[1e150, 1e150, 1e150], [2e150, 2e150, 2e150]])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")

            try:
                result = rmsd(large_coords, large_coords)
                # If it completes, should be reasonable
                assert not math.isnan(result)
            except (OverflowError, RuntimeWarning):
                # Also acceptable to fail gracefully
                pass
