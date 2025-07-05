#!/usr/bin/env python3
"""
Vector3d tests for CABS-flex.
Tests vector operations and numerical stability.
"""

import numpy as np

from CABS.structures.vector3d import Vector3d


class TestVector3dBasics:
    """Test basic Vector3d functionality."""

    def test_vector_creation(self):
        """Test vector creation."""
        # Default constructor
        v1 = Vector3d()
        assert v1.x == 0.0
        assert v1.y == 0.0
        assert v1.z == 0.0

        # Constructor with values
        v2 = Vector3d(1.0, 2.0, 3.0)
        assert v2.x == 1.0
        assert v2.y == 2.0
        assert v2.z == 3.0

        # Constructor from array-like
        v3 = Vector3d([4.0, 5.0, 6.0])
        assert v3.x == 4.0
        assert v3.y == 5.0
        assert v3.z == 6.0

    def test_vector_arithmetic(self):
        """Test vector arithmetic operations."""
        v1 = Vector3d(1.0, 2.0, 3.0)
        v2 = Vector3d(4.0, 5.0, 6.0)

        # Addition
        v3 = v1 + v2
        assert v3.x == 5.0
        assert v3.y == 7.0
        assert v3.z == 9.0

        # Subtraction
        v4 = v2 - v1
        assert v4.x == 3.0
        assert v4.y == 3.0
        assert v4.z == 3.0

        # Scalar multiplication
        v5 = v1 * 2.0
        assert v5.x == 2.0
        assert v5.y == 4.0
        assert v5.z == 6.0

        # Scalar division
        v6 = v1 / 2.0
        assert v6.x == 0.5
        assert v6.y == 1.0
        assert v6.z == 1.5

    def test_vector_properties(self):
        """Test vector properties."""
        v = Vector3d(3.0, 4.0, 0.0)

        # Length calculation
        length = v.length()
        assert abs(length - 5.0) < 1e-10

        # Squared length
        length_sq = v.mod2()
        assert abs(length_sq - 25.0) < 1e-10

    def test_vector_normalization(self):
        """Test vector normalization."""
        v = Vector3d(3.0, 4.0, 0.0)
        normalized = v.norm()

        # Should be unit vector
        assert abs(normalized.length() - 1.0) < 1e-10

        # Direction should be preserved
        assert abs(normalized.x - 0.6) < 1e-10
        assert abs(normalized.y - 0.8) < 1e-10
        assert abs(normalized.z - 0.0) < 1e-10

    def test_vector_dot_product(self):
        """Test dot product."""
        v1 = Vector3d(1.0, 0.0, 0.0)
        v2 = Vector3d(0.0, 1.0, 0.0)
        v3 = Vector3d(1.0, 1.0, 0.0)

        # Perpendicular vectors
        dot1 = v1.dot(v2)
        assert abs(dot1) < 1e-10

        # Parallel vectors
        dot2 = v1.dot(v1)
        assert abs(dot2 - 1.0) < 1e-10

        # 45-degree vectors
        dot3 = v1.dot(v3)
        assert abs(dot3 - 1.0) < 1e-10

    def test_vector_cross_product(self):
        """Test cross product."""
        v1 = Vector3d(1.0, 0.0, 0.0)
        v2 = Vector3d(0.0, 1.0, 0.0)

        cross = v1.cross(v2)
        expected = Vector3d(0.0, 0.0, 1.0)

        assert abs(cross.x - expected.x) < 1e-10
        assert abs(cross.y - expected.y) < 1e-10
        assert abs(cross.z - expected.z) < 1e-10

    def test_vector_equality(self):
        """Test vector equality comparison."""
        v1 = Vector3d(1.0, 2.0, 3.0)
        v2 = Vector3d(1.0, 2.0, 3.0)
        v3 = Vector3d(1.1, 2.0, 3.0)

        assert v1 == v2
        assert v1 != v3

    def test_vector_string_representation(self):
        """Test vector string representation."""
        v = Vector3d(1.234, 2.567, 3.890)
        s = str(v)
        assert isinstance(s, str)
        assert "1.234" in s or "1.23" in s  # Depending on formatting
        assert "2.567" in s or "2.57" in s
        assert "3.890" in s or "3.89" in s


class TestVector3dNumerical:
    """Test numerical behavior of Vector3d."""

    def test_tiny_vectors(self):
        """Test operations with very small vectors."""
        tiny = Vector3d(1e-15, 1e-15, 1e-15)

        # Should handle tiny numbers
        length = tiny.length()
        assert length > 0
        assert not np.isnan(length)
        assert not np.isinf(length)

        # Squared length should be even smaller
        length_sq = tiny.mod2()
        assert length_sq > 0
        assert length_sq < length

    def test_large_vectors(self):
        """Test operations with very large vectors."""
        large = Vector3d(1e12, 1e12, 1e12)

        length = large.length()
        assert not np.isnan(length)
        assert not np.isinf(length)
        assert length > 1e12

    def test_normalization_edge_cases(self):
        """Test normalization edge cases."""
        # Very small but non-zero vector
        tiny = Vector3d(1e-14, 1e-14, 1e-14)
        try:
            normalized = tiny.norm()
            # If it succeeds, should be unit vector
            assert abs(normalized.length() - 1.0) < 1e-6
        except ZeroDivisionError:
            # Acceptable for extremely small vectors
            pass

        # Zero vector should raise error
        zero = Vector3d(0.0, 0.0, 0.0)
        try:
            zero.norm()
            assert False, "Zero vector normalization should raise exception"
        except ZeroDivisionError:
            pass  # Expected

    def test_precision_preservation(self):
        """Test that precision is preserved in operations."""
        v1 = Vector3d(1.0, 0.0, 0.0)
        v2 = Vector3d(0.0, 1.0, 0.0)

        # Multiple operations should preserve precision
        result = v1 + v2
        result = result * 2.0
        result = result / 2.0
        result = result - v2

        # Should get back to v1
        assert abs(result.x - 1.0) < 1e-14
        assert abs(result.y - 0.0) < 1e-14
        assert abs(result.z - 0.0) < 1e-14


class TestVector3dToNumpy:
    """Test Vector3d to numpy conversion."""

    def test_to_numpy_conversion(self):
        """Test conversion to numpy array."""
        v = Vector3d(1.0, 2.0, 3.0)
        arr = v.to_numpy()

        assert isinstance(arr, np.ndarray)
        assert arr.shape == (3,)
        np.testing.assert_allclose(arr, [1.0, 2.0, 3.0])

    def test_numpy_array_operations(self):
        """Test operations with numpy arrays."""
        v = Vector3d(1.0, 2.0, 3.0)
        arr = v.to_numpy()

        # Should be able to do numpy operations
        doubled = arr * 2
        np.testing.assert_allclose(doubled, [2.0, 4.0, 6.0])

        # Check dot product consistency
        v2 = Vector3d(4.0, 5.0, 6.0)
        arr2 = v2.to_numpy()

        numpy_dot = np.dot(arr, arr2)
        vector_dot = v.dot(v2)

        assert abs(numpy_dot - vector_dot) < 1e-10


class TestVector3dSpecialCases:
    """Test special cases and edge conditions."""

    def test_in_place_operations(self):
        """Test in-place operations."""
        v = Vector3d(1.0, 2.0, 3.0)
        original_id = id(v)

        # In-place addition
        v += Vector3d(1.0, 1.0, 1.0)
        assert id(v) == original_id  # Should be same object
        assert v.x == 2.0
        assert v.y == 3.0
        assert v.z == 4.0

    def test_vector_as_boolean(self):
        """Test vector truth value."""
        zero_v = Vector3d(0.0, 0.0, 0.0)
        non_zero_v = Vector3d(1.0, 0.0, 0.0)

        # Test if vectors can be used in boolean context
        # (This depends on implementation)
        if hasattr(zero_v, "__bool__") or hasattr(zero_v, "__nonzero__"):
            assert not zero_v
            assert non_zero_v

    def test_vector_copy_behavior(self):
        """Test vector copying behavior."""
        v1 = Vector3d(1.0, 2.0, 3.0)
        v2 = Vector3d(v1.x, v1.y, v1.z)

        # Should be equal but different objects
        assert v1 == v2
        assert id(v1) != id(v2)

        # Modifying one shouldn't affect the other
        v1.x = 999.0
        assert v2.x == 1.0


if __name__ == "__main__":
    # Simple test runner for when pytest is not available
    import sys

    test_classes = [
        TestVector3dBasics,
        TestVector3dNumerical,
        TestVector3dToNumpy,
        TestVector3dSpecialCases,
    ]

    total_tests = 0
    passed_tests = 0

    for test_class in test_classes:
        instance = test_class()
        methods = [method for method in dir(instance) if method.startswith("test_")]

        for method_name in methods:
            total_tests += 1
            try:
                method = getattr(instance, method_name)
                method()
                print(f"✓ {test_class.__name__}.{method_name}")
                passed_tests += 1
            except Exception as e:
                print(f"✗ {test_class.__name__}.{method_name}: {e}")

    print(f"\nResults: {passed_tests}/{total_tests} tests passed")
    if passed_tests == total_tests:
        print("All tests passed!")
        sys.exit(0)
    else:
        print("Some tests failed!")
        sys.exit(1)
