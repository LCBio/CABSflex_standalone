#!/usr/bin/env python3
"""
Test runner for CABS-flex test suite.
Can run with or without pytest.
"""

import importlib
from pathlib import Path
import sys
import traceback


def run_tests_without_pytest():
    """Run tests without pytest dependency."""
    test_modules = [
        "test_core_functionality",
        "test_vector3d_comprehensive",
        "test_numerical_utils",
        "test_integration_new",
    ]

    total_tests = 0
    passed_tests = 0
    failed_tests = []

    print("Running CABS-flex test suite (without pytest)")
    print("=" * 50)

    for module_name in test_modules:
        print(f"\nRunning {module_name}...")
        try:
            module = importlib.import_module(module_name)

            # Find test classes
            test_classes = []
            for name in dir(module):
                obj = getattr(module, name)
                if (
                    isinstance(obj, type)
                    and name.startswith("Test")
                    and hasattr(obj, "__init__")
                ):
                    test_classes.append(obj)

            # Run tests in each class
            for test_class in test_classes:
                instance = test_class()
                methods = [
                    method for method in dir(instance) if method.startswith("test_")
                ]

                for method_name in methods:
                    total_tests += 1
                    try:
                        method = getattr(instance, method_name)
                        method()
                        print(f"  ✓ {test_class.__name__}.{method_name}")
                        passed_tests += 1
                    except Exception as e:
                        print(f"  ✗ {test_class.__name__}.{method_name}: {e}")
                        failed_tests.append(
                            f"{module_name}.{test_class.__name__}.{method_name}"
                        )

        except Exception as e:
            print(f"  ERROR importing {module_name}: {e}")
            traceback.print_exc()

    print("\n" + "=" * 50)
    print(f"Results: {passed_tests}/{total_tests} tests passed")

    if failed_tests:
        print("\nFailed tests:")
        for test in failed_tests:
            print(f"  - {test}")

    if passed_tests == total_tests:
        print("🎉 All tests passed!")
        return 0
    else:
        print(f"❌ {len(failed_tests)} tests failed!")
        return 1


def run_tests_with_pytest():
    """Run tests with pytest."""
    import pytest

    # Run pytest on the tests directory
    exit_code = pytest.main(["tests/", "-v", "--tb=short", "--disable-warnings"])
    return exit_code


def main():
    """Main test runner."""
    # Check if pytest is available
    try:
        import pytest

        use_pytest = True
    except ImportError:
        use_pytest = False

    # Check command line arguments
    if "--no-pytest" in sys.argv:
        use_pytest = False
    elif "--pytest" in sys.argv:
        if not use_pytest:
            print("Error: pytest not available but --pytest was requested")
            return 1

    # Add tests directory to path
    tests_dir = Path(__file__).parent
    sys.path.insert(0, str(tests_dir))

    if use_pytest:
        print("Running with pytest...")
        return run_tests_with_pytest()
    else:
        print("Running without pytest...")
        return run_tests_without_pytest()


if __name__ == "__main__":
    sys.exit(main())
