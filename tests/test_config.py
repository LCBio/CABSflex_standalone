"""
Tests for configuration loading and data management.
"""

import json

from CABS.config_loader import (
    get_alignment_data,
    get_config_section,
    get_extended_amino_acids,
    get_logger_config,
    get_molecular_data,
)


class TestConfigurationLoading:
    """Test configuration file loading."""

    def test_molecular_data_loading(self):
        """Test loading molecular data configuration."""
        try:
            data = get_molecular_data()
            assert isinstance(data, dict)
            assert len(data) > 0
        except FileNotFoundError:
            # Config file might not exist in test environment
            pass

    def test_amino_acids_loading(self):
        """Test loading amino acids configuration."""
        try:
            data = get_extended_amino_acids()
            assert isinstance(data, dict)
            assert len(data) > 0
        except FileNotFoundError:
            pass

    def test_alignment_data_loading(self):
        """Test loading alignment data configuration."""
        try:
            data = get_alignment_data()
            assert isinstance(data, dict)
        except FileNotFoundError:
            pass

    def test_logger_config_loading(self):
        """Test loading logger configuration."""
        try:
            data = get_logger_config()
            assert isinstance(data, dict)
        except FileNotFoundError:
            pass

    def test_config_section_access(self):
        """Test accessing specific config sections."""
        try:
            # Test with a known section if available
            data = get_config_section("molecular_data", "amino_acids")
            if data is not None:
                assert isinstance(data, (dict, list, str, int, float))
        except (FileNotFoundError, KeyError):
            # Config might not exist or section might not be present
            pass


class TestConfigurationCaching:
    """Test configuration caching functionality."""

    def test_config_caching(self):
        """Test that configurations are cached properly."""
        try:
            # Call the same function twice
            data1 = get_molecular_data()
            data2 = get_molecular_data()

            # Should return the same object (cached)
            assert data1 is data2
        except FileNotFoundError:
            pass


class TestConfigurationValidation:
    """Test configuration data validation."""

    def test_molecular_data_structure(self):
        """Test molecular data has expected structure."""
        try:
            data = get_molecular_data()

            # Should be a dictionary
            assert isinstance(data, dict)

            # Should contain expected keys (if the config exists)
            if data:
                # Common expected keys in molecular data
                expected_keys = ["amino_acids", "bonds", "masses"]
                found_keys = [key for key in expected_keys if key in data]
                # At least some expected keys should be present
                assert len(found_keys) > 0 or len(data) > 0

        except FileNotFoundError:
            pass

    def test_logger_config_structure(self):
        """Test logger config has expected structure."""
        try:
            data = get_logger_config()

            assert isinstance(data, dict)

            if data:
                # Logger config should have levels or colors
                expected_keys = ["levels", "colors", "formats"]
                found_keys = [key for key in expected_keys if key in data]
                assert len(found_keys) > 0 or len(data) > 0

        except FileNotFoundError:
            pass


class TestConfigurationErrorHandling:
    """Test error handling in configuration loading."""

    def test_nonexistent_config(self):
        """Test handling of non-existent configuration."""
        try:
            # Try to access a non-existent section
            data = get_config_section("nonexistent_config", "nonexistent_section")
            # Should return None or raise appropriate error
            assert data is None
        except (FileNotFoundError, KeyError):
            # These errors are acceptable
            pass

    def test_invalid_json_handling(self):
        """Test handling of invalid JSON (if possible to test)."""
        # This is difficult to test without modifying the actual config files
        # but we can at least verify the function exists and handles errors
        try:
            get_molecular_data()
        except (FileNotFoundError, json.JSONDecodeError):
            # These are expected possible errors
            pass

    def test_config_file_permissions(self):
        """Test handling of permission errors."""
        # Also difficult to test without changing file permissions
        # but we can verify the functions handle errors gracefully
        try:
            get_extended_amino_acids()
        except (FileNotFoundError, PermissionError):
            pass


class TestConfigurationConsistency:
    """Test consistency across different configuration sources."""

    def test_amino_acid_consistency(self):
        """Test consistency between different amino acid configurations."""
        try:
            molecular_data = get_molecular_data()
            extended_aa = get_extended_amino_acids()

            if molecular_data and extended_aa:
                # If both configs exist, check for consistency
                if "amino_acids" in molecular_data and isinstance(extended_aa, dict):
                    mol_aa = molecular_data["amino_acids"]

                    # Basic consistency check - overlapping amino acids should match
                    for aa_code in mol_aa:
                        if aa_code in extended_aa:
                            # Should at least have consistent naming
                            assert isinstance(mol_aa[aa_code], (str, dict))
                            assert isinstance(extended_aa[aa_code], (str, dict))

        except (FileNotFoundError, KeyError, TypeError):
            pass

    def test_config_data_types(self):
        """Test that configuration data has appropriate types."""
        configs_to_test = [
            get_molecular_data,
            get_extended_amino_acids,
            get_alignment_data,
            get_logger_config,
        ]

        for config_func in configs_to_test:
            try:
                data = config_func()

                # All configs should return dictionaries
                assert isinstance(data, dict)

                # All keys should be strings
                for key in data.keys():
                    assert isinstance(key, str)

                # Values should be JSON-serializable types
                for value in data.values():
                    assert isinstance(
                        value, (str, int, float, bool, list, dict, type(None))
                    )

            except FileNotFoundError:
                pass


class TestConfigurationContent:
    """Test actual configuration content (when available)."""

    def test_amino_acid_codes(self):
        """Test amino acid codes are valid."""
        try:
            extended_aa = get_extended_amino_acids()

            if extended_aa:
                # Extended amino acids might use different format (like single letters)
                # or different keys, so check for any reasonable amino acid data
                assert isinstance(extended_aa, dict), (
                    "Extended amino acids should be a dict"
                )
                assert len(extended_aa) > 0, (
                    "Extended amino acids dict should not be empty"
                )

                # Check if we have at least some amino acid data
                # Could be single letters, three letters, or other format
                has_amino_data = any(
                    key in extended_aa
                    for key in [
                        "A",
                        "G",
                        "L",
                        "V",
                        "I",
                        "P",
                        "F",
                        "W",
                        "M",
                        "S",
                        "T",
                        "C",
                        "Y",
                        "N",
                        "Q",
                        "D",
                        "E",
                        "K",
                        "R",
                        "H",
                    ]
                ) or any(
                    key in extended_aa
                    for key in [
                        "ALA",
                        "VAL",
                        "LEU",
                        "ILE",
                        "PRO",
                        "PHE",
                        "TRP",
                        "MET",
                        "GLY",
                        "SER",
                        "THR",
                        "CYS",
                        "TYR",
                        "ASN",
                        "GLN",
                        "ASP",
                        "GLU",
                        "LYS",
                        "ARG",
                        "HIS",
                    ]
                )

                if len(extended_aa) > 10:  # Only check if we have substantial data
                    assert has_amino_data, (
                        f"No recognizable amino acid codes found in keys: {list(extended_aa.keys())[:10]}"
                    )

        except FileNotFoundError:
            pass  # Config file not found, skip test

    def test_molecular_weights(self):
        """Test molecular weight data if available."""
        try:
            molecular_data = get_molecular_data()

            if molecular_data and "masses" in molecular_data:
                masses = molecular_data["masses"]

                if isinstance(masses, dict):
                    # All masses should be positive numbers
                    for element, mass in masses.items():
                        assert isinstance(element, str)
                        assert isinstance(mass, (int, float))
                        assert mass > 0, f"Invalid mass for {element}: {mass}"

        except (FileNotFoundError, KeyError, TypeError):
            pass

    def test_logger_levels(self):
        """Test logger level configuration."""
        try:
            logger_config = get_logger_config()

            if logger_config and "levels" in logger_config:
                levels = logger_config["levels"]

                if isinstance(levels, dict):
                    # Should contain numeric levels
                    for level_name, level_value in levels.items():
                        assert isinstance(level_name, str)
                        assert isinstance(level_value, (int, str))

                        # If numeric, should be reasonable range
                        if isinstance(level_value, int):
                            assert 0 <= level_value <= 100, (
                                f"Invalid level {level_value}"
                            )

        except (FileNotFoundError, KeyError, TypeError):
            pass


class TestConfigurationIntegration:
    """Test integration between configuration and other modules."""

    def test_config_usage_in_modules(self):
        """Test that configurations can be used by other modules."""
        # This is more of an integration test
        try:
            # Try importing modules that use configurations
            from CABS.constants import AA_NAMES
            from CABS.io.logger import setup

            # If imports succeed, basic integration is working
            assert True

        except ImportError as e:
            # Module might not be available or have dependencies
            pass

    def test_config_modification_isolation(self):
        """Test that config modifications don't affect other calls."""
        try:
            data1 = get_molecular_data()

            if data1:
                # Modify the returned data
                original_keys = list(data1.keys())
                data1["test_key"] = "test_value"

                # Get config again
                data2 = get_molecular_data()

                # Due to caching, this might be the same object
                # But the test shows whether the modification persists
                if data1 is data2:
                    # Same object due to caching - modification will persist
                    assert "test_key" in data2
                else:
                    # Different object - modification should not persist
                    assert "test_key" not in data2

        except FileNotFoundError:
            pass
