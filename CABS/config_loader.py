"""
Configuration loader for CABS package.

This module provides centralized loading of configuration data from JSON files.
It caches loaded configurations to avoid repeated file I/O operations.
"""

from functools import lru_cache
import json
from typing import Any, Dict, Union
import warnings

try:
    from importlib.resources import as_file, files
except ImportError:
    # Suppress pkg_resources deprecation warning until we fully migrate to importlib.resources
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=DeprecationWarning, module="pkg_resources"
        )
        from pkg_resources import resource_filename


@lru_cache(maxsize=None)
def _load_json_config(config_name: str) -> Dict[str, Any]:
    """
    Load JSON configuration file with caching.

    Args:
        config_name: Name of the JSON configuration file (without .json extension)

    Returns:
        Dictionary containing the configuration data

    Raises:
        FileNotFoundError: If the configuration file is not found
        json.JSONDecodeError: If the file contains invalid JSON
    """
    config_file = f"{config_name}.json"

    try:
        # Try modern importlib.resources first
        try:
            with as_file(files("CABS") / "data" / config_file) as data_file:
                with open(data_file, encoding="utf-8") as f:
                    return json.load(f)
        except (ImportError, AttributeError):
            # Fallback to pkg_resources with warning suppression
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", category=DeprecationWarning, module="pkg_resources"
                )
                data_file = resource_filename("CABS", f"data/{config_file}")
            with open(data_file, encoding="utf-8") as f:
                return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Configuration file '{config_file}' not found in CABS/data/"
        )
    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(
            f"Invalid JSON in configuration file '{config_file}': {e.msg}", e.doc, e.pos
        )


def get_molecular_data() -> Dict[str, Any]:
    """Get molecular data configuration."""
    return _load_json_config("molecular_data")


def get_extended_amino_acids() -> Dict[str, Any]:
    """Get extended amino acids configuration."""
    return _load_json_config("extended_amino_acids")


def get_alignment_data() -> Dict[str, Any]:
    """Get alignment data configuration."""
    return _load_json_config("alignment_data")


def get_logger_config() -> Dict[str, Any]:
    """Get logger configuration."""
    return _load_json_config("logger_config")


def get_config_section(config_name: str, section: str) -> Any:
    """
    Get a specific section from a configuration file.

    Args:
        config_name: Name of the JSON configuration file (without .json extension)
        section: Section name within the configuration

    Returns:
        The requested configuration section

    Raises:
        KeyError: If the section is not found in the configuration
    """
    config = _load_json_config(config_name)
    if section not in config:
        raise KeyError(
            f"Section '{section}' not found in configuration '{config_name}'"
        )
    return config[section]


def clear_cache() -> None:
    """Clear the configuration cache. Useful for testing or reloading configs."""
    _load_json_config.cache_clear()


# Convenience functions for commonly used configurations
def get_sidecnt() -> Dict[str, list]:
    """Get side chain coordinates data."""
    return get_config_section("molecular_data", "sidecnt")


def get_cabs_ss() -> Dict[str, int]:
    """Get CABS secondary structure mapping."""
    return get_config_section("molecular_data", "cabs_ss")


def get_cabs_ss_reverse() -> Dict[str, str]:
    """Get reverse CABS secondary structure mapping."""
    return get_config_section("molecular_data", "cabs_ss_reverse")


def get_aa_names() -> Dict[str, str]:
    """Get amino acid names mapping (single to three letter)."""
    return get_config_section("molecular_data", "aa_names")


def get_aa_sub_names() -> Dict[str, str]:
    """Get amino acid substitution names mapping (three to single letter)."""
    return get_config_section("molecular_data", "aa_sub_names")


def get_aa_sub_names_extended() -> Dict[str, str]:
    """Get extended amino acid substitution names mapping."""
    return get_config_section("extended_amino_acids", "aa_sub_names_extended")


def get_blosum62_indices() -> Dict[str, int]:
    """Get BLOSUM62 matrix indices."""
    return get_config_section("alignment_data", "blosum62_indices")


def get_blosum62_matrix() -> list:
    """Get BLOSUM62 substitution matrix."""
    return get_config_section("alignment_data", "blosum62_matrix")


def get_log_colors() -> Dict[str, str]:
    """Get log color codes."""
    return get_config_section("logger_config", "log_colors")


def get_log_levels() -> Dict[str, str]:
    """Get log level labels."""
    return get_config_section("logger_config", "log_levels")


def get_color_prefix() -> Dict[str, str]:
    """Get colored log level prefixes."""
    return get_config_section("logger_config", "color_prefix")


def get_type_dispatch() -> Dict[str, str]:
    """Get type dispatch mapping."""
    return get_config_section("logger_config", "type_dispatch")


def get_cabs_constants() -> Dict[str, Any]:
    """Get CABS constants configuration."""
    return _load_json_config("cabs_constants")


# Convenience functions for CABS constants
def get_cabs_lattice_defaults() -> Dict[str, Union[float, list]]:
    """Get CABS lattice default values."""
    return get_config_section("cabs_constants", "cabs_lattice_defaults")


def get_sc_modeling_thresholds() -> Dict[str, float]:
    """Get side chain modeling thresholds."""
    return get_config_section("cabs_constants", "sc_modeling_thresholds")


def get_pdb_output_options() -> Dict[str, str]:
    """Get PDB output options mapping."""
    return get_config_section("cabs_constants", "pdb_output_options")


def get_bfac_output_options() -> Dict[str, str]:
    """Get beta factor output options mapping."""
    return get_config_section("cabs_constants", "bfac_output_options")


def get_csv_output_options() -> Dict[str, str]:
    """Get CSV output options mapping."""
    return get_config_section("cabs_constants", "csv_output_options")


def get_valid_letters() -> Dict[str, str]:
    """Get valid letters for different output types."""
    return get_config_section("cabs_constants", "valid_letters")


def get_protein_restraints_modes() -> list:
    """Get protein restraints modes."""
    return get_config_section("cabs_constants", "protein_restraints_modes")


def get_protein_category_modes() -> list:
    """Get protein category modes."""
    return get_config_section("cabs_constants", "protein_category_modes")


def get_unleashed_aliases() -> list:
    """Get unleashed aliases."""
    return get_config_section("cabs_constants", "unleashed_aliases")


def get_allowed_aa_methods() -> list:
    """Get allowed amino acid reconstruction methods."""
    return get_config_section("cabs_constants", "allowed_aa_methods")


def get_cabs_files() -> list:
    """Get CABS simulation files list."""
    return get_config_section("cabs_constants", "cabs_files")


def get_dssp_ss_mapping() -> Dict[str, str]:
    """Get DSSP secondary structure mapping."""
    return get_config_section("cabs_constants", "dssp_ss_mapping")


def get_default_colors() -> list:
    """Get default colors for plots."""
    return get_config_section("cabs_constants", "default_colors")


def get_contact_map_constants() -> Dict[str, Union[int, float]]:
    """Get contact map constants."""
    return get_config_section("cabs_constants", "contact_map_constants")


def get_peptide_replacements() -> Dict[str, str]:
    """Get peptide replacement patterns."""
    return get_config_section("cabs_constants", "peptide_replacements")


def get_output_directories() -> Dict[str, str]:
    """Get output directory names."""
    return get_config_section("cabs_constants", "output_directories")


def get_default_filenames() -> Dict[str, str]:
    """Get default file names."""
    return get_config_section("cabs_constants", "default_filenames")


def get_calculation_constants() -> Dict[str, Union[int, float]]:
    """Get calculation constants and thresholds."""
    return get_config_section("cabs_constants", "calculation_constants")


def get_math_constants() -> Dict[str, Union[int, float]]:
    """Get mathematical constants."""
    return get_config_section("cabs_constants", "math_constants")


def get_string_patterns() -> Dict[str, str]:
    """Get string patterns and templates."""
    return get_config_section("cabs_constants", "string_patterns")


def get_system_constants() -> Dict[str, str]:
    """Get system and platform constants."""
    return get_config_section("cabs_constants", "system_constants")


def get_model_constants() -> Dict[str, str]:
    """Get model file constants."""
    return get_config_section("cabs_constants", "model_constants")


def get_file_extensions() -> Dict[str, Union[str, list]]:
    """Get file extensions configuration."""
    return get_config_section("cabs_constants", "file_extensions")


def get_default_peptide_ss() -> str:
    """Get default peptide secondary structure."""
    return get_config_section("cabs_constants", "default_peptide_ss")


def get_error_messages() -> Dict[str, str]:
    """Get error messages."""
    return get_config_section("cabs_constants", "error_messages")


def get_default_values() -> Dict[str, Union[int, float, str, None]]:
    """Get default parameter values."""
    return get_config_section("cabs_constants", "default_values")


def get_nsp3_constants() -> Dict[str, str]:
    """Get NetSurfP-3.0 constants."""
    return get_config_section("cabs_constants", "nsp3_constants")


def get_config_header() -> str:
    """Get configuration file header template."""
    return get_config_section("cabs_constants", "config_header")

def get_cg2all_env_prefix() -> str:
    """Reads custom path information saved by the installer script."""
    return get_config_section("cabs_paths", "cg2all_env_prefix")
