"""
CABS: Protein structure prediction and docking simulation package.

A modern Python 3 implementation of the CABS coarse-grained protein model
for protein folding, flexibility analysis, and protein-protein docking simulations.
"""

from pathlib import Path
from typing import List

__version__ = "3.0.2"
__title__ = "CABS"
__author__ = "Laboratory of Computational Biology"
__email__ = "k.wroblewski7@uw.edu.pl"
__license__ = "MIT"

# Global list for cleanup of temporary files and directories
_JUNK: List[Path] = []


def cleanup_junk() -> None:
    """Clean up temporary files and directories."""
    for path in _JUNK:
        try:
            if path.is_dir():
                import shutil

                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink(missing_ok=True)
        except Exception:
            pass  # Ignore cleanup errors
    _JUNK.clear()


# Ensure cleanup on module import
import atexit

atexit.register(cleanup_junk)
