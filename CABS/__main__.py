"""
Modern main entry point for CABS application with proper Python 3 support.
"""

import argparse
import sys
import importlib.util
import traceback as _tr
from pathlib import Path
from shutil import rmtree
from typing import List, Optional, Union

# Create pre-parser for initial argument handling
pre_parser = argparse.ArgumentParser(
    description="CABS: Protein structure prediction and docking simulation tool",
    add_help=False  # Disable automatic -h/--help to handle it manually
)
pre_parser.add_argument(
    "job_type", 
    nargs='?',  # Make job_type optional
    choices=["dock", "flex"], 
    help="Specify the simulation type (dock for protein-protein docking, flex for flexibility analysis)"
)
pre_parser.add_argument(
    "-c", "--config", 
    help="Path to configuration file with simulation parameters"
)
pre_parser.add_argument(
    "--version", 
    action="store_true", 
    help="Show version information and exit"
)
pre_parser.add_argument(
    "-h", "--help", 
    action="store_true", 
    help="Show detailed help information and exit"
)

# Import CABS modules
try:
    import CABS.core.job as job
    import CABS.io.optparser as optparser
    import CABS.io.logger as logger
    from CABS import __version__, _JUNK
except ImportError as e:
    print(f"Error: Could not import CABS modules: {e}", file=sys.stderr)
    print("Make sure you're running from the correct directory and CABS is properly installed.", file=sys.stderr)
    sys.exit(1)


def run(cmd_line: List[str], job_type: Optional[str] = None) -> None:
    """
    Main execution function for CABS simulations.
    
    Args:
        cmd_line: Command line arguments
        job_type: Override job type (dock or flex)
    """
    if job_type:
        cmd_line.insert(0, job_type)
    
    try:
        pre_args, remains = pre_parser.parse_known_args(cmd_line)
        determined_job_type = pre_args.job_type or job_type
    except SystemExit:
        # argparse calls sys.exit on error, we catch it to handle gracefully
        return

    # If no job type is determined, show general help or combined help
    if not determined_job_type:
        if pre_args.help:
            # Show combined help for both dock and flex
            print("CABS: Protein structure prediction and docking simulation tool")
            print(f"Version: {__version__}")
            print("\nCABS provides two main simulation modes:")
            print("\n" + "="*60)
            print("CABSflex - Flexibility Analysis")
            print("="*60)
            print(optparser.flex_parser.format_help())
            print("\n" + "="*60)
            print("CABSdock - Protein-Protein Docking")
            print("="*60)
            print(optparser.dock_parser.format_help())
        else:
            # Show brief usage
            print("CABS: Protein structure prediction and docking simulation tool")
            print(f"Version: {__version__}")
            print("\nUsage:")
            print("  CABSflex [options] <input.pdb>    # For flexibility analysis")
            print("  CABSdock [options] <input.pdb>    # For protein-protein docking")
            print("\nFor detailed help:")
            print("  python -m CABS --help           # Show all options")
            print("  CABSflex --help                 # Show CABSflex options only")
            print("  CABSdock --help                 # Show CABSdock options only")
        sys.exit(0)

    # Configure parsers and tasks based on job type
    if determined_job_type == 'dock':
        parser = optparser.dock_parser
        task = job.DockTask
        usage = optparser.dock_usage
    elif determined_job_type == 'flex':    
        parser = optparser.flex_parser
        task = job.FlexTask
        usage = optparser.flex_usage
    else:
        print(f"Error: Invalid job type '{determined_job_type}'. Use 'dock' or 'flex'.", file=sys.stderr)
        sys.exit(1)

    module_name = f'CABS{determined_job_type}'
    
    # Handle version request
    if pre_args.version:
        print(f"CABS version {__version__}")
        sys.exit(0)

    # Handle help request for specific job type
    if pre_args.help:
        print(parser.format_help())
        sys.exit(0)

    # Handle configuration file
    if pre_args.config:
        config_path = Path(pre_args.config)
        if not config_path.exists():
            logger.exit_program(
                module_name,
                f"Configuration file '{config_path}' does not exist."
            )
        try:
            remains = optparser.ConfigFileParser(str(config_path)).args + remains
        except Exception as e:
            logger.exit_program(
                module_name,
                f"Error reading configuration file '{config_path}': {e}"
            )

    # Show usage if no arguments provided
    if not remains:
        print(usage)
        sys.exit(0)

    # Parse arguments and create job
    try:
        config = vars(parser.parse_args(remains))
        cabs_job = task(**config)
    except Exception as e:
        logger.exit_program(
            module_name,
            f"Error creating job: {e}",
            traceback=_tr.format_exc(),
            exc=e
        )

    # Execute the job
    try:
        cabs_job.run()
    except KeyboardInterrupt:
        logger.critical(module_name, "Interrupted by user.")
        sys.exit(130)  # Standard exit code for Ctrl+C
    except Exception as e:
        msg = " ".join(str(arg) for arg in e.args) if e.args else str(e)
        logger.exit_program(
            module_name,
            msg,
            traceback=_tr.format_exc(),
            exc=e
        )
    finally:
        # Clean up resources
        logger.close_log()
        for file_path in _JUNK:
            try:
                if Path(file_path).is_dir():
                    rmtree(file_path, ignore_errors=True)
                else:
                    Path(file_path).unlink(missing_ok=True)
            except Exception:
                pass  # Ignore cleanup errors


def run_dock(cmd_line: Optional[List[str]] = None) -> None:
    """
    Entry point for CABS docking simulations.
    
    Args:
        cmd_line: Command line arguments, defaults to sys.argv[1:]
    """
    if cmd_line is None:
        cmd_line = sys.argv[1:]
    run(cmd_line, job_type='dock')


def run_flex(cmd_line: Optional[List[str]] = None) -> None:
    """
    Entry point for CABS flexibility analysis.
    
    Args:
        cmd_line: Command line arguments, defaults to sys.argv[1:]
    """
    if cmd_line is None:
        cmd_line = sys.argv[1:]
    run(cmd_line, job_type='flex')


if __name__ == '__main__':
    run(sys.argv[1:])