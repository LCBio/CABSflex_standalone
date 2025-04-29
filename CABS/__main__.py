import argparse
import sys
import importlib.util
import traceback as _tr
from shutil import rmtree

pre_parser = argparse.ArgumentParser(
    description="Custom command-line utility",
    add_help=False  # Disable automatic -h/--help
)
pre_parser.add_argument("job_type", choices=["dock", "flex"], help="Specify the command (dock or flex)")
pre_parser.add_argument("-c", "--config", help="Specify the configuration file")
pre_parser.add_argument("--version", action="store_true", help="Show the version")
pre_parser.add_argument("-h", "--help", action="store_true", help="Show the help")

try:
    import CABS.job
    import CABS.optparser as optparser
    import CABS.logger as logger
    from CABS import __version__, _JUNK
except ImportError:
    cabs_spec = importlib.util.spec_from_file_location("CABS", "./CABS/__init__.py")
    cabs_module = importlib.util.module_from_spec(cabs_spec)
    cabs_spec.loader.exec_module(cabs_module)

    optparser = cabs_module.optparser
    logger = cabs_module.logger
    __version__ = cabs_module.__version__
    _JUNK = cabs_module._JUNK


def run(cmd_line: str, job_type=None):
    if job_type: cmd_line.insert(0, job_type)
    pre_args, remains = pre_parser.parse_known_args(cmd_line)
    job_type = pre_parser.parse_known_args(cmd_line)[0].job_type

    if job_type == 'dock':
        parser = optparser.dock_parser
        task = CABS.job.DockTask
        usage = optparser.dock_usage
    elif job_type == 'flex':    
        parser = optparser.flex_parser
        task = CABS.job.FlexTask
        usage = optparser.flex_usage
    else:
        raise ValueError("Invalid job type specified. Use 'dock' or 'flex'.")

    module_name = 'CABS' + job_type
    
    if pre_args.version:
        print(__version__)
        sys.exit(0)

    if pre_args.help:
        print(parser.format_help())
        sys.exit(0)

    if pre_args.config:
        try:
            remains = optparser.ConfigFileParser(pre_args.config).args + remains
        except IOError:
            logger.exit_program(
                module_name,
                f"Config file: '{pre_args.config}' does not exist."
            )

    if not remains:
        print(usage)
        sys.exit(0)

    config = vars(parser.parse_args(remains))
    job = task(**config)

    try:
        job.run()
    except KeyboardInterrupt:
        logger.critical(module_name, "Interrupted by user.")
    except Exception as e:
        msg = " ".join(str(arg) for arg in e.args)
        logger.exit_program(
            module_name,
            msg,
            traceback=_tr.format_exc(),
            exc=e
        )
    finally:
        logger.close_log()
        for _file in _JUNK:
            rmtree(_file, ignore_errors=True)


def run_dock(cmd_line=sys.argv[1:]):
    run(cmd_line, job_type='dock')


def run_flex(cmd_line=sys.argv[1:]):
    run(cmd_line, job_type='flex')

if __name__ == '__main__':
    run(sys.argv[1:])