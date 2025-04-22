import argparse
import sys
import importlib.util
import traceback as _tr
from shutil import rmtree

pre_parser = argparse.ArgumentParser(
    description="Custom command-line utility",
    add_help=False  # Disable automatic -h/--help
)
pre_parser.add_argument("cabs_cmd", choices=["dock", "flex"], help="Specify the command (dock or flex)")
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
    if not job_type: job_type = pre_parser.parse_known_args()[0].cabs_cmd
    pre_args, remains = pre_parser.parse_known_args(cmd_line)

    if job_type not in ['dock', 'flex']:
        raise ValueError(f"Invalid command: {job_type}")
    
    parsers = {
        'dock': optparser.dock_parser,
        'flex': optparser.flex_parser
    }
    tasks = {
        'dock': CABS.job.DockTask,
        'flex': CABS.job.FlexTask
    }
    usage = {
        'dock': optparser.dock_usage,
        'flex': optparser.flex_usage
    }

    parser = parsers[job_type]
    task = tasks[job_type]
    usage = usage[job_type]
    module_name = 'CABS' + job_type


    if pre_args.version:
        print(__version__)
        sys.exit(0)
    elif pre_args.help:
        _help = parser.format_help()
        print(_help)
        sys.exit(0)
    elif pre_args.config:
        try:
            remains = optparser.ConfigFileParser(pre_args.config).args + remains
        except IOError:
            logger.exit_program(
                module_name,
                f'Config file: \'{pre_args.config}\' does not exist.'
            )
    elif not len(remains):
        print(usage)
        sys.exit(0)

    config = vars(parser.parse_args(remains))

    job = task(**config)
    try:
        job.run()
    except KeyboardInterrupt:
        logger.critical(module_name, 'Interrupted by user.')
    except Exception as e:
        msg = str(e)
        for a in e.args:
            try:
                msg += ' ' + str(a)
            except ValueError:
                pass
        logger.exit_program(module_name, msg, traceback=_tr.format_exc(), exc=e)
    finally:    
        logger.close_log()
        for _file in _JUNK:
            rmtree(_file, ignore_errors=True)


def run_dock(cmd_line=sys.argv[1:]):
    run(cmd_line, job_type='dock')


def run_flex(cmd_line=sys.argv[1:]):
    run(cmd_line, job_type='flex')

if __name__ == '__main__':
    # not sure if its necessary to prepars args before run(), it is done again in run()
    # args, remains = pre_parser.parse_known_args()
    run(sys.argv[1:])
