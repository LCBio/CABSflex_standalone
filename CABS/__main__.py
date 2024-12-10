import argparse
import sys
import importlib.util
import traceback as _tr
from shutil import rmtree


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


def run(cabs_cmd: str, cmd_line: str):
    if cabs_cmd not in ['dock', 'flex']:
        raise ValueError(f"Invalid command: {cabs_cmd}")

    module_name = 'CABS' + cabs_cmd
    parser = getattr(optparser, cabs_cmd + '_parser')

    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('-c', '--config')
    pre_parser.add_argument('--version', action='store_true')
    pre_parser.add_argument('-h', '--help', action='store_true')

    pre_args, remains = pre_parser.parse_known_args(cmd_line)
    # return pre_args, remains
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
        print(getattr(optparser, cabs_cmd + '_usage'))
        sys.exit(0)

    config = vars(parser.parse_args(remains))

    task = getattr(CABS.job, cabs_cmd.title() + 'Task')  
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
    run('dock', cmd_line)


def run_flex(cmd_line=sys.argv[1:]):
    run('flex', cmd_line)



def pre_parser():
    pre_parser = argparse.ArgumentParser(description="Custom command-line utility")
    pre_parser.add_argument("cabs_cmd", choices=["dock", "flex"], help="Specify the command (dock or flex)")
    pre_parser.add_argument("-c", "--config", help="Specify the configuration file")
    pre_parser.add_argument("--version", action="store_true", help="Show the version")
    pre_parser.add_argument("-h", "--help", action="store_true", help="Show the help")
    args = pre_parser.parse_args()
    return args


if __name__ == '__main__':
    args = pre_parser()
    run(args.cabs_cmd, cmd_line=sys.argv[1:])
