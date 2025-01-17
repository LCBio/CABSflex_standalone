import os
import sys
from threading import Thread
from threading import Event
from sys import stderr
from time import time, strftime, gmtime, sleep
import textwrap
import platform

_name = 'Logger'

colors = {
    'blue': '\033[94m',
    'yellow': '\033[93m',
    'green': '\033[92m',
    'red': '\033[91m',
    'light_blue': '\033[96m',
    'purple': '\033[95m',
    'end': '\033[0m',
}

log_levels = {
    0: '[CRITICAL]',
    1: '[WARNING]',
    2: '[INFO]',
    3: '[OUT FILES]',
    4: '[DEBUG]'
}

color_prefix = {
    0: colors['red'] + '[CRITICAL]' + colors['end'],
    1: colors['yellow'] + '[WARNING]' + colors['end'],
    2: colors['blue'] + '[INFO]' + colors['end'],
    3: colors['purple'] + '[OUT FILES]' + colors['end'],
    4: colors['green'] + '[DEBUG]' + colors['end'],
}

_init_time = time()
_log_level = 2
_color = True
_stream = sys.stderr
_line_format = '%-20s %-19s%-75s %s\n'
_first_line_format = '%-20s %-19s%-75s \n'
_middle_line_format = '%-22s%-75s \n'
_last_line_format = '%-22s%-75s %s\n'
_line_break = 76
_remote = False
_prefix = color_prefix
_save_sec_str = False
_save_dssp = False
_save_bfac = False
_save_restraints = False
_save_plddt = False
_save_category = False
_save_contact = False

def setup(log_level=2, remote=False, work_dir='', save_sec_str=False, save_dssp=False, save_bfac=False, save_restraints=False,
          save_plddt=False, save_category=False, save_contact=False):
    global _log_level, _color, _stream, _remote, _line_break, _prefix, _save_sec_str, _save_dssp, _save_bfac,\
        _save_restraints,_save_plddt, _save_category, _save_contact
    global _line_format, _middle_line_format, _first_line_format, _last_line_format
    _remote = remote
    if _remote or not sys.stderr.isatty():
        _color = False
        _line_format = '%-12s %-10s%-75s %s\n'
        _first_line_format = '%-12s %-10s%-75s \n'
        _middle_line_format = '%-22s %-75s \n'
        _last_line_format = '%-22s %-75s %s\n'
        _prefix = log_levels
    if _remote:
        _log_path = os.path.join(work_dir, 'CABS.log')
        try:
            _stream = open(_log_path, 'a+')
            _stream.write('#' * 110 + "\n")
        except IOError:
            try:
                os.makedirs(work_dir)
                _stream = open(_log_path, 'a+')
            except OSError:
                warning(
                    module_name=_name,
                    msg='Could not open a log file at %s. Writing to standard error instead.' % _log_path
                )
                raise
    if platform.system() == "Windows":
        _color = False
        _prefix = log_levels
    _log_level = log_level
    info(_name, 'Verbosity set to: ' + str(log_level) + ' - ' + log_levels[log_level])
    _save_sec_str = save_sec_str
    _save_dssp = save_dssp
    _save_bfac = save_bfac
    _save_restraints = save_restraints
    _save_plddt = save_plddt
    _save_category = save_category
    _save_contact = save_contact


def close_log():
    if _stream is not sys.stderr:
        _stream.close()


def log_files():
    """
    :return: True if verbosity is high enough to save extra output (LOG FILE level)
    """
    return _log_level >= 3


def output_sec_str():
    """
    :return: True if flag --sec-str-output was set
    """
    return _save_sec_str


def output_dssp():
    """
    :return: True if flag --dssp-output was set
    """
    return _save_dssp


def output_bfac():
    """
    :return: True if flag --bfac-output was set
    """
    return _save_bfac


def output_restraints():
    """
    :return: True if flag --restraints-output was set
    """
    return _save_restraints


def output_plddt():
    """
    :return: True if flag --plddt-output was set
    """
    return _save_plddt


def output_category():
    """
    :return: True if flag --category-output was set
    """
    return _save_category


def output_contact():
    """
    :return: True if flag --contact-output was set
    """
    return _save_contact


def coloring(color_name='light_blue', msg=''):
    if _color:
        return colors[color_name] + msg + colors['end']
    return msg


def log(module_name='MISC', msg='Processing ', l_level=2, out=None):
    if out is None:
        out = _stream
    if l_level <= _log_level:
        t = gmtime(time() - _init_time)
        if len(msg) < _line_break:
            msg = _line_format % (
                _prefix[l_level],
                coloring(msg=module_name + ':', color_name='light_blue'),
                msg,
                strftime('(%H:%M:%S)', t)
            )
            out.write(msg)
            out.flush()
        else:
            if isinstance(msg, bytes):
                msg = msg.decode('utf-8')
            else:
                msg = msg
            lines = textwrap.wrap(msg, width=_line_break - 1)
            first_line = _first_line_format % (
                _prefix[l_level],
                coloring(msg=module_name + ':', color_name='light_blue'),
                lines[0]
            )
            out.write(first_line)
            for lineNumber in range(1, len(lines) - 1):
                line = _middle_line_format % (' ', lines[lineNumber])
                out.write(line)
            final_line = _last_line_format % (' ', lines[-1], strftime('(%H:%M:%S)', t))
            out.write(final_line)
            out.flush()


def critical(module_name='_name', msg=''):
    log(module_name=module_name, msg=msg, l_level=0)


def warning(module_name='_name', msg=''):
    log(module_name=module_name, msg=msg, l_level=1)


def info(module_name='_name', msg=''):
    log(module_name=module_name, msg=msg, l_level=2)


def log_file(module_name='_name', msg=''):
    log(module_name=module_name, msg=msg, l_level=3)


def debug(module_name='_name', msg=''):
    log(module_name=module_name, msg=msg, l_level=4)


def to_file(filename='', content='', msg='', allow_err=True, traceback=True):
    """
    :param filename: path for the file to be saved
    :param content: a string to be saved (be careful not to pass a string that is too long
    :param msg: optional: a message to be logged
    :param allow_err: if True a warning will be logged on OSError, if False program exit call will be made
    :param traceback: if True an Exception will be raised on exit call
    :return:
    """
    if filename and content:
        try:
            if os.path.isfile(filename):
                log_file(module_name=_name, msg='Overwriting %s' % filename)
            with open(filename, 'w') as f:
                f.write(content)
        except IOError:
            if allow_err:
                warning(module_name=_name, msg='IOError while writing to: %s' % filename)
            else:
                exit_program(module_name=_name, msg='IOError while writing to: %s' % filename, traceback=traceback)
    if msg:
        log_file(module_name=_name, msg=msg)


def exit_program(module_name=_name, msg='Shutting down', traceback=None, exc=None):
    """
    Exits the program, depending on options might do it quietly, raise, or print traceback.

    :param module_name: string with  the calling module's name
    :param msg: string, message to be printed when the program exits
    :param traceback: bool, if log level is high traceback will be printed
    :param exc: a specific exception passed by the caller
    :return: None
    """

    if exc:
        _msg = '%s: %s' % (msg, exc)
    else:
        _msg = msg
    critical(module_name=module_name, msg=_msg)
    if _log_level > 3 and traceback:
        _stream.write(traceback)
    sys.exit(1)


class ProgressBar:
    """This class assumes a manual call to done() will be made to exit the bar"""
    WIDTH = 65
    FORMAT = '%-20s %-19s[%s] %.1f%%\r'
    BAR0 = ' '
    BAR1 = '#'

    def __init__(self, total=100, module_name='', job_name='', out=stderr, delay=0, start_msg=''):
        if _log_level >= 2 and not _remote:
            self.stream = out
        else:
            self.stream = open(os.devnull, 'w')
        self.total = total
        self.current = 0.
        self.job_name = job_name
        self.is_done = False
        self.module_name = module_name
        self.prefix = _prefix[2]
        if start_msg:
            self.stream.write(coloring(msg=start_msg) + '\n')
        if self.job_name:
            log(module_name=self.module_name, msg=self.job_name + ' running...', out=self.stream)
        self.start_time = time()
        self.update()
        sleep(delay)

    def write(self):
        percent = 1.0 * self.current / self.total
        num = int(self.WIDTH * percent)
        percent = round(100. * percent, 1)
        bar = self.BAR1 * num + self.BAR0 * (self.WIDTH - num)
        self.stream.write(self.FORMAT % (self.prefix, coloring(msg=self.module_name + ':'), bar, percent))
        self.stream.flush()

    def update(self, state=-1.):
        if state < 0:
            self.current += 1.
        else:
            self.current = state
        if self.current >= self.total:
            return True
        self.write()
        return False

    def finish(self):
        if self.current < self.total:
            for i in range(int(self.total - self.current)):
                self.update()
                sleep(0.001)

    def done(self, show_time=True):
        if not self.is_done:
            self.finish()
            self.stream.write(' ' * 80 + '\r')
            if show_time:
                t = gmtime(time() - self.start_time)
                log(
                    module_name=self.module_name,
                    msg=self.job_name + ' done in %s' % strftime('%H:%M:%S', t),
                    out=self.stream
                )
            self.stream.flush()
            self.is_done = True


class CabsObserver(Thread):

    def __init__(self, interval=0.5, progress_file=None, job_name='CABS simulation', msg=''):
        Thread.__init__(self)
        self.exit_event = Event()
        self.interval = interval
        self.progress_bar = ProgressBar(module_name='CABS', job_name=job_name, start_msg=msg)
        self.progress_file = progress_file
        self.daemon = True  # In case main program ends abruptly
        self.start()

    def exit(self):
        self.progress_bar.done()
        self.exit_event.set()

    def run(self):
        while not self.exit_event.is_set():
            if self.progress_bar.update(self.status()):
                self.exit()
            sleep(self.interval)

    def status(self):
        try:
            with open(self.progress_file, 'rb') as f:
                progress = float(f.read())
        except (IOError, ValueError):
            progress = 0
        return progress
