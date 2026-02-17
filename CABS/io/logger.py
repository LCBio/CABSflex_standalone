"""Logging module for CABS with Python 3 support."""

from enum import Enum
import os
from pathlib import Path
import platform
import sys
from sys import stderr
import textwrap
from threading import Event, Thread
from time import gmtime, sleep, strftime, time
from typing import Optional, TextIO, Union

from CABS.config_loader import get_color_prefix, get_log_colors, get_log_levels

__all__ = [
    "CabsObserver",
    "ProgressBar",
    "critical",
    "debug",
    "exit_program",
    "info",
    "log",
    "log_file",
    "setup",
    "to_file",
    "warning",
]

_name = "Logger"


class LogLevel(Enum):
    """Enumeration for log levels."""

    CRITICAL = 0
    WARNING = 1
    INFO = 2
    OUT_FILES = 3
    DEBUG = 4


class Colors(Enum):
    """ANSI color codes for terminal output."""

    BLUE = "\033[94m"
    YELLOW = "\033[93m"
    GREEN = "\033[92m"
    RED = "\033[91m"
    LIGHT_BLUE = "\033[96m"
    PURPLE = "\033[95m"
    END = "\033[0m"


# Load configuration data
colors = get_log_colors()
log_levels = get_log_levels()
color_prefix = get_color_prefix()

_init_time = time()
_log_level = 2
_color = True
_stream = sys.stderr
_line_format = "{:<20} {:<19}{:<75} {}\n"
_first_line_format = "{:<20} {:<19}{:<75} \n"
_middle_line_format = "{:<22}{:<75} \n"
_last_line_format = "{:<22}{:<75} {}\n"
_line_break = 76
_remote = False
_prefix = color_prefix
_save_dssp = False
_save_ss = False
_save_restraints = False
_progress_bar = True


def setup(
    log_level: int = 2,
    remote: bool = False,
    work_dir: Union[str, Path] = "",
    save_dssp: bool = False,
    save_ss: bool = False,
    save_restraints: bool = False,
    progress_bar: bool = True,
) -> None:
    """
    Initialize the logging system with specified parameters.

    Args:
        log_level: Logging verbosity level (0-4)
        remote: Whether running in remote mode
        work_dir: Working directory for log files
        save_dssp: Whether to save DSSP output
        save_ss: Whether to save secondary structure output
        save_restraints: Whether to save restraints output
        progress_bar: Whether to show progress bar
    """
    global _log_level, _color, _stream, _remote, _line_break, _prefix, _save_dssp, _save_ss, _save_restraints, _progress_bar
    global _line_format, _middle_line_format, _first_line_format, _last_line_format

    _remote = remote
    work_dir = Path(work_dir) if work_dir else Path.cwd()

    if _remote or not sys.stderr.isatty():
        _color = False
        _line_format = "{:<12} {:<10}{:<75} {}\n"
        _first_line_format = "{:<12} {:<10}{:<75} \n"
        _middle_line_format = "{:<22}{:<75} \n"
        _last_line_format = "{:<22}{:<75} {}\n"
        _prefix = log_levels

    if _remote:
        _log_path = work_dir / "CABS.log"
        try:
            _stream = _log_path.open("a+", encoding="utf-8")
            _stream.write("#" * 110 + "\n")
        except OSError:
            try:
                work_dir.mkdir(parents=True, exist_ok=True)
                _stream = _log_path.open("a+", encoding="utf-8")
            except OSError:
                warning(
                    module_name=_name,
                    msg=f"Could not open a log file at {_log_path}. Writing to standard error instead.",
                )
                raise

    if platform.system() == "Windows":
        _color = False
        _prefix = log_levels

    _log_level = log_level
    info(_name, f"Verbosity set to: {log_level} - {log_levels[str(log_level)]}")
    _save_dssp = save_dssp
    _save_ss = save_ss
    _save_restraints = save_restraints
    _progress_bar = progress_bar


def close_log() -> None:
    """Close the log file if it's not stderr."""
    if _stream is not sys.stderr:
        _stream.close()


def log_files() -> bool:
    """
    Check if verbosity is high enough to save extra output (LOG FILE level).

    Returns:
        True if log level is 3 or higher
    """
    return _log_level >= LogLevel.OUT_FILES.value


def output_dssp() -> bool:
    """
    Check if DSSP output flag was set.

    Returns:
        True if DSSP output should be saved
    """
    return _save_dssp


def output_ss() -> bool:
    """
    Check if secondary structure output flag was set.

    Returns:
        True if SS output should be saved
    """
    return _save_ss


def output_restraints() -> bool:
    """
    Check if restraints output flag was set.

    Returns:
        True if restraints output should be saved
    """
    return _save_restraints


def coloring(color_name: str = "light_blue", msg: str = "") -> str:
    """
    Apply color formatting to a message if colors are enabled.

    Args:
        color_name: Name of the color to apply
        msg: Message to colorize

    Returns:
        Colorized message or plain message if colors disabled
    """
    if _color:
        return colors[color_name] + msg + colors["end"]
    return msg


def log(
    module_name: str = "MISC",
    msg: str = "Processing",
    l_level: int = 2,
    out: Optional[TextIO] = None,
) -> None:
    """
    Log a message with specified level and formatting.

    Args:
        module_name: Name of the module issuing the log
        msg: Message to log
        l_level: Log level (0-4)
        out: Output stream to write to
    """
    if out is None:
        out = _stream
    if l_level <= _log_level:
        t = gmtime(time() - _init_time)
        if len(msg) < _line_break:
            msg_str = _line_format.format(
                _prefix[str(l_level)],
                coloring(msg=f"{module_name}:", color_name="light_blue"),
                msg,
                strftime("(%H:%M:%S)", t),
            )
            out.write(msg_str)
            out.flush()
        else:
            if isinstance(msg, bytes):
                msg = msg.decode("utf-8")
            lines = textwrap.wrap(msg, width=_line_break - 1)
            first_line = _first_line_format.format(
                _prefix[str(l_level)],
                coloring(msg=f"{module_name}:", color_name="light_blue"),
                lines[0],
            )
            out.write(first_line)
            for line_number in range(1, len(lines) - 1):
                line = _middle_line_format.format(" ", lines[line_number])
                out.write(line)
            final_line = _last_line_format.format(
                " ", lines[-1], strftime("(%H:%M:%S)", t)
            )
            out.write(final_line)
            out.flush()


def critical(module_name: str = "_name", msg: str = "") -> None:
    """Log a critical message."""
    log(module_name=module_name, msg=msg, l_level=LogLevel.CRITICAL.value)


def warning(module_name: str = "_name", msg: str = "") -> None:
    """Log a warning message."""
    log(module_name=module_name, msg=msg, l_level=LogLevel.WARNING.value)


def info(module_name: str = "_name", msg: str = "") -> None:
    """Log an info message."""
    log(module_name=module_name, msg=msg, l_level=LogLevel.INFO.value)


def log_file(module_name: str = "_name", msg: str = "") -> None:
    """Log a file-related message."""
    log(module_name=module_name, msg=msg, l_level=LogLevel.OUT_FILES.value)


def debug(module_name: str = "_name", msg: str = "") -> None:
    """Log a debug message."""
    log(module_name=module_name, msg=msg, l_level=LogLevel.DEBUG.value)


def to_file(
    filename: Union[str, Path] = "",
    content: str = "",
    msg: str = "",
    allow_err: bool = True,
    traceback: bool = True,
) -> None:
    """
    Write content to a file with error handling and logging.

    Args:
        filename: Path for the file to be saved
        content: String content to be saved
        msg: Optional message to be logged
        allow_err: If True, log warning on error; if False, exit program
        traceback: If True, raise exception on exit call
    """
    if filename:
        file_path = Path(filename)
        try:
            if file_path.exists():
                log_file(module_name=_name, msg=f"Overwriting {file_path}")
            with file_path.open("w", encoding="utf-8") as f:
                f.write(content if content else "")
        except OSError:
            if allow_err:
                warning(module_name=_name, msg=f"IOError while writing to: {file_path}")
            else:
                exit_program(
                    module_name=_name,
                    msg=f"IOError while writing to: {file_path}",
                    traceback=traceback,
                )
    if msg:
        log_file(module_name=_name, msg=msg)


def exit_program(
    module_name: str = _name,
    msg: str = "Shutting down",
    traceback: Optional[str] = None,
    exc: Optional[Exception] = None,
) -> None:
    """
    Exit the program with appropriate logging and error handling.

    Args:
        module_name: Name of the calling module
        msg: Message to be printed when the program exits
        traceback: Traceback string to print if log level is high enough
        exc: Specific exception passed by the caller
    """
    if exc:
        _msg = f"{msg}: {exc}"
    else:
        _msg = msg
    critical(module_name=module_name, msg=_msg)
    if _log_level > LogLevel.OUT_FILES.value and traceback:
        _stream.write(traceback)
    sys.exit(1)


class ProgressBar:
    """Progress bar for long-running operations with modern Python features."""

    WIDTH: int = 65
    FORMAT: str = "{:<20} {:<19}[{}] {:.1f}%\r"
    BAR0: str = " "
    BAR1: str = "#"

    def __init__(
        self,
        total: int = 100,
        module_name: str = "",
        job_name: str = "",
        out: TextIO = stderr,
        delay: float = 0,
        start_msg: str = "",
    ) -> None:
        """
        Initialize progress bar.

        Args:
            total: Total number of items to process
            module_name: Name of the module using the progress bar
            job_name: Name of the job being tracked
            out: Output stream for progress display
            delay: Initial delay before starting
            start_msg: Message to display at start
        """
        if _log_level >= LogLevel.INFO.value and not _remote and _progress_bar:
            self.stream = out
        else:
            self.stream = open(os.devnull, "w")
        self.total = total
        self.current = 0.0
        self.job_name = job_name
        self.is_done = False
        self.module_name = module_name
        self.prefix = _prefix[str(LogLevel.INFO.value)]

        if start_msg:
            self.stream.write(coloring(msg=start_msg) + "\n")
        if self.job_name:
            log(
                module_name=self.module_name,
                msg=f"{self.job_name} running...",
                out=self.stream,
            )
        self.start_time = time()
        self.update()
        sleep(delay)

    def write(self) -> None:
        """Write the current progress bar state."""
        percent = 1.0 * self.current / self.total
        num = int(self.WIDTH * percent)
        percent = round(100.0 * percent, 1)
        bar = self.BAR1 * num + self.BAR0 * (self.WIDTH - num)
        self.stream.write(
            self.FORMAT.format(
                self.prefix, coloring(msg=f"{self.module_name}:"), bar, percent
            )
        )
        self.stream.flush()

    def update(self, state: float = -1.0) -> bool:
        """
        Update progress bar state.

        Args:
            state: New state value, or -1 to increment by 1

        Returns:
            True if progress is complete
        """
        if state < 0:
            self.current += 1.0
        else:
            self.current = state
        if self.current >= self.total:
            return True
        self.write()
        return False

    def finish(self) -> None:
        """Complete the progress bar to 100%."""
        if self.current < self.total:
            for _ in range(int(self.total - self.current)):
                self.update()
                sleep(0.001)

    def done(self, show_time: bool = True) -> None:
        """
        Mark progress as complete and display completion message.

        Args:
            show_time: Whether to show elapsed time
        """
        if not self.is_done:
            self.finish()
            self.stream.write(" " * 80 + "\r")
            if show_time:
                t = gmtime(time() - self.start_time)
                log(
                    module_name=self.module_name,
                    msg=f"{self.job_name} done in {strftime('%H:%M:%S', t)}",
                    out=self.stream,
                )
            self.stream.flush()
            self.is_done = True


class CabsObserver(Thread):
    """Thread-based observer for monitoring CABS simulation progress."""

    def __init__(
        self,
        interval: float = 0.5,
        progress_file: Optional[Union[str, Path]] = None,
        job_name: str = "CABS simulation",
        msg: str = "",
    ) -> None:
        """
        Initialize CABS observer thread.

        Args:
            interval: Update interval in seconds
            progress_file: Path to file containing progress information
            job_name: Name of the job being monitored
            msg: Initial message to display
        """
        super().__init__()
        self.exit_event = Event()
        self.interval = interval
        self.progress_bar = ProgressBar(
            module_name="CABS", job_name=job_name, start_msg=msg
        )
        self.progress_file = Path(progress_file) if progress_file else None
        self.daemon = True  # In case main program ends abruptly
        self.start()

    def exit(self) -> None:
        """Signal the observer to exit and complete the progress bar."""
        self.progress_bar.done()
        self.exit_event.set()

    def run(self) -> None:
        """Main thread loop for monitoring progress."""
        while not self.exit_event.is_set():
            if self.progress_bar.update(self.status()):
                self.exit()
            sleep(self.interval)

    def status(self) -> float:
        """
        Read current status from progress file.

        Returns:
            Current progress value, or 0 if unable to read
        """
        try:
            if self.progress_file and self.progress_file.exists():
                with self.progress_file.open("rb") as f:
                    progress = float(f.read())
            else:
                progress = 0.0
        except (OSError, ValueError):
            progress = 0.0
        return progress
