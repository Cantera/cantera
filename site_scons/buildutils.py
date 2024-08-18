from __future__ import annotations
import json
import os
import sys
import textwrap
import re
import subprocess
import difflib
import time
import shutil
import enum
from pathlib import Path
from packaging.version import parse as parse_version
from packaging.specifiers import SpecifierSet
import logging
from typing import TYPE_CHECKING
from collections.abc import Mapping as MappingABC
from SCons.Variables import PathVariable, EnumVariable, BoolVariable
from SCons.Script import Dir

try:
    import numpy as np
except ImportError:
    np = None

__all__ = ("Option", "PathOption", "BoolOption", "EnumOption", "Configuration",
           "logger", "remove_directory", "remove_file", "test_results",
           "add_RegressionTest", "get_command_output", "listify", "which",
           "ConfigBuilder", "multi_glob", "get_spawn", "quoted", "add_system_include",
           "get_pip_install_location", "compiler_flag_list", "setup_python_env",
           "checkout_submodule", "check_for_python", "make_relative_path_absolute",
           "check_sundials", "config_error", "run_preprocessor")

if TYPE_CHECKING:
    from typing import Iterable, TypeVar, Union, List, Dict, Tuple, Optional, \
        Iterator, Sequence
    import SCons.Environment
    import SCons.Node.FS
    import SCons.Variables
    import SCons.SConf

    SCEnvironment = SCons.Environment.Environment
    SConfigure = SCons.SConf.SConfBase
    TextOrSequence = Union[str, List[str], Tuple[str]]
    LFSNode = List[Union[SCons.Node.FS.File, SCons.Node.FS.Dir]]
    SCVariables = SCons.Variables.Variables
    TPathLike = TypeVar("TPathLike", Path, str)


class Option:
    """Object corresponding to SCons configuration option.

    The purpose of this class is to collect system-dependent default parameters for
    SCons configuration options, which are made available to both help text and
    reST documentation. The class works in conjunction with the `Configuration` class
    to select and convert system-dependent parameters to SCons configuration options.
    """

    def __init__(self,
            name: str,
            description: str,
            default: "Union[str, bool, Dict[str, Union[str, bool]]]",
            choices: "Optional[Union[Iterable[str], SCVariables]]" = None):
        self.name = name
        self.description = Option._deblank(description)
        self.default = default
        self.choices = choices

        self._wrapper = textwrap.TextWrapper(width=80)
        self.set_wrapper_indent(4)

    def to_scons(
            self,
            env: "Optional[SCEnvironment]" = None
            ) -> "Union[SCVariables, Tuple[str, str, Union[str, Dict[str, str]]]]":
        """Convert option to SCons variable"""
        default = self.default
        if isinstance(default, str) and "$" in default and env is not None:
            default = env.subst(default)
        elif not isinstance(default, (str, bool)):
            raise TypeError(f"Invalid defaults option with type '{type(default)}'")

        if self.choices is None:
            out = self.name, self.description, default
        else:
            out = self.name, self.description, default, self.choices

        try:
            if isinstance(self, BoolOption):
                return BoolVariable(*out)
            elif isinstance(self, PathOption):
                return PathVariable(*out)
            elif isinstance(self, EnumOption):
                return EnumVariable(*out)
            else:
                return out

        except Exception as err:
            logger.error(
                f"Error converting '{self.name}' to SCons variable:\n{out}")
            raise err

    @property
    def wrapper(self) -> "textwrap.TextWrapper":
        """Line wrapper for text output"""
        return self._wrapper

    def set_wrapper_indent(self, indent: int = 4) -> None:
        """Update indent used for line wrapping"""
        self._wrapper.initial_indent = indent * " "
        self._wrapper.subsequent_indent = indent * " "

    @staticmethod
    def _deblank(string: str) -> str:
        """Remove whitespace before and after line breaks"""
        out = [s.strip() for s in string.split("\n")]
        if not len(out[-1]):
            out = out[:-1]
        out = "\n".join(out)
        return out

    def _build_title(self, backticks: bool = True, indent: int = 3) -> str:
        """Build title describing option and defaults"""
        # First line: "* option-name: [ 'choice1' | 'choice2' ]"

        def decorate(key: str, tick: bool = True) -> str:
            key = f"'{key}'" if tick else key
            return f"``{key}``" if backticks else key

        # format choices
        if isinstance(self, PathOption):
            choices = f"path/to/{self.name}"
            choices = f"{decorate(choices, False)}"
        elif isinstance(self.choices, (list, tuple)):
            choices = list(self.choices)
            for yes_no in ["y", "n"]:
                if yes_no in choices:
                    # ensure consistent order
                    choices.remove(yes_no)
                    choices = [yes_no] + choices
            choices = " | ".join([decorate(c) for c in choices])
        elif isinstance(self, BoolOption) or isinstance(self.default, bool):
            choices = f"{decorate('yes')} | {decorate('no')}"
        else:
            choices = f"{decorate('string', False)}"

        # assemble title
        bullet = f"{'*':<{indent}}"
        return f"{bullet}{decorate(self.name, False)}: [ {choices} ]\n"

    def _build_description(self, backticks: bool = True, indent: int = 3) -> str:
        """Assemble description block (help text)"""

        if not backticks:
            # Help text, wrapped and indented
            self.set_wrapper_indent(indent)
            out = self.wrapper.wrap(re.sub(r"\s+", " ", self.description))
            return "\n".join(out) + "\n"

        # assemble description
        linebreak = "\n" + " " * indent
        description = linebreak.join(self.description.split("\n"))
        pat = r'"([a-zA-Z0-9\-\+\*$_.,: =/\'\\]+)"'
        double_quoted = []
        for item in re.findall(pat, description):
            # enclose double-quoted strings in '``'
            found = f'"{item}"'
            double_quoted += [found]
            replacement = f"``{found}``"
            description = description.replace(found, replacement)
        pat = r"\'([a-zA-Z0-9\-\+\*$_.,:=/\\]+)\'"
        for item in re.findall(pat, description):
            # replace "'" for single-quoted words by '``'; do not replace "'" when
            # whitespace is enclosed or if word is part of double-quoted string
            if any([item in dq for dq in double_quoted]):
                continue
            found = f"'{item}'"
            replacement = found.replace("'", "``")
            description = description.replace(found, replacement)
        pat = r"\*([^\*]+)"
        asterisks = re.findall(pat, description)
        if len(asterisks) == 1:
            # catch unbalanced '*', for example in '*nix'
            found = f"*{asterisks[0]}"
            replacement = fr"\{found}"
            description = description.replace(found, replacement)

        return f"{'':<{indent}}{description}\n"

    @staticmethod
    def _build_defaults(
            defaults: "Union[str, bool, Dict[str, Union[str, bool]]]",
            backticks: bool = True,
            title: str = "default",
            indent: int = 3,
            hanging: int = 0) -> str:
        """Assemble defaults from dictionary"""

        def decorate(key: str, tick: bool = True) -> str:
            key = f"'{key}'" if tick else key
            return f"``{key}``" if backticks else key

        dash = f"{'-':<{indent}}"
        level1 = f"{dash:>{2 * indent + hanging}}"
        level2 = f"{dash:>{3 * indent + hanging}}"

        # Add default platform-specific value
        if isinstance(defaults, bool):
            default = "yes" if defaults else "no"
            return f"{level1}{title}: {decorate(default)}\n"

        if not isinstance(defaults, dict):
            return f"{level1}{title}: {decorate(defaults)}\n"

        # First line of default options
        compilers = {"cl": "MSVC", "gcc": "GCC", "clang": "Clang"}
        toolchains = {"mingw", "intel", "msvc"}
        if any(d in compilers for d in defaults):
            comment = "compiler"
        elif any(d in toolchains for d in defaults):
            comment = "toolchain"
        else:
            comment = "platform"
        out = [f"{level1}{title}: {comment} dependent\n"]

        # Compiler/toolchain options
        if comment in ["compiler", "toolchain"]:
            for key, value in defaults.items():
                if isinstance(value, bool):
                    value = "yes" if value else "no"
                value = decorate(value)
                if key in compilers:
                    key = compilers[key]
                if key == "default":
                    out += [f"{level2}Otherwise: {value}"]
                else:
                    out += [f"{level2}If using {decorate(key, False)}: {value}"]
            return "\n".join(out) + "\n"

        # Platform options
        for key, value in defaults.items():
            if isinstance(value, bool):
                value = "yes" if value else "no"
            value = decorate(value)
            if key == "default":
                out += [f"{level2}Otherwise: {value}"]
            else:
                out += [f"{level2}{key}: {value}"]
        return "\n".join(out) + "\n"

    def to_rest(self, dev: bool = False, indent: int = 3) -> str:
        """Convert description of option to restructured text (reST)"""
        # assemble output
        tag = "sconsopt-" + self.name.replace("_", "-").lower()
        if dev:
            tag += "-dev"
        out = f".. _{tag}:\n\n"
        out += f"{self._build_title(backticks=True)}"
        out += f"{self._build_description(indent=indent)}\n"
        out += f"{self._build_defaults(self.default, backticks=True)}"

        return out

    def help(self, env: "Optional[SCEnvironment]" = None) -> str:
        """Convert option help for command line interface (CLI) output"""
        # assemble output
        out = f"{self._build_title(backticks=False, indent=2)}"
        out += self._build_description(backticks=False, indent=4)
        out += self._build_defaults(self.default, backticks=False, indent=2, hanging=2)

        if env is None:
            return out

        # Get the actual value in the current environment
        actual = env[self.name] if self.name in env else None

        # Fix the representation of Boolean options to match the default values
        if actual in ["True", "False"]:
            actual = True if actual == "True" else False

        # Print the value if it differs from the default
        if actual is not None and str(actual) != str(self.default):
            out += self._build_defaults(
                actual, backticks=False, title="actual", indent=2, hanging=2)

        return out


class PathOption(Option):
    """Object corresponding to SCons PathVariable"""
    pass


class BoolOption(Option):
    """Object corresponding to SCons BoolVariable"""
    pass


class EnumOption(Option):
    """Object corresponding to SCons EnumVariable"""
    pass


class Configuration:
    """
    Class enabling selection of options based on a dictionary of `Option` objects
    that allows for a differentiation between platform/compiler dependent options.

    In addition, the class facilitates the generation of formatted help text and
    reST documentation.
    """

    header = [
        textwrap.dedent(
            """
                **************************************************
                *   Configuration options for building Cantera   *
                **************************************************

        The following options can be passed to SCons to customize the Cantera
        build process. They should be given in the form:

            scons build option1=value1 option2=value2

        Variables set in this way will be stored in the 'cantera.conf' file and reused
        automatically on subsequent invocations of scons. Alternatively, the
        configuration options can be entered directly into 'cantera.conf' before
        running 'scons build'. The format of this file is:

            option1 = 'value1'
            option2 = 'value2'
        """
        )
    ]

    def __init__(self):
        self.options: "Dict[str, Option]" = {}
        self.exported: "List[str]" = []

    def add(self, options: "Union[Option, Sequence[Option]]") -> None:
        """Add new options"""
        if isinstance(options, Option):
            options = [options]
        for item in options:
            if not isinstance(item, Option):
                raise TypeError(f"Invalid option with type '{type(item)}'")
            self.options[item.name] = item

    def list_options(self):
        """Create formatted list of available configuration options"""
        options = sorted(self.options.keys())

        # get formatting information
        n_columns = 3
        n_rows = len(options) // n_columns
        n_extra = len(options) % n_columns
        if n_extra:
            # add one additional row to fit all options
            n_rows += 1
            options.extend([""] * (n_columns - n_extra))

        # get width of individual columns
        max_len = [
            max([len(options[row + col * n_rows]) + 3
                for row in range(n_rows)])
            for col in range(n_columns)]

        # format options as alphabetically ordered columns
        out = []
        for row in range(n_rows):
            line = []
            for col in range(n_columns):
                line.append(f"{options[row + col * n_rows]:{max_len[col]}}")
            out.append("".join(line))
        return "\n".join(out)

    def __getitem__(self, key: str) -> "Option":
        """Make class subscriptable"""
        return self.options[key]

    def __iter__(self) -> "Iterator[str]":
        """Returns itself as an iterator."""
        for k in self.options.keys():
            yield k

    def select(self, key: str = "default") -> None:
        """Select attribute dictionary entries corresponding to *key*"""
        for attr, item in self.options.items():

            if isinstance(item.default, dict) and key in item.default:
                item.default = item.default[key]

    def to_scons(
            self,
            keys: "Union[str, Sequence[str]]" = "",
            env: "Optional[SCEnvironment]" = None
        ) -> "List[SCVariables]":
        """
        Convert options to SCons variables.

        To avoid redundant SCons options, each variable is only exported once.
        """
        if keys == "":
            keys = list(self.options.keys())
        elif isinstance(keys, str):
            keys = [keys]

        keys = [key for key in keys if key not in self.exported]
        out = [self.options[key].to_scons(env) for key in keys]
        self.exported += keys

        return out

    def to_rest(self, option: "Optional[str]" = None, dev: bool = False) -> str:
        """Convert description of configuration options to restructured text (reST)"""
        if option in self.options and option is not None:
            return "\n" + self.options[option].to_rest(dev=dev)
        elif option is not None:
            raise KeyError(f"Unknown option '{option}'")

        message = ["Options List"]
        message.append(f"{'':^<{len(message[-1])}}")
        message.append("")
        for key in self.options.keys():
            message.append(self.options[key].to_rest(dev=dev))

        return "\n".join(message)

    def help(
            self,
            option: "Optional[str]" = None,
            env: "Optional[SCEnvironment]" = None) -> str:
        """Convert configuration help for command line interface (CLI) output"""
        if option in self.options and option is not None:
            return "\n" + self.options[option].help(env=env)
        elif option is not None:
            raise KeyError(f"Unknown option '{option}'")

        message = self.header
        message.append(f"{'':->80}")
        message.append("")

        for value in self.options.values():
            message.append(value.help(env=env))

        return "\n".join(message)


# Add custom logging levels
LOGGING_STATUS_NUM = 32
logging.addLevelName(LOGGING_STATUS_NUM, "STATUS")
LOGGING_FAILED_NUM = 42
logging.addLevelName(LOGGING_FAILED_NUM, "FAILED")


class LevelAdapter(logging.LoggerAdapter):
    """This adapter processes the ``print_level`` keyword-argument to log functions.

    In the default Logger functions, it is not possible to add extra keyword arguments
    to modify behavior. This Adapter allows the Logger functions (.debug(), .info(),
    etc.) to include the keyword argument ``print_level``, which takes a Boolean value
    to determine whether or not the level of the message should be shown with the rest
    of the message. The actual message formatting is handled elsewhere. Example::

       >>> logger.info("Message", print_level=True)
       INFO: Message
       >>> logger.error("Message", print_level=False)
       Message
    """
    def __init__(self, logger):
        self.logger = logger
        self.logger.status = self.status
        self.logger.failed = self.failed

    def status(self, message, *args, **kws):
        # custom logger adapted from https://stackoverflow.com/questions/2183233
        if self.isEnabledFor(LOGGING_STATUS_NUM):
            msg, kwargs = self.process(message, kws)
            # logger takes its '*args' as 'args'
            self._log(LOGGING_STATUS_NUM, msg, args, **kwargs)

    def failed(self, message, *args, **kws):
        # custom logger adapted from https://stackoverflow.com/questions/2183233
        if self.isEnabledFor(LOGGING_FAILED_NUM):
            msg, kwargs = self.process(message, kws)
            # logger takes its '*args' as 'args'
            self._log(LOGGING_FAILED_NUM, msg, args, **kwargs)

    def process(self, msg, kwargs):
        """Pop the value of ``print_level`` into the ``extra`` dictionary.

        Key-value pairs in the "extra" dictionary are set as attributes on the
        ``LogRecord`` instance.
        """
        if "print_level" in kwargs:
            print_level = kwargs.pop("print_level")
            if "extra" in kwargs:
                kwargs["extra"].update(print_level=print_level)
            else:
                kwargs["extra"] = {"print_level": print_level}
        return msg, kwargs


# Modified from https://stackoverflow.com/a/42823461
class BraceLogRecord(logging.LogRecord):
    """Format a log record using brace syntax {} instead of %."""

    def getMessage(self) -> str:
        msg = str(self.msg)
        if self.args:
            if isinstance(self.args, MappingABC):
                msg = msg.format_map(self.args)
            else:
                msg = msg.format(*self.args)
        return msg


class OutputFormatter(logging.Formatter):
    """Format log output depending on whether the level should be shown.

    Intended to be used with the LevelAdapter class to allow the ``print_level``
    keyword argument to be added to Logger method calls (``.info()`` and others). The
    ``print_level`` Boolean value is used to determine whether or not the level of the
    logging message should be printed. By default, the level is shown. Example::

       >>> logger.info("Message", print_level=False)
       Message
       >>> logger.info("Message")
       INFO: Message
    """

    def format(self, record: logging.LogRecord) -> str:
        if record.exc_info or record.exc_text:
            raise ValueError("This formatter does not support exceptions")
        elif record.stack_info:
            raise ValueError("This formatter does not support stack traces")

        no_level_style = "{message}"
        level_style = "{levelname}: " + no_level_style
        record.message = record.getMessage()
        if getattr(record, "print_level", True):
            s = level_style.format(**record.__dict__)
        else:
            s = no_level_style.format(**record.__dict__)
        return s


# Modified from https://stackoverflow.com/a/36338212
class LevelFilter(logging.Filter):
    """Filter out log messages above or below preset cutoffs.

    Log levels in Python correspond to integers, with the lowest, DEBUG, set to
    10 and the highest, CRITICAL, set to 50. This filter causes a log handler to
    reject messages that are above or below numerical cutoffs. Example::

    >>> # Handles log levels from debug up to, but not including, error
    >>> handler.addFilter(LevelFilter(logging.DEBUG, logging.ERROR))
    >>> # Handles log levels from warning up to and including critical
    >>> handler.addFilter(LevelFilter(logging.WARNING, logging.CRITICAL+1))
    """

    def __init__(self, low: int, high: int) -> None:
        self._low = low
        self._high = high
        super().__init__()

    def filter(self, record: logging.LogRecord) -> bool:
        if self._low <= record.levelno < self._high:
            return True
        return False


logging.setLogRecordFactory(BraceLogRecord)
logger = logging.getLogger("cantera")
logger.setLevel(logging.INFO)
f = OutputFormatter()
stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setFormatter(f)
stdout_handler.addFilter(LevelFilter(logging.DEBUG, logging.ERROR))
logger.addHandler(stdout_handler)

stderr_handler = logging.StreamHandler(sys.stderr)
stderr_handler.setFormatter(f)
stderr_handler.addFilter(LevelFilter(logging.ERROR, logging.CRITICAL + 1))
logger.addHandler(stderr_handler)

logger = LevelAdapter(logger)


class TestResult(enum.IntEnum):
    """Represent the passing/failing result of a test.

    To be used instead of a bare integer for clarity.
    """

    PASS = 0
    FAIL = 1


class DefineDict:
    """
    A dictionary-like object which generates appropriate preprocessor
    define statements from its dict of variable / value
    pairs. Variables whose value is None or that are not in the dict
    are left undefined.
    """

    def __init__(self, data: dict) -> None:
        self.data = data
        self.undefined = set()

    def __getitem__(self, key: str) -> str:
        if key not in self.data or self.data[key] is None:
            self.undefined.add(key)
            return f"/* #undef {key!s} */"
        else:
            return f"#define {key!s} {self.data[key]!s}"


class ConfigBuilder:
    """
    Used along with DefineDict to generate a customized config.h file
    from a config.h.in file using the variables given in 'defines'.
    """

    def __init__(self, defines: dict) -> None:
        self.defines = DefineDict(defines)

    def __call__(self, target: "LFSNode", source: "LFSNode", env):
        """
        Note that all three arguments are required by SCons although only the first
        two are used. All of them must be keyword arguments.
        """
        for s, t in zip(source, target):
            config_h_in = Path(str(s)).read_text()
            config_h = Path(str(t))

            config_h.write_text(config_h_in.format_map(self.defines))
            self.print_config(str(t))

    def print_config(self, filename: str):
        message = [f"Generating {filename!s} with the following settings:"]

        for key, val in sorted(self.defines.data.items()):
            if val is not None:
                message.append(f"    {key!s:<35} {val}")
        for key in sorted(self.defines.undefined):
            message.append(f"    {key!s:<35} *undefined*")

        logger.info("\n".join(message))


class TestResults:
    """
    A class that stores information about all the regression tests
    that are defined and which ones have passed / failed in order to
    print a summary at the end of the build process.
    """

    def __init__(self):
        self.tests = {}
        self.passed = {}
        self.failed = {}

    def print_report(self, target, source, env):
        """Print the test results report.

        Note that the three arguments are not used here but are required by SCons,
        and they must be keyword arguments.
        """
        message = [textwrap.dedent(
            f"""
            {' Testing Summary ':*^88}

            Tests passed: {sum(self.passed.values())!s}
            Up-to-date tests skipped: {len(self.tests)!s}
            Tests failed: {sum(self.failed.values())!s}
            """)]
        if self.failed:
            message.append("Failed tests:")
            for failed in self.failed:
                message.append(f"    - {failed}")
            message.append("")
        message.append(f"{'*' * 88}\n")
        message = "\n".join(message)
        logger.status(message, print_level=False)
        if self.failed:
            logger.failed("One or more tests failed.")
            sys.exit(1)


test_results = TestResults()


def regression_test(target: "LFSNode", source: "LFSNode", env: "SCEnvironment"):
    """
    Run a regression test comparing the output of a test program with
    existing "blessed" output files.

    target - The name of an output file that will be generated if
    the test is successful.

    source - A list containing the name of the program to run and
    (optionally) a list of file names to be passed as command line
    arguments.

    The test also relies on several parameters that are passed in via
    variables in the SCons environment:

    env['test_command_options'] - non-file command line options
    to be passed to to the test program
    """
    # unpack:
    program = source[0]
    if len(source) > 1:
        clargs = [s.name for s in source[1:]]
    else:
        clargs = []

    # Name to use for the output file
    blessed_name = env["test_blessed_file"]
    if blessed_name is not None and "blessed" in blessed_name:
        output_name = Path(blessed_name.replace("blessed", "output"))
    else:
        output_name = Path("test_output.txt")

    dir = Path(target[0].dir.abspath)

    comparisons = env["test_comparisons"]
    if blessed_name is not None:
        comparisons.append((Path(blessed_name), output_name))

    # Remove pre-existing output files to prevent tests passing based on outdated
    # output files
    for _, output in comparisons:
        outpath = dir.joinpath(output)
        if outpath.exists():
            outpath.unlink()

    # command line options
    clopts = env["test_command_options"].split()

    # Run the test program
    ret = subprocess.run(
        [program.abspath] + clopts + clargs,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        cwd=dir, env=env["ENV"], universal_newlines=True,
    )
    if ret.returncode:
        logger.error("FAILED (program exit code:{})", ret.returncode)
    dir.joinpath(output_name).write_text(ret.stdout)

    diff = 0

    # Compare output files
    for blessed, output in comparisons:
        if not dir.joinpath(output).is_file():
            logger.error(f"Output file '{output}' not found", print_level=False)
            logger.error("FAILED", print_level=False)
            diff |= TestResult.FAIL
            continue
        logger.info(f"Comparing '{blessed}' with '{output}'", print_level=False)
        d = compare_files(env, dir.joinpath(blessed), dir.joinpath(output))
        if d:
            logger.error("FAILED", print_level=False)
        diff |= d

    for blessed, output in env["test_profiles"]:
        if not dir.joinpath(output).is_file():
            logger.error(f"Output file '{output}' not found", print_level=False)
            logger.error("FAILED", print_level=False)
            diff |= TestResult.FAIL
            continue
        logger.info(f"Comparing '{blessed}' with '{output}'", print_level=False)
        d = compare_profiles(env, dir.joinpath(blessed), dir.joinpath(output))
        if d:
            logger.error("FAILED", print_level=False)
        diff |= d

    del test_results.tests[env["active_test_name"]]

    passed_file = Path(target[0].abspath)
    if diff or ret.returncode:
        if passed_file.exists():
            passed_file.unlink()

        test_results.failed[env["active_test_name"]] = 1
        if env["fast_fail_tests"]:
            sys.exit(1)
    else:
        logger.info("PASSED", print_level=False)
        passed_file.write_text(time.asctime())
        test_results.passed[env["active_test_name"]] = 1


def compare_files(env: "SCEnvironment", file1: Path, file2: Path) -> TestResult:
    """
    Compare the contents of two files, using a method chosen based on
    their file extensions.
    """
    if file1.suffix == ".csv" and file2.suffix == ".csv":
        return compare_csv_files(env, file1, file2)
    else:
        return compare_text_files(env, file1, file2)


def compare_text_files(env: "SCEnvironment", file1: Path, file2: Path) -> TestResult:
    """
    Compare the contents of two text files while:
       - ignoring trailing whitespace
       - ignoring any lines starting with strings specified in the
         variable env['test_ignoreLines'].
       - comparing floating point numbers only up to the printed precision
    """
    text1 = [line.rstrip() for line in file1.read_text().split("\n")
             if not line.startswith(tuple(env["test_ignoreLines"]))]
    text2 = [line.rstrip() for line in file2.read_text().split("\n")
             if not line.startswith(tuple(env["test_ignoreLines"]))]

    # Try to compare the files without testing the floating point numbers
    diff = list(difflib.unified_diff(text1, text2))
    if not diff:
        return TestResult.PASS

    atol = env["test_csv_threshold"]
    rtol = env["test_csv_tolerance"]

    # Replace nearly-equal floating point numbers with exactly equivalent
    # representations to avoid confusing difflib
    float_regex = re.compile(r"(\s*)([+-]{0,1}\d+\.{0,1}\d*([eE][+-]{0,1}\d*){0,1})")
    for i, (line1, line2) in enumerate(zip(text1, text2)):
        if line1 == line2:
            continue

        # group(1) is the left space padding
        # group(2) is the number
        floats1 = [(m.group(1), m.group(2)) for m in float_regex.finditer(line1)]
        floats2 = [(m.group(1), m.group(2)) for m in float_regex.finditer(line2)]

        # If the lines don't contain the same number of numbers,
        # we're not going to pass the diff comparison no matter what
        if len(floats1) != len(floats2):
            continue

        # if the lines don't have the same non-numeric text,
        # we're not going to pass the diff comparison
        if float_regex.sub("", line1).strip() != float_regex.sub("", line2).strip():
            continue

        all_match = True
        for float_1, float_2 in zip(floats1, floats2):
            if float_1 == float_2:
                # String representations match, so replacement is unnecessary
                continue

            try:
                num1 = float(float_1[1])
                num2 = float(float_2[1])
            except ValueError:
                # Something went wrong -- one of the strings isn't actually a number,
                # so just ignore this line and let the test fail
                pass
            else:
                precision = max(get_precision(float_1[1]), get_precision(float_2[1]))
                atol = atol + pow(10, precision) * 1.1
                abserr = abs(num1 - num2)
                relerr = abserr / (0.5 * abs(num1 + num2) + atol)
                if abserr > atol and relerr > rtol:
                    logger.error(
                        "Values differ: {:14g} {:14g}; "
                        "rel. err = {:.3e}; abs. err = {:.3e}",
                        num1, num2, relerr, abserr, print_level=False,
                    )
                    all_match = False
                    break

        # All the values are sufficiently close, so replace the string
        # so that the diff of this line will succeed
        if all_match:
            text2[i] = line1

    # Try the comparison again
    diff = list(difflib.unified_diff(text1, text2))
    if diff:
        message = [f"Found differences between {file1!s} and {file2!s}:", ">>>"]
        message.extend(diff)
        message.append("<<<")
        logger.error("\n".join(message), print_level=False)
        return TestResult.FAIL

    return TestResult.PASS


def get_precision(number: str) -> int:
    """Return the precision of the least significant digit in a number.

    Return an integer representing the power of 10 of the least significant digit in
    ``number``, which must be a string.

    Patterns to consider:
    123 -> 0
    123.45 -> -2
    123.45e6 -> 4
    123e4 -> 4
    """

    number = number.lower()
    if "e" in number:
        number, exponent = number.split("e")
        exponent = int(exponent)
    else:
        exponent = 0

    if "." in number:
        digits = -len(number.split(".")[1])
    else:
        digits = 0

    return exponent + digits


def compare_profiles(
    env: "SCEnvironment", ref_file: Path, sample_file: Path
) -> TestResult:
    """
    Compare two 2D arrays of spatial or time profiles. Each data set should
    contain the time or space coordinate in the first column and data series
    to be compared in successive columns.

    The coordinates in each data set do not need to be the same: The data from
    the second data set will be interpolated onto the coordinates in the first
    data set before being compared. This means that the range of the "sample"
    data set should be at least as long as the "reference" data set.

    After interpolation, each data point must satisfy a combined relative and absolute
    error criterion specified by `rtol` and `atol`.
    """
    if not np:
        logger.warning("Skipping profile comparison because numpy is not available")
        return TestResult.PASS

    atol = env["test_csv_threshold"]
    rtol = env["test_csv_tolerance"]
    xtol = env["test_csv_tolerance"]

    reference = np.genfromtxt(ref_file, delimiter=",").T
    sample = np.genfromtxt(sample_file, delimiter=",").T
    if reference.shape[0] != sample.shape[0]:
        logger.error(
            "The output array does not have the same number of variabls as the "
            "reference array."
        )
        return TestResult.FAIL

    # trim header columns if present
    if np.isnan(sample[0, 0]) and np.isnan(reference[0, 0]):
        reference = reference[:, 1:]
        sample = sample[:, 1:]
    if np.isnan(reference).any() or np.isnan(sample).any():
        logger.error(
            "The output array and reference array have different headers "
            "or contain non-numeric data."
        )
        return TestResult.FAIL

    n_vars = reference.shape[0]
    n_times = reference.shape[1]

    bad = []
    template = "{0:10.4e}  {1:5d}  {2:14.7e}  {3:14.7e}  {4:9.3e}  {5:9.3e}  {6:9.3e}"
    header = ["Failed series comparisons:"]
    header.append("{:10s}  {:5s}  {:14s}  {:14s}  {:9s}  {:9s}  {:9s}".format(
        "coordinate", "comp.", "reference val.", "test value", "abs. err", "rel. err",
        "pos. err"
    ))
    header.append(f"{10*'-'}  -----  {14*'-'}  {14*'-'}  {9*'-'}  {9*'-'}  {9*'-'}")
    ref_ptp = np.ptp(reference, axis=1)
    ref_max = np.abs(reference).max(axis=1)
    sample_ptp = np.ptp(sample, axis=1)
    sample_max = np.abs(sample).max(axis=1)
    scale = np.maximum(
        np.maximum(ref_ptp[1:], ref_max[1:]),
        np.maximum(sample_ptp[1:], sample_max[1:])
    ).reshape(n_vars - 1, -1)
    ref_diff = np.diff(reference)
    slope = ref_diff[1:, :] / ref_diff[0, :] * ref_ptp[0]
    slope = np.hstack((np.zeros((n_vars - 1, 1)), slope))
    comp = np.zeros((n_vars - 1, n_times))
    for i, row in enumerate(sample[1:]):
        comp[i, :] = np.interp(reference[0, :], sample[0, :], row)

    abserr = np.abs(reference[1:] - comp)
    relerr = abserr / (scale + atol)
    # error that can be accounted for by shifting the profile along
    # the time / spatial coordinate
    xerr = abserr / (np.abs(slope) + atol)
    if np.any(abserr > atol) and np.any(relerr > rtol) and np.any(xerr > xtol):
        it = np.nditer((abserr, relerr, xerr), flags=["multi_index"])
        for a, r, x in it:
            i, j = it.multi_index
            bad.append((reference[0, j], i, reference[i, j], comp[i, j], a, r, x))

    footer = []
    maxrows = 10
    if len(bad) > maxrows:
        bad.sort(key=lambda row: -row[5])
        footer += [f"Plus {len(bad) - maxrows} more points exceeding error thresholds."]
        bad = bad[:maxrows]

    if bad:
        logger.error(
            "\n".join(header + [template.format(*row) for row in bad] + footer),
            print_level=False,
        )
        return TestResult.FAIL
    else:
        return TestResult.PASS


def compare_csv_files(env: "SCEnvironment", file1: Path, file2: Path) -> TestResult:
    """
    Compare the contents of two .csv file to see if they are
    similar. Similarity is defined according to tolerances stored in
    the environment as:

        env['test_csv_threshold']
        env['test_csv_tolerance']

    The comparison for each variable is:

        |a-b|/(max(|a|,|b|) + threshold) < tolerance

    Returns 0 if all variables in the files are similar and 1 if the
    files are dissimilar. Lines containing non-numeric data are
    automatically ignored.
    """
    if not np:
        logger.warning("Skipping profile comparison because numpy is not available")
        return TestResult.PASS

    # decide how many header lines to skip
    f = Path(file1).read_text().split("\n")
    header_rows = 0
    good_chars = set("0123456789.+-eE, ")
    for line in f:
        if not set(line).issubset(good_chars):
            header_rows += 1
        else:
            break

    try:
        data1 = np.genfromtxt(file1, skip_header=header_rows, delimiter=",")
        data2 = np.genfromtxt(file2, skip_header=header_rows, delimiter=",")
    except (IOError, StopIteration) as e:
        logger.error(f"Could not read data files: {file1}; {file2}", exc_info=e)
        return TestResult.FAIL

    threshold = env["test_csv_threshold"]
    try:
        denom = np.maximum(np.abs(data2), np.abs(data1)) + threshold
        relerror = np.abs(data2 - data1) / denom
        maxerror = np.nanmax(relerror.flat)
    except (ValueError, TypeError) as e:
        logger.error("Could not compute error.", exc_info=e)
        return TestResult.FAIL

    tol = env["test_csv_tolerance"]
    if maxerror < tol:  # Threshold based on printing 6 digits in the CSV file
        return TestResult.PASS

    n_fail = np.sum(relerror > tol)
    n_tot = relerror.size
    message = [
        "Files differ.",
        f"{n_fail:d} / {n_tot:d} elements above specified tolerance ({tol:f})",
        "  row   col   reference       test            rel. error",
        "  ----  ----  --------------  --------------  ----------",
    ]
    it = np.nditer([relerror, data1, data2], flags=["multi_index"])
    for rele, ref, test in it:
        if rele > tol:
            r = it.multi_index[0] + header_rows + 1
            c = it.multi_index[1] + 1
            message.append(
                f"  {r:4d}  {c:4d}  {ref:14.7e}  {test:14.7e}  {rele:10.4e}"
            )
    logger.error("\n".join(message))
    return TestResult.FAIL


def regression_test_message(target, source, env: "SCEnvironment") -> str:
    """
    Determines the message printed by SCons while building a
    RegressionTest target.

    Note that the first two arguments are not used but are required by SCons and they
    must be keyword arguments.
    """
    return f"""* Running test '{env["active_test_name"]}'..."""


def add_RegressionTest(env: "SCEnvironment") -> None:
    """
    Add "RegressionTest" as a Builder in the specified Scons Environment.
    """
    env["BUILDERS"]["RegressionTest"] = env.Builder(
        action=env.Action(regression_test, regression_test_message)
    )


def compiler_flag_list(
        flags: "Union[str, Iterable]",
        compiler: str,
        excludes: "Optional[Iterable]" = (),
    ) -> "List[str]":
    """
    Separate concatenated compiler flags in ``flags``.

    ``compiler`` is either ``"cl"`` for MSVC or anything else for a different
    compiler.

    Entries starting with the regular expression patterns in ``excludes`` are omitted.
    """
    if not isinstance(flags, str):
        flags = " ".join(flags)

    if compiler == "cl":
        # Options can start with "/", or "$"
        expr = r"""(?:^|\ +)           # start of string or leading whitespace
                   ([/\$].+?)          # capture start of option
                   (?=\ +[-/\$]|\ *$)  # start of next option or end of string
                """
    else:
        # Options can start with "-"
        expr = r"""(?:^|\ +)           # start of string or leading whitespace
                   (-.+?)              # capture start of option
                   (?=\ +-|\ *$)  # start of next option or end of string
                """

    # split concatenated entries
    flags = re.findall(expr, flags, re.VERBOSE)

    # Remove duplicates and excluded items
    excludes = tuple(excludes)
    cc_flags = []
    for flag in flags:
        if flag in cc_flags:
                continue
        if not any(re.match(exclude, flag) for exclude in excludes):
            cc_flags.append(flag)

    return cc_flags


def add_system_include(env, include, mode='append'):
    # Add a file to the include path as a "system" include directory, which will
    # suppress warnings stemming from code that isn't part of Cantera, and reduces
    # time spent scanning these files for changes.
    if mode == 'append':
        add = env.Append
    elif mode == 'prepend':
        add = env.Prepend
    else:
        raise ValueError("mode must be 'append' or 'prepend'")

    if env['CC'] == 'cl':
        add(CPPPATH=include)
    else:
        if isinstance(include, (list, tuple)):
            for inc in include:
                add(CXXFLAGS=('-isystem', inc))
        else:
            add(CXXFLAGS=('-isystem', include))


def quoted(s: str) -> str:
    """Return the given string wrapped in double quotes."""
    return f'"{s}"'


def multi_glob(env: "SCEnvironment", subdir: str, *args: str):
    """Use SCons Glob to find nodes in a subdirectory using many file extensions.

    Each argument in ``args`` is assumed to be file extension,
    unless the arg starts with a ``'^'``, in which case the remainder
    of the argument is taken to be a complete pattern.
    """
    matches = []
    for ext in args:
        if ext.startswith("^"):
            matches += env.Glob(Path(subdir).joinpath(ext[1:]))
        else:
            matches += env.Glob(Path(subdir).joinpath(f"*.{ext}"))
    return matches


def which(program: str) -> "Optional[str]":
    """Replicates the functionality of the 'which' shell command."""
    for ext in ("", ".exe", ".bat"):
        fpath = Path(program + ext)
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = Path(path).joinpath(fpath)
            if exe_file.exists() and os.access(exe_file, os.X_OK):
                return str(exe_file)
    return None


def listify(value: "Union[str, Iterable]") -> "List[str]":
    """
    Convert an option specified as a string to a list, using spaces as
    delimiters. Passes lists and tuples transparently.
    """
    if isinstance(value, str):
        return value.split()
    else:
        # Already a sequence. Return as a list
        return list(value)


def remove_file(name: "TPathLike") -> None:
    """Remove file (if it exists) and print a log message."""
    path_name = Path(name)
    if path_name.exists():
        logger.info(f"Removing file '{name!s}'")
        path_name.unlink()


def remove_directory(name: "TPathLike") -> None:
    """Remove directory recursively and print a log message."""
    path_name = Path(name)
    if path_name.exists() and path_name.is_dir():
        logger.info(f"Removing directory '{name!s}'")
        shutil.rmtree(path_name)


def ipdb():
    """
    Break execution and drop into an IPython debug shell at the point
    where this function is called.
    """
    from IPython.core.debugger import Pdb
    from IPython.core import getipython

    ip = getipython.get_ipython()
    def_colors = ip.colors
    Pdb(def_colors).set_trace(sys._getframe().f_back)


def get_spawn(env: "SCEnvironment"):
    """
    A replacement for env['SPAWN'] on Windows that can deal with very long
    commands, namely those generated when linking. This is only used when
    compiling with MinGW, as SCons automatically uses a tempfile for the
    MSVC link command.

    Pass the return value of this function as the SPAWN keyword argument to
    the Library target, for example:

        env.SharedLibrary(..., SPAWN=get_spawn(env))

    Adapted from https://github.com/SCons/scons/wiki/LongCmdLinesOnWin32
    """

    if "cmd.exe" not in env["SHELL"] or env.subst("$CXX") == "cl":
        return env["SPAWN"]

    def our_spawn(sh: str, escape: str, cmd: str, args: str, environ: "Dict[str, str]"):
        newargs = " ".join(args[1:])
        cmdline = cmd + " " + newargs
        startupinfo = subprocess.STARTUPINFO()  # type: ignore
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW  # type: ignore
        proc = subprocess.Popen(cmdline,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                startupinfo=startupinfo,
                                shell=False,
                                env=environ)
        _, err = proc.communicate()
        rv = proc.wait()
        if rv:
            logger.error(err)
        return rv

    return our_spawn


def get_command_output(cmd: str, *args: str, ignore_errors=False):
    """
    Run a command with arguments and return its output.
    """
    environ = dict(os.environ)
    if "PYTHONHOME" in environ:
        # Can cause problems when trying to run a different Python interpreter
        del environ["PYTHONHOME"]
    kwargs = {}
    if ignore_errors:
        kwargs["stderr"] = subprocess.DEVNULL
    data = subprocess.run(
        [cmd] + list(args),
        env=environ,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        check=not ignore_errors,
        **kwargs,
    )
    return data.stdout.strip()

_python_info = None
def setup_python_env(env):
    """Set up an environment for compiling Python extension modules"""

    global _python_info
    if _python_info is None:
        # Get information needed to build the Python module
        script = textwrap.dedent("""\
        from sysconfig import *
        import numpy
        import json
        import site
        import sys
        vars = get_config_vars()
        vars["plat"] = get_platform()
        vars["numpy_include"] = numpy.get_include()
        vars["site_packages"] = [d for d in site.getsitepackages() if d.endswith("-packages")]
        vars["user_site_packages"] = site.getusersitepackages()
        vars["abiflags"] = getattr(sys, "abiflags", "")
        print(json.dumps(vars))
        """)
        _python_info = json.loads(get_command_output(env["python_cmd"], "-c", script))

    info = _python_info
    module_ext = info["EXT_SUFFIX"]
    inc = info["INCLUDEPY"]
    pylib = info.get("LDLIBRARY")
    prefix = info["prefix"]
    py_version_short = parse_version(info["py_version_short"])
    py_version_full = parse_version(info["py_version"])
    py_version_nodot = info["py_version_nodot"]
    plat = info['plat'].replace('-', '_').replace('.', '_')
    numpy_include = info["numpy_include"]
    env.Prepend(CPPPATH=Dir('#include'))
    if env["system_sundials"] == "n":
        env.Prepend(CPPPATH=Dir('#include/cantera/ext'))

    add_system_include(env, (inc, numpy_include), 'prepend')
    env.Prepend(LIBS=env['cantera_shared_libs'])

    # Fix the module extension for Windows from the sysconfig library.
    # See https://github.com/python/cpython/pull/22088 and
    # https://bugs.python.org/issue39825
    if (py_version_full < parse_version("3.8.7")
        and env["OS"] == "Windows"
        and module_ext == ".pyd"
    ):
        module_ext = f".cp{py_version_nodot}-{info['plat'].replace('-', '_')}.pyd"

    env["py_module_ext"] = module_ext
    env["py_version_nodot"] = py_version_nodot
    env["py_version_short"] = info["py_version_short"]
    env["py_plat"] = plat
    env["py_base"] = info["installed_base"]
    env["site_packages"] = info["site_packages"]
    env["user_site_packages"] = info["user_site_packages"]
    if env["OS"] != "Windows":
        env["py_libpath"] = [info["LIBPL"], info["LIBDIR"]]
        py_lib = "python" + info["py_version_short"] + info["abiflags"]
    else:
        env["py_libpath"] = [info["installed_base"] + "\\libs"]
        py_lib = "python" + py_version_nodot
    env["py_libs"] = [py_lib] + [lib[2:] for lib in info.get("LIBS", "").split()
                                 if lib.startswith("-l")]

    # Don't print deprecation warnings for internal Python changes.
    # Only applies to Python 3.8. The field that is deprecated in Python 3.8
    # and causes the warnings to appear will be removed in Python 3.9 so no
    # further warnings should be issued.
    if env["HAS_CLANG"] and py_version_short == parse_version("3.8"):
        env.Append(CXXFLAGS='-Wno-deprecated-declarations')

    if env['OS'] == 'Darwin':
        env.Append(LINKFLAGS='-undefined dynamic_lookup')
    elif env['OS'] == 'Windows':
        env.Append(LIBPATH=prefix + '/libs')
        if env['toolchain'] == 'mingw':
            env.Append(LIBS=f"python{py_version_nodot}")
            if env['OS_BITS'] == 64:
                env.Append(CPPDEFINES='MS_WIN64')

    if not env.get("require_numpy_1_7_API", True):
        env.Append(CPPDEFINES="NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION")


    return env

def get_pip_install_location(
    python_cmd: str,
    user: bool = False,
    prefix: str | None = None,
    root: str | None = None
) -> dict[str, str]:
    """Determine the location where pip will install files.

    This relies on pip's internal API so it may break in future versions.
    Unfortunately, I don't really see another way to determine this information
    reliably.
    """
    # These need to be quoted if they're not None, even if they're a falsey value
    # like the empty string. Otherwise, we want the literal None value.
    prefix = quoted(prefix) if prefix is not None else None
    root = quoted(root) if root is not None else None
    install_script = textwrap.dedent(f"""
        from pip import __version__ as pip_version
        from packaging.version import parse as parse_version
        import pip
        import json
        pip_version = parse_version(pip_version)
        if pip_version < parse_version("10.0.0"):
            from pip.locations import distutils_scheme
            scheme = distutils_scheme("Cantera", user={user}, root={root},
                                      prefix={prefix})
        else:
            from pip._internal.locations import get_scheme
            scheme = get_scheme("Cantera", user={user}, root={root},
                                prefix={prefix})

        if not isinstance(scheme, dict):
            scheme = {{k: getattr(scheme, k) for k in dir(scheme)
                       if not k.startswith("_")}}
        print(json.dumps(scheme))
    """)
    return json.loads(get_command_output(python_cmd, "-c", install_script))


def checkout_submodule(name: str, submodule_path: str):
    if not os.path.exists(".git"):
        logger.error(f"{name} is missing. Extract package in {submodule_path}.")
        sys.exit(1)

    try:
        code = subprocess.call(["git", "submodule", "update", "--init",
                                "--recursive", submodule_path])
    except Exception:
        code = -1
    if code:
        logger.error(f"{name} submodule checkout failed.\n"
                      "Try manually checking out the submodule by running:\n\n"
                     f"    git submodule update --init --recursive {submodule_path}\n")
        sys.exit(1)


def check_for_python(
    env: "SCEnvironment",
    command_line_targets: List[str]
) -> Dict[str, Union[str | bool]]:
    """Check for compatible versions of Python and Python-specific dependencies.

    Args:
        env (SCons.Environment): The SCons construction environment.
        command_line_targets (list[str]): The list of targets passed on the command line
            by the user.

    Returns:
        dict[str, str | bool]: Dictionary with either one or two keys:
            * ``"python_package"``: String with the value ``"y"`` or ``"n"`` to
              indicate whether or not the Python package will be built. This key will
              always be present.
            * ``"require_numpy_1_7_API"``: Boolean indicating whether Cython will use
              NumPy 1.7 API support. This is purely a compile-time constant and does not
              affect which versions of NumPy Cantera supports. This key may or may not
              be present in the return dictionary.

    Raises:
        config_error: If the user specifies that the Python package should be built but
            no compatible versions of Python, Cython, and NumPy can be found.
    """
    # Pytest is required only to test the Python module
    check_for_pytest = "test" in command_line_targets or any(
        target.startswith(("test-python", "test-help")) for target in command_line_targets
    )

    # Check for the minimum ruamel.yaml version at install and test
    # time. The check happens at install and test time because ruamel.yaml is
    # only required to run the Python interface, not to build it.
    check_for_ruamel_yaml = check_for_pytest or any(
        target in command_line_targets
        for target in ["install", "test"]
    )

    def check_module(name):
        return textwrap.dedent(f"""\
            try:
                import {name}
                versions["{name}"] = {name}.__version__
            except ImportError as {name}_err:
                err += str({name}_err) + "\\n"
        """)

    script = textwrap.dedent("""\
        import sys
        import json
        versions = {}
        versions["python"] = "{v.major}.{v.minor}".format(v=sys.version_info)
        err = ""
    """)
    script += check_module("numpy")
    script += check_module("Cython")

    if check_for_ruamel_yaml:
        script += textwrap.dedent("""\
            try:
                from ruamel import yaml
                versions["ruamel.yaml"] = yaml.__version__
            except ImportError as ru_err:
                try:
                    import ruamel_yaml as yaml
                    versions["ruamel.yaml"] = yaml.__version__
                except ImportError as ru_err_2:
                    err += str(ru_err) + "\\n"
                    err += str(ru_err_2) + "\\n"
        """)
    if check_for_pytest:
        script += check_module("pytest")

    script += textwrap.dedent("""\
        print("versions:", json.dumps(versions))
        if err:
            print(err)
    """)

    warn_no_python = False
    try:
        info = get_command_output(env["python_cmd"], "-c", script).splitlines()
    except OSError as err:
        logger.debug(f"Error checking for Python:\n{err}")
        warn_no_python = True
    except subprocess.CalledProcessError as err:
        logger.debug(f"Error checking for Python:\n{err} {err.output}")
        warn_no_python = True

    if warn_no_python:
        if env["python_package"] == "default":
            logger.warning(
                "Not building the Python package because the Python interpreter "
                f"{env['python_cmd']!r} could not be found.")
            return {"python_package": "n"}
        else:
            logger.error(
                f"Could not execute the Python interpreter {env['python_cmd']!r}")
            sys.exit(1)

    versions = {}
    for line in info:
        if line.startswith("versions:"):
            versions = {
                k: parse_version(v)
                for k, v in json.loads(line.split(maxsplit=1)[1]).items()
            }
            break

    if len(info) > 1:
        msg = ["Unexpected output while checking Python dependency versions:"]
        msg.extend(line for line in info if not line.startswith("versions:"))
        logger.warning("\n| ".join(msg))

    if env["python_package"] == "y":
        logger_method = logger.error
        exit_on_error = True
    else:
        logger_method = logger.warning
        exit_on_error = False

    python_version = versions.get("python")
    if not python_version or python_version < env["python_min_version"]:
        logger_method(
            f"Python version is incompatible. Found {python_version} but "
            f"{env['python_min_version']} or newer is required. In order to install "
            "Cantera without Python support, specify 'python_package=n'.")
        if exit_on_error:
            sys.exit(1)
        return {"python_package": "n"}
    elif python_version >= env["python_max_version"]:
        logger.warning(
            f"Python {python_version} is not supported for Cantera "
            f"{env['cantera_version']}. Python versions {env['python_max_version']} and "
            "newer are untested and may result in unexpected behavior. Proceed "
            "with caution.")

    numpy_version = versions.get("numpy")
    if not numpy_version:
        logger_method("NumPy not found. Not building the Python package.")
        if exit_on_error:
            sys.exit(1)
        return {"python_package": "n"}
    elif numpy_version not in env["numpy_version_spec"]:
        logger_method(
            f"NumPy is an incompatible version: Found {numpy_version} but "
            f"{env['numpy_version_spec']} is required.")
        if exit_on_error:
            sys.exit(1)
        return {"python_package": "n"}
    else:
        logger.info(f"Using NumPy version {numpy_version}")

    cython_version = versions.get("Cython")
    if not cython_version:
        logger_method("Cython not found. Not building the Python package.")
        if exit_on_error:
            sys.exit(1)
        return {"python_package": "n"}
    elif cython_version not in env["cython_version_spec"]:
        logger_method(
            f"Cython is an incompatible version: Found {cython_version} but "
            f"{env['cython_version_spec']} is required.")
        if exit_on_error:
            sys.exit(1)
        return {"python_package": "n"}
    elif cython_version < parse_version("3.0.0"):
        logger.info(
            f"Using Cython version {cython_version} (uses legacy NumPy API)"
        )
        require_numpy_1_7_API = True
    else:
        logger.info(f"Using Cython version {cython_version}")
        require_numpy_1_7_API = False

    if check_for_ruamel_yaml:
        ruamel_yaml_version = versions.get("ruamel.yaml")
        if not ruamel_yaml_version:
            logger.error(
                f"ruamel.yaml was not found. {env['ruamel_version_spec']} "
                "is required.")
            sys.exit(1)
        elif ruamel_yaml_version not in env["ruamel_version_spec"]:
            logger.error(
                "ruamel.yaml is an incompatible version: Found "
                f"{ruamel_yaml_version}, but {env['ruamel_version_spec']} "
                "is required.")
            sys.exit(1)
        else:
            logger.info(f"Using ruamel.yaml version {ruamel_yaml_version}")

    if check_for_pytest:
        pytest_version = versions.get("pytest")
        if not pytest_version:
            logger.error(
                f"pytest was not found. {env['pytest_version_spec']} "
                "is required.")
            sys.exit(1)
        elif pytest_version not in env["pytest_version_spec"]:
            logger.error(
                "pytest is an incompatible version: Found "
                f"{pytest_version}, but {env['pytest_version_spec']} "
                "is required.")
            sys.exit(1)
        else:
            logger.info(f"Using pytest version {pytest_version}")

    return {"python_package": "y", "require_numpy_1_7_API": require_numpy_1_7_API}


def make_relative_path_absolute(path_to_check: Union[str, Path]) -> str:
    """If a path is absolute, return it in POSIX format.
    If a path is relative, assume it's relative to the source root, convert it to an
    absolute path, and return the converted path in POSIX format.
    """
    pth = Path(path_to_check)
    if not pth.is_absolute():
        pth = Path(Dir("#" + path_to_check).abspath)

    return pth.as_posix()


def run_preprocessor(
    conf: "SConfigure",
    includes: "TextOrSequence",
    text: "TextOrSequence",
    defines: "TextOrSequence" = ()
) -> Tuple[int, str]:
    """Run the C preprocessor and return the last line of the processed source.

    This function can be used to extract ``#define``ed values from header files to use,
    for example, to determine whether a header is present or the version of an external
    dependency.

    A source file is constructed by concatenating the ``includes``, ``text``, and
    ``defines``. That file is passed to the C preprocessor selected SCons for the
    environment. If the preprocessor errors during execution, this function returns the
    preprocessor return code and an empty string. If preprocessing is successful, the
    return code of the preprocessor and the last line of the preprocessor output are
    returned.

    Args:
        conf (SCons.SConf.SConfBase): An instance of the SConf configuration object.
        includes (str | list[str] | tuple[str]): Names of headers to include in the
            constructed source file.
        text (str | list[str] | tuple[str]): The text of the source file to include in
            the constructed source file.
        defines (str | list[str] | tuple[str]): Names of variables to ``#define`` in
            the constructed source file. Optional, defaults to an empty tuple.

    Returns:
        tuple[int, str]: The return code of the preprocessor and a string of the
            preprocessor output. If the preprocessor encounters an error, the string
            will be empty.
    """
    if not isinstance(includes, (tuple, list)):
        includes = [includes]
    if not isinstance(text, (tuple, list)):
        text = [text]
    if not isinstance(defines, (tuple, list)):
        defines = [defines]
    if "msvc" in conf.env["toolchain"]:
        # Yes this needs 4 slashes on each side. The first pair are escaped by Python,
        # the second pair are escaped by the shell to leave a single backslash to be
        # interpreted by the Windows shell as a directory separator.
        preprocessor_flags = ['/P', '/Fi".\\\\.sconf_temp\\\\"']
    else:
        preprocessor_flags = ["-E"]
    conf.env.Prepend(CXXFLAGS=preprocessor_flags)
    source = ["#define " + d for d in defines]
    source.extend("#include " + ii for ii in includes)
    source.extend(text)
    retcode = conf.TryCompile(text="\n".join(source), extension=".cpp")
    for flag in preprocessor_flags:
        conf.env["CXXFLAGS"].remove(flag)
    if retcode:
        retval = None
        # On MSVC, the `/P` flag produces an output file named for the input file,
        # that is, `conftest_<hash of contents>_0.i`. However, SCons assumes that the
        # `.obj` file will still be the target and sets that as the `conf.lastTarget`
        # with the filename: `conftest_<hash of contents>_0_<hash of action>.obj`.
        # Since we need the contents of the `.i` file, we need this string munging on
        # the `.obj` filename to find the right file. If SCons ever decides to change
        # how they name conftest files, this will probably break.
        if "msvc" in conf.env["toolchain"]:
            fname = conf.lastTarget.name.rsplit("_", maxsplit=1)[0]
            content = Path(".sconf_temp", fname).with_suffix(".i").read_text().splitlines()
        else:
            content = conf.lastTarget.get_text_contents().splitlines()
        # Go from bottom to top of the file, since the preprocessor is going to spit
        # out lots of irrelevant lines.
        for line in reversed(content):
            # The MSVC compiler tends to produce extra output at the end of the `.i`
            # file that we don't want.
            if not line.strip() or line.strip().startswith("conftest"):
                continue
            retval = line.strip()
            break
        if retval is None:
            raise ValueError("Could not find version. See config.log.")
        return retcode, retval
    else:
        return retcode, ""


def check_sundials(conf: "SConfigure", sundials_version: str) -> Dict[str, Union[str, int, bool]]:
    """Check for the version of SUNDIALS and whether SUNDIALS was built with BLAS/LAPACK

    Args:
        conf (SCons.SConf.SConfBase): An instance of the SConf configuration object.
        sundials_version (str): A string with the version of SUNDIALS. The expected
            format is ``"X Y Z"`` including the quote symbols.

    Returns:
        Dict[str, str]: A dictionary with three keys:
            * ``"system_sundials"``: String ``"y"`` or ``"n"`` indicating whether a
              compatible version of SUNDIALS was found in the system directories.
            * ``"sundials_version"``: String with the parsed version of SUNDIALS with
              periods between the version components.
            * ``"sundials_blas_lapack"``: Boolean or integer indicating whether SUNDIALS was built
              with BLAS/LAPACK support. Always ``False`` if ``"system_sundials"`` is
              ``"n"``.

    Raises:
        ``config_error``: If the user specified ``system_sundials==y`` but a compatible
            version could not be found.
    """
    sundials_ver = parse_version(
        ".".join(sundials_version.strip().replace('"', "").split())
    )
    should_exit_with_error = conf.env["system_sundials"] == "y"
    if sundials_ver < parse_version("3.0") or sundials_ver >= parse_version("8.0"):
        if should_exit_with_error:
            config_error(f"Sundials version must be >=3.0,<8.0. Found {sundials_ver}.")
        return {"system_sundials": "n", "sundials_version": "", "has_sundials_lapack": 0}
    elif sundials_ver > parse_version("7.0.0"):
        logger.warning(f"Sundials version {sundials_ver} has not been tested.")

    cvode_checks = {
        SpecifierSet(">=3.0,<4.0"): ("CVodeCreate(CV_BDF, CV_NEWTON);", ["sundials_cvodes"]),
        SpecifierSet(">=4.0,<6.0"): ("CVodeCreate(CV_BDF);", ["sundials_cvodes"]),
        SpecifierSet(">=6.0,<7.0"): (
            "SUNContext ctx; SUNContext_Create(0, &ctx);", ["sundials_cvodes"]
        ),
        SpecifierSet(">=7.0,<8.0"): (
            "SUNContext ctx; SUNContext_Create(SUN_COMM_NULL, &ctx);", ["sundials_core"]
        )
    }
    for version_spec, (cvode_call, libs) in cvode_checks.items():
        if sundials_ver in version_spec:
            ret = conf.CheckLibWithHeader(
                libs,
                header="cvodes/cvodes.h",
                language='C++',
                call=cvode_call,
                autoadd=False,
            )
            # CheckLibWithHeader returns True to indicate success
            if not ret:
                if should_exit_with_error:
                    config_error(
                        "Could not link to the Sundials library. Did you set the "
                        "include/library paths?"
                    )
                return {"system_sundials": "n", "sundials_version": "", "has_sundials_lapack": 0}
            break

    logger.info(f"Using system installation of Sundials version {sundials_ver}.")

    # Determine whether or not Sundials was built with BLAS/LAPACK
    if sundials_ver < parse_version("5.5"):
        # In Sundials 2.6-5.4, SUNDIALS_BLAS_LAPACK is either defined or undefined
        has_sundials_lapack = conf.CheckDeclaration('SUNDIALS_BLAS_LAPACK',
                '#include "sundials/sundials_config.h"', 'C++')
    elif sundials_ver <= parse_version("6.6.0"):
        # In Sundials 5.5-6.6.0, two defines are included specific to the
        # SUNLINSOL packages indicating whether SUNDIALS has been built with LAPACK
        lapackband = conf.CheckDeclaration(
            "SUNDIALS_SUNLINSOL_LAPACKBAND",
            '#include "sundials/sundials_config.h"',
            "C++",
        )
        lapackdense = conf.CheckDeclaration(
            "SUNDIALS_SUNLINSOL_LAPACKDENSE",
            '#include "sundials/sundials_config.h"',
            "C++",
        )
        has_sundials_lapack = lapackband and lapackdense
    else:
        # In Sundials 6.6.1, the SUNDIALS_BLAS_LAPACK_ENABLED macro was introduced
        has_sundials_lapack = conf.CheckDeclaration("SUNDIALS_BLAS_LAPACK_ENABLED",
                '#include "sundials/sundials_config.h"', 'c++')

    if not has_sundials_lapack and conf.env['use_lapack']:
        logger.warning(
            "External BLAS/LAPACK has been specified for Cantera but SUNDIALS was built "
            "without this support. Cantera will use the slower default solver "
            "implementations included with SUNDIALS. You can resolve this warning by "
            "installing or building SUNDIALS with BLAS/LAPACK support."
        )

    return {
        "system_sundials": "y",
        "sundials_version": str(sundials_ver),
        "has_sundials_lapack": has_sundials_lapack
    }


def config_error(message: str) -> None:
    """Log an error message to the console and exit the build with code 1."""
    if logger.getEffectiveLevel() == logging.DEBUG:
        logger.error(message)
        debug_message = [
            f"\n{' Contents of config.log: ':*^80}\n",
            Path("config.log").read_text().strip(),
            f"\n{' End of config.log ':*^80}",
        ]
        logger.debug("\n".join(debug_message), print_level=False)
    else:
        error_message = [message]
        error_message.append("\nSee 'config.log' for details.")
        logger.error("\n".join(error_message))
    sys.exit(1)
