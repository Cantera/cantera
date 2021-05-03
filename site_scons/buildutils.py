import os
import sys
import platform
import textwrap
import re
import subprocess
import difflib
import time
import shutil
import enum
from pathlib import Path
import logging
from typing import TYPE_CHECKING
import collections

try:
    import numpy as np
except ImportError:
    np = None

__all__ = ("logger", "remove_directory", "remove_file", "test_results",
           "add_RegressionTest", "get_command_output", "listify", "which",
           "ConfigBuilder", "multi_glob", "get_spawn", "help", "quoted")

if TYPE_CHECKING:
    from typing import Iterable, Union, List, TypeVar, Dict
    import SCons.Environment
    import SCons.Node.FS
    import SCons.Variables

    SCEnvironment = SCons.Environment.Environment
    LFSNode = List[Union[SCons.Node.FS.File, SCons.Node.FS.Dir]]
    SCVariables = SCons.Variables.Variables
    TPathLike = TypeVar("TPathLike", Path, str)


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
            if isinstance(self.args, collections.Mapping):
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
        values = {
            "passed": sum(self.passed.values()),
            "failed": sum(self.failed.values()),
            "skipped": len(self.tests),
        }
        message = textwrap.dedent(
            """
            *****************************
            ***    Testing Summary    ***
            *****************************

            Tests passed: {passed!s}
            Up-to-date tests skipped: {skipped!s}
            Tests failed: {failed!s}
            """
        ).format_map(values)
        if self.failed:
            message = (message + "Failed tests:" +
                       "".join("\n    - " + n for n in self.failed) +
                       "\n")
        message = message + "*****************************"
        if self.failed:
            logger.error("One or more tests failed.\n" + message, print_level=False)
            sys.exit(1)
        else:
            logger.info(message, print_level=False)


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

    # command line options
    clopts = env["test_command_options"].split()

    # Run the test program
    dir = Path(target[0].dir.abspath)
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
    comparisons = env["test_comparisons"]
    if blessed_name is not None:
        comparisons.append((Path(blessed_name), output_name))

    for blessed, output in comparisons:
        logger.info(f"Comparing '{blessed}' with '{output}'", print_level=False)
        d = compare_files(env, dir.joinpath(blessed), dir.joinpath(output))
        if d:
            logger.error("FAILED", print_level=False)
        diff |= d

    for blessed, output in env["test_profiles"]:
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
    ref_ptp = reference.ptp(axis=1)
    ref_max = np.abs(reference).max(axis=1)
    sample_ptp = sample.ptp(axis=1)
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


def which(program: str) -> bool:
    """Replicates the functionality of the 'which' shell command."""
    for ext in ("", ".exe", ".bat"):
        fpath = Path(program + ext)
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = Path(path).joinpath(fpath)
            if exe_file.exists() and os.access(exe_file, os.X_OK):
                return True
    return False


def help(env: "SCEnvironment", options: "SCVariables") -> None:
    """Print help about configuration options and exit.

    Print a nicely formatted description of a SCons configuration
    option, its permitted values, default value, and current value
    if different from the default.
    """

    message = [
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

                **************************************************"""
        )
    ]

    option_wrapper = textwrap.TextWrapper(
        initial_indent=4 * " ",
        subsequent_indent=4 * " ",
        width=72,
    )
    for opt in options.options:
        # Extract the help description from the permitted values. Original format
        # is in the format: "Help text (value1|value2)" for EnumVariable and
        # BoolVariable types or "Help text" for other Variables
        if opt.help.endswith(")"):
            help, values = opt.help.rsplit("(", maxsplit=1)
            values = values.rstrip(")").strip().replace("|", " | ")
            if not values:
                values = "string"
        else:
            help = opt.help
            values = "string"

        # First line: "* option-name: [ choice1 | choice2 ]"
        lines = [f"* {opt.key}: [ {values} ]"]

        # Help text, wrapped and indented 4 spaces
        lines.extend(option_wrapper.wrap(re.sub(r"\s+", " ", help)))

        # Fix the representation of Boolean options, which are stored as
        # Python bool types
        default = opt.default
        if default is True:
            default = "yes"
        elif default is False:
            default = "no"

        lines.append(f"    - default: {default!r}")

        # Get the actual value in the current environment
        if opt.key in env:
            actual = env.subst(f"${opt.key!s}")
        else:
            actual = None

        # Fix the representation of Boolean options to match the default values
        if actual == "True":
            actual = "yes"
        elif actual == "False":
            actual = "no"

        # Print the value if it differs from the default
        if actual != default:
            lines.append(f"    - actual: {actual!r}")
        message.append("\n".join(lines))

    logger.info("\n\n".join(message), print_level=False)


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

    def our_spawn(sh: str, escape: str, cmd: str, args: str, environ: Dict[str, str]):
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


def get_command_output(cmd: str, *args: str):
    """
    Run a command with arguments and return its output.
    """
    environ = dict(os.environ)
    if "PYTHONHOME" in environ:
        # Can cause problems when trying to run a different Python interpreter
        del environ["PYTHONHOME"]
    data = subprocess.run(
        [cmd] + list(args),
        env=environ,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )
    return data.stdout.strip()


# Monkey patch for SCons Cygwin bug
# See https://github.com/SCons/scons/issues/2664
if "cygwin" in platform.system().lower():
    import SCons.Node.FS
    SCons.Node.FS._my_normcase = lambda x: x
