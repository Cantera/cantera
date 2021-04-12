import os
from os.path import join as pjoin
from os.path import normpath
import sys
import platform
import textwrap
import re
import subprocess
import difflib
import time
import types
import shutil
import itertools
import enum

import SCons.Errors
import SCons
import SCons.Node.FS
from pkg_resources import parse_version
import distutils.sysconfig


from pathlib import Path
import logging

from typing import Union

build_logger = logging.getLogger("build")
build_logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("{levelname!s}: {message!s}", style="{"))
build_logger.addHandler(handler)

output_logger = logging.getLogger("output")
output_logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("{message!s}", style="{"))
output_logger.addHandler(handler)

try:
    import numpy as np
except ImportError:
    np = None


class TestResult(enum.IntEnum):
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

    def __call__(self, target, source, env):
        for s, t in zip(source, target):
            config_h_in = Path(str(s)).read_text()
            config_h = Path(str(t))

            config_h.write_text(config_h_in.format_map(self.defines))
            self.print_config(str(t))

    def print_config(self, filename):
        message = [f"Generating {filename!s} with the following settings:"]

        for key, val in sorted(self.defines.data.items()):
            if val is not None:
                message.append(f"    {key!s:<35} {val}")
        for key in sorted(self.defines.undefined):
            message.append(f"    {key!s:<35} *undefined*")

        build_logger.info("\n".join(message))


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
        values = {
            "passed": sum(self.passed.values()),
            "failed": sum(self.failed.values()),
            "skipped": len(self.tests),
        }
        message = textwrap.dedent("""\
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
            output_logger.error("One or more tests failed.\n" + message)
            sys.exit(1)
        else:
            output_logger.info(message)


test_results = TestResults()


def regression_test(target, source, env):
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
        build_logger.error(f"FAILED (program exit code:{ret.returncode})")
    dir.joinpath(output_name).write_text(ret.stdout)

    diff = 0
    # Compare output files
    comparisons = env["test_comparisons"]
    if blessed_name is not None:
        comparisons.append((Path(blessed_name), output_name))

    for blessed, output in comparisons:
        output_logger.info(f"Comparing '{blessed}' with '{output}'")
        d = compare_files(env, dir.joinpath(blessed), dir.joinpath(output))
        if d:
            output_logger.error("FAILED")
        diff |= d

    for blessed, output in env["test_profiles"]:
        output_logger.info(f"Comparing '{blessed}' with '{output}'")
        d = compare_profiles(env, dir.joinpath(blessed), dir.joinpath(output))
        if d:
            output_logger.error("FAILED")
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
        output_logger.info("PASSED")
        passed_file.write_text(time.asctime())
        test_results.passed[env["active_test_name"]] = 1


def compare_files(env, file1: Path, file2: Path):
    """
    Compare the contents of two files, using a method chosen based on
    their file extensions.
    """
    if file1.suffix == ".csv" and file2.suffix == ".csv":
        return compare_csv_files(env, file1, file2)
    else:
        return compare_text_files(env, file1, file2)


def compare_text_files(env, file1: Path, file2: Path):
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
                        f"Values differ: {num1:14g} {num2:14g}; "
                        f"rel. err = {relerr:.3e}; abs. err = {abserr:.3e}"
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
        output_logger.error("\n".join(message))
        return TestResult.FAIL

    return TestResult.PASS


def get_precision(number: str) -> int:
    """Return the integer representing the power of 10 of the least significant digit in the number represented as a string."""
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


def compare_profiles(env, ref_file, sample_file):
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
        build_logger.warning("Skipping profile comparison because numpy is not available")
        return TestResult.PASS

    atol = env["test_csv_threshold"]
    rtol = env["test_csv_tolerance"]
    xtol = env["test_csv_tolerance"]

    reference = np.genfromtxt(ref_file, delimiter=",").T
    sample = np.genfromtxt(sample_file, delimiter=",").T
    if reference.shape[0] != sample.shape[0]:
        build_logger.error(
            "The output array does not have the same number of variabls as the "
            "reference array."
        )
        return TestResult.FAIL

    # trim header columns if present
    if np.isnan(sample[0, 0]) and np.isnan(reference[0, 0]):
        reference = reference[:, 1:]
        sample = sample[:, 1:]
    if np.isnan(reference).any() or np.isnan(sample).any():
        build_logger.error(
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
            bad.append((reference[0, j], i, reference[i,j], comp[i,j], a, r, x))

    # Fix line lengths here
    footer = []
    maxrows = 10
    if len(bad) > maxrows:
        bad.sort(key=lambda row: -row[5])
        footer += [f"Plus {len(bad) - maxrows} more points exceeding error thresholds."]
        bad = bad[:maxrows]

    if bad:
        output_logger.error("\n".join(header + [template.format(*row) for row in bad] + footer))
        return TestResult.FAIL
    else:
        return TestResult.PASS


def compare_csv_files(env, file1: Path, file2: Path):
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
        build_logger.warning("Skipping profile comparison because numpy is not available")
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
        build_logger.error(f"Could not read data files: {file1}; {file2}", exc_info=e)
        return TestResult.FAIL

    threshold = env["test_csv_threshold"]
    try:
        denom = np.maximum(np.abs(data2), np.abs(data1)) + threshold
        relerror = np.abs(data2 - data1) / denom
        maxerror = np.nanmax(relerror.flat)
    except (ValueError, TypeError) as e:
        build_logger.error("Could not compute error.", exc_info=e)
        return TestResult.FAIL

    tol = env["test_csv_tolerance"]
    if maxerror < tol: # Threshold based on printing 6 digits in the CSV file
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
    build_logger.error("\n".join(message))
    return TestResult.FAIL


def regression_test_message(target, source, env):
    """
    Determines the message printed by SCons while building a
    RegressionTest target.
    """
    return """* Running test '%s'...""" % env['active_test_name']


def add_RegressionTest(env):
    """
    Add "RegressionTest" as a Builder in the specified Scons Environment.
    """
    env['BUILDERS']['RegressionTest'] = env.Builder(
        action=env.Action(regression_test, regression_test_message))


def quoted(s):
    """ Returns the given string wrapped in double quotes."""
    return '"%s"' % s


def mglob(env, subdir, *args):
    """
    Each arg in args is assumed to be file extension,
    unless the arg starts with a '^', in which case the remainder
    of the arg is taken to be a complete pattern.
    """
    matches = []
    for ext in args:
        if ext.startswith('^'):
            matches += env.Glob(pjoin(subdir, ext[1:]))
        else:
            matches += env.Glob(pjoin(subdir, '*.%s' % ext))
    return matches


def psplit(s):
    """
    Split a path given as a string into a list.
    This is the inverse of os.path.join.
    """
    s = s.strip('/\\')
    head, tail = os.path.split(s)
    path = [tail]
    while head:
        head, tail = os.path.split(head)
        path.append(tail)

    path.reverse()
    return path


def subdirs(path):
    """ Get the subdirectories of a specified directory """
    for subdir in os.listdir(path):
        dirpath = pjoin(path, subdir)
        if os.path.isdir(dirpath):
            yield subdir

def stripDrive(s):
    """
    Remove a Windows drive letter specification from a path.
    """
    if len(s) > 1 and s[1] == ':':
        return s[2:]
    else:
        return s


def which(program):
    """ Replicates the functionality of the 'which' shell command """
    def is_exe(fpath):
        for ext in ('', '.exe', '.bat'):
            if os.path.exists(fpath + ext):
                return os.access(fpath + ext, os.X_OK)
        return False

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

optionWrapper = textwrap.TextWrapper(initial_indent='    ',
                                   subsequent_indent='    ',
                                   width=72)


def formatOption(env, opt):
    """
    Print a nicely formatted description of a SCons configuration
    option, its permitted values, default value, and current value
    if different from the default.
    """
    # Extract the help description from the permitted values. Original format
    # is in the format: "Help text ( value1|value2 )" or "Help text"
    if opt.help.endswith(')'):
        parts = opt.help.split('(')
        help = '('.join(parts[:-1])
        values = parts[-1][:-1].strip().replace('|', ' | ')
        if values == '':
            values = 'string'
    else:
        help = opt.help
        values = 'string'

    # Fix the representation of boolean options, which are stored as
    # Python bools, but need to be passed by the user as strings
    default = opt.default
    if default is True:
        default = 'yes'
    elif default is False:
        default = 'no'

    # First line: "* option-name: [ choice1 | choice2 ]"
    lines = ['* %s: [ %s ]' % (opt.key, values)]

    # Help text, wrapped and idented 4 spaces
    lines.extend(optionWrapper.wrap(re.sub(r'\s+', ' ',help)))

    # Default value
    lines.append('    - default: %r' % default)

    # Get the actual value in the current environment
    if opt.key in env:
        actual = env.subst('${%s}' % opt.key)
    else:
        actual = None

    # Fix the representation of boolean options
    if actual == 'True':
        actual = 'yes'
    elif actual == 'False':
        actual = 'no'

    # Print the value if it differs from the default
    if actual != default:
        lines.append('    - actual: %r' % actual)
    lines.append('')

    return lines


def listify(value):
    """
    Convert an option specified as a string to a list, using spaces as
    delimiters. Passes lists transparently.
    """
    if isinstance(value, (list, tuple)):
        # Already a sequence. Return as a list
        return list(value)
    else:
        # assume `value` is a string
        return value.split()

def remove_file(name: Union[Path, str]) -> None:
    """Remove file (if it exists) and print a log message."""
    path_name = Path(name)
    if path_name.exists():
        build_logger.info(f"Removing file '{name!s}'")
        path_name.unlink()

def remove_directory(name: Union[Path, str]) -> None:
    """Remove directory recursively and print a log message."""
    path_name = Path(name)
    if path_name.exists() and path_name.is_dir():
        build_logger.info(f"Removing directory '{name!s}'")
        shutil.rmtree(path_name)

def ipdb():
    """
    Break execution and drop into an IPython debug shell at the point
    where this function is called.
    """
    from IPython.core.debugger import Pdb
    from IPython.core import ipapi

    ip = ipapi.get()
    def_colors = ip.colors
    Pdb(def_colors).set_trace(sys._getframe().f_back)


def getSpawn(env):
    """
    A replacement for env['SPAWN'] on Windows that can deal with very long
    commands, namely those generated when linking. This is only used when
    compiling with MinGW, as SCons automatically uses a tempfile for the
    MSVC link command.

    Pass the return value of this function as the SPAWN keyword argument to
    the Library target, e.g.:

        env.SharedLibrary(..., SPAWN=getSpawn(env))

    Adapted from http://www.scons.org/wiki/LongCmdLinesOnWin32
    """

    if 'cmd.exe' not in env['SHELL'] or env.subst('$CXX') == 'cl':
        return env['SPAWN']

    try:
        useShowWindow = subprocess.STARTF_USESHOWWINDOW
    except AttributeError:
        useShowWindow = subprocess._subprocess.STARTF_USESHOWWINDOW

    def ourSpawn(sh, escape, cmd, args, environ):
        newargs = ' '.join(args[1:])
        cmdline = cmd + " " + newargs
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= useShowWindow
        proc = subprocess.Popen(cmdline,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                startupinfo=startupinfo,
                                shell=False,
                                env=environ)
        data, err = proc.communicate()
        rv = proc.wait()
        if rv:
            print("=====")
            print(err)
            print("=====")
        return rv

    return ourSpawn


def getCommandOutput(cmd, *args):
    """
    Run a command with arguments and return its output.
    """
    environ = dict(os.environ)
    if 'PYTHONHOME' in environ:
        # Can cause problems when trying to run a different Python interpreter
        del environ['PYTHONHOME']
    data = subprocess.check_output([cmd] + list(args), env=environ)
    if sys.version_info.major == 3:
        return data.strip().decode('utf-8')
    else:
        return data.strip()

# Monkey patch for SCons Cygwin bug
# See http://scons.tigris.org/issues/show_bug.cgi?id=2664
if 'cygwin' in platform.system().lower():
    SCons.Node.FS._my_normcase = lambda x: x
