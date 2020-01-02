from __future__ import print_function
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

import SCons.Errors
import SCons
import SCons.Node.FS
from distutils.version import LooseVersion, StrictVersion
import distutils.sysconfig

class DefineDict(object):
    """
    A dictionary-like object which generates appropriate preprocessor
    define statements from its dict of variable / value
    pairs. Variables whose value is None or that are not in the dict
    are left undefined.
    """
    def __init__(self, data):
        self.data = data
        self.undefined = set()

    def __getitem__(self, key):
        if key not in self.data:
            self.undefined.add(key)
            return '/* #undef %s */' % key
        elif self.data[key] is None:
            return '/* #undef %s */' % key
        else:
            return '#define %s %s' % (key, self.data[key])


class ConfigBuilder(object):
    """
    Used along with DefineDict to generate a customized config.h file
    from a config.h.in file using the variables given in 'defines'.
    """
    def __init__(self, defines):
        self.defines = DefineDict(defines)

    def __call__(self, source, target, env):
        for s, t in zip(source, target):
            config_h_in = open(str(s), "r")
            config_h = open(str(t), "w")

            config_h.write(config_h_in.read() % self.defines)
            config_h_in.close()
            config_h.close()
            self.print_config(str(t))

    def print_config(self, filename):
        print('Generating %s with the following settings:' % filename)
        for key, val in sorted(self.defines.data.items()):
            if val is not None:
                print("    %-35s %s" % (key, val))
        for key in sorted(self.defines.undefined):
            print("    %-35s %s" % (key, '*undefined*'))


class TestResults(object):
    """
    A class that stores information about all the regression tests
    that are defined and which ones have passed / failed in order to
    print a summary at the end of the build process.
    """
    def __init__(self):
        self.tests = {}
        self.passed = {}
        self.failed = {}

    def printReport(self, target, source, env):
        if self.failed:
            failures = ('Failed tests:' +
                        ''.join('\n    - ' + n for n in self.failed) +
                        '\n')
        else:
            failures = ''
        print("""
*****************************
***    Testing Summary    ***
*****************************

Tests passed: %(passed)s
Up-to-date tests skipped: %(skipped)s
Tests failed: %(failed)s
%(failures)s
*****************************""" % dict(
            passed=sum(self.passed.values()),
            failed=sum(self.failed.values()),
            skipped=len(self.tests),
            failures=failures))

        if self.failed:
            raise SCons.Errors.BuildError(self, 'One or more tests failed.')


testResults = TestResults()


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
    blessedName = env['test_blessed_file']
    if blessedName is not None and 'blessed' in blessedName:
        outputName = blessedName.replace('blessed', 'output')
    else:
        outputName = 'test_output.txt'

    # command line options
    clopts = env['test_command_options'].split()

    # Run the test program
    dir = str(target[0].dir.abspath)
    with open(pjoin(dir,outputName), 'w') as outfile:
        code = subprocess.call([program.abspath] + clopts + clargs,
                               stdout=outfile, stderr=outfile,
                               cwd=dir, env=env['ENV'])

    if code:
        print('FAILED (program exit code:{0})'.format(code))

    diff = 0
    # Compare output files
    comparisons = env['test_comparisons']
    if blessedName is not None:
        comparisons.append((blessedName,outputName))

    for blessed,output in comparisons:
        print("""Comparing '%s' with '%s'""" % (blessed, output))
        d = compareFiles(env, pjoin(dir, blessed), pjoin(dir, output))
        if d:
            print('FAILED')
        diff |= d

    del testResults.tests[env['active_test_name']]

    if diff or code:
        if os.path.exists(target[0].abspath):
            os.path.unlink(target[0].abspath)

        testResults.failed[env['active_test_name']] = 1
        if env["fast_fail_tests"]:
            sys.exit(1)
    else:
        print('PASSED')
        with open(target[0].path, 'w') as passed_file:
            passed_file.write(time.asctime()+'\n')
        testResults.passed[env['active_test_name']] = 1


def compareFiles(env, file1, file2):
    """
    Compare the contents of two files, using a method chosen based on
    their file extensions.
    """
    if file1.endswith('csv') and file2.endswith('csv'):
        return compareCsvFiles(env, file1, file2)
    else:
        return compareTextFiles(env, file1, file2)


def compareTextFiles(env, file1, file2):
    """
    Compare the contents of two text files while:
       - ignoring trailing whitespace
       - ignoring any lines starting with strings specified in the
         variable env['test_ignoreLines'].
       - comparing floating point numbers only up to the printed precision
    """
    text1 = [line.rstrip() for line in open(file1).readlines()
             if not line.startswith(tuple(env['test_ignoreLines']))]
    text2 = [line.rstrip() for line in open(file2).readlines()
             if not line.startswith(tuple(env['test_ignoreLines']))]

    # Try to compare the files without testing the floating point numbers
    diff = list(difflib.unified_diff(text1, text2))
    if not diff:
        return 0

    atol = env['test_csv_threshold']
    rtol = env['test_csv_tolerance']

    # Replace nearly-equal floating point numbers with exactly equivalent
    # representations to avoid confusing difflib
    reFloat = re.compile(r'(\s*)([+-]{0,1}\d+\.{0,1}\d*([eE][+-]{0,1}\d*){0,1})')
    for i in range(min(len(text1), len(text2))):
        line1 = text1[i]
        line2 = text2[i]
        if line1 == line2:
            continue

        # group(1) is the left space padding
        # group(2) is the number
        floats1 = [(m.group(1),m.group(2))
                   for m in list(reFloat.finditer(line1))]
        floats2 = [(m.group(1),m.group(2))
                   for m in list(reFloat.finditer(line2))]

        # If the lines don't contain the same number of numbers,
        # we're not going to pass the diff comparison no matter what
        if len(floats1) != len(floats2):
            continue

        # if the lines don't have the same non-numeric text,
        # we're not going to pass the diff comparison
        if reFloat.sub('', line1).strip() != reFloat.sub('', line2).strip():
            continue

        allMatch = True
        for j in range(len(floats1)):
            if floats1[j] == floats2[j]:
                # String representations match, so replacement is unnecessary
                continue

            try:
                delta = max(getPrecision(floats1[j][1]), getPrecision(floats2[j][1]))
                num1 = float(floats1[j][1])
                num2 = float(floats2[j][1])
                abserr = abs(num1-num2)
                relerr = abserr / (0.5 * abs(num1 + num2) + atol)
                if abserr > (1.1*delta + atol) and relerr > rtol:
                    print('Values differ: {0: 14g} {1: 14g}; rel. err = {2:.3e}; abs. err = {3:.3e}'.format(num1, num2, relerr, abserr))
                    allMatch = False
                    break
            except Exception as e:
                # Something went wrong -- one of the strings isn't actually a number,
                # so just ignore this line and let the test fail
                pass

        # All the values are sufficiently close, so replace the string
        # so that the diff of this line will succeed
        if allMatch:
            text2[i] = line1

    # Try the comparison again
    diff = list(difflib.unified_diff(text1, text2))
    if diff:
        print('Found differences between %s and %s:' % (file1, file2))
        print('>>>')
        print('\n'.join(diff))
        print('<<<')
        return 1

    return 0


def getPrecision(x):
    """
    Return the number corresponding to the least significant digit of
    the number represented by the string 'x'.
    """
    x = x.lower()
    # Patterns to consider:
    # 123
    # 123.45
    # 123.45e6
    # 123e4
    if 'e' in x:
        x, exponent = x.split('e')
        exponent = int(exponent)
    else:
        exponent = 0

    decimalPt = x.find('.')
    if decimalPt == -1:
        decimalPt = len(x) - 1

    precision = decimalPt + exponent - len(x) + 1

    return 10**precision


def compareCsvFiles(env, file1, file2):
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
    try:
        import numpy as np
    except ImportError:
        print('WARNING: skipping .csv diff because numpy is not installed')
        return 0

    # decide how many header lines to skip
    f = open(file1)
    headerRows = 0
    goodChars = set('0123456789.+-eE, ')
    for line in f:
        if not set(line[:-1]).issubset(goodChars):
            headerRows += 1
        else:
            break

    try:
        data1 = np.genfromtxt(file1, skip_header=headerRows, delimiter=',')
        data2 = np.genfromtxt(file2, skip_header=headerRows, delimiter=',')
    except (IOError, StopIteration) as e:
        print(e)
        return 1

    try:
        relerror = (np.abs(data2-data1) /
                    (np.maximum(np.abs(data2), np.abs(data1)) +
                     env['test_csv_threshold']))
        maxerror = np.nanmax(relerror.flat)
    except ValueError as e:
        print(e)
        return 1

    tol = env['test_csv_tolerance']
    if maxerror > tol: # Threshold based on printing 6 digits in the CSV file
        print("Files differ. %i / %i elements above specified tolerance (%f)" %
              (np.sum(relerror > tol), relerror.size, tol))
        print('  row   col   reference     test          rel. error')
        print('  ----  ----  ------------  ------------  ----------')
        for i,j in itertools.product(*map(range, relerror.shape)):
            if relerror[i,j] > tol:
                row = i + headerRows + 1
                col = j + 1
                print('  % 4i  % 4i  % 12f  % 12f  % 10f' %
                      (row, col, data1[i,j], data2[i,j], relerror[i,j]))
        return 1
    else:
        return 0


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

def removeFile(name):
    """ Remove file (if it exists) and print a log message """
    if os.path.exists(name):
        print('Removing file "%s"' % name)
        os.remove(name)

def removeDirectory(name):
    """ Remove directory recursively and print a log message """
    if os.path.exists(name):
        print('Removing directory "%s"' % name)
        shutil.rmtree(name)

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
