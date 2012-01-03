import os
from os.path import join as pjoin
import sys
import platform
import textwrap
import re
import subprocess
import difflib
import time
import types

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
            config_h_in = file(str(s), "r")
            config_h = file(str(t), "w")

            config_h.write(config_h_in.read() % self.defines)
            config_h_in.close()
            config_h.close()
            self.print_config(str(t))

    def print_config(self, filename):
        print 'Generating %s with the following settings:' % filename
        for key, val in sorted(self.defines.data.iteritems()):
            if val is not None:
                print "    %-35s %s" % (key, val)
        for key in sorted(self.defines.undefined):
            print "    %-35s %s" % (key, '*undefined*')


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
        total = len(self.passed) + len(self.failed)
        print """
**********************************
*** Regression Testing Summary ***
**********************************

Tests passed: %(passed)s
Tests failed: %(failed)s
Up-to-date tests skipped: %(skipped)s

**********************************""" % dict(
            passed=len(self.passed),
            failed=len(self.failed),
            skipped=len(self.tests))

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
                               cwd=dir)

    diff = 0
    # Compare output files
    comparisons = env['test_comparisons']
    if blessedName is not None:
        comparisons.append((blessedName,outputName))

    for blessed,output in comparisons:
        print """Comparing '%s' with '%s'""" % (blessed, output)
        diff |= compareFiles(env, pjoin(dir, blessed), pjoin(dir, output))

    del testResults.tests[env['active_test_name']]

    if diff or code:
        print 'FAILED'

        if os.path.exists(target[0].abspath):
            os.path.unlink(target[0].abspath)

        testResults.failed[env['active_test_name']] = 1
    else:
        print 'PASSED'
        open(target[0].path, 'w').write(time.asctime()+'\n')
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
    Compare the contents of two text files, ignoring trailing
    whitespace and any lines starting with strings specified in the
    variable env['test_ignoreLines'].
    """
    text1 = [line.rstrip() for line in open(file1).readlines()
             if not line.startswith(tuple(env['test_ignoreLines']))]
    text2 = [line.rstrip() for line in open(file2).readlines()
             if not line.startswith(tuple(env['test_ignoreLines']))]

    diff = list(difflib.unified_diff(text1, text2))
    if diff:
        'Found differences between %s and %s:' % (file1, file2)
        print '>>>'
        print '\n'.join(diff)
        print '<<<'
        return 1

    return 0


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
        print 'WARNING: skipping .csv diff because numpy is not installed'
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
        data1 = np.genfromtxt(file1, skiprows=headerRows, delimiter=',')
        data2 = np.genfromtxt(file2, skiprows=headerRows, delimiter=',')
    except (IOError, StopIteration) as e:
        print e
        return 1

    relerror = (np.abs(data2-data1) /
                (np.maximum(np.abs(data2), np.abs(data1)) +
                 env['test_csv_threshold']))
    maxerror = np.nanmax(relerror.flat)
    tol = env['test_csv_tolerance']
    if maxerror > tol: # Threshold based on printing 6 digits in the CSV file
        print ("Files differ. %i / %i elements above specified tolerance" %
               (np.sum(relerror > tol), relerror.size))
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
    head, tail = os.path.split(s)
    path = [tail]
    while head:
        head, tail = os.path.split(head)
        path.append(tail)

    path.reverse()
    return path


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
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

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
    option, it's permitted values, default value, and current value
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
    Convert an option specified as a string to a list.  Allow both
    comma and space as delimiters. Passes lists transparently.
    """
    if isinstance(value, types.StringTypes):
        return value.replace(',', ' ').split()
    else:
        # Already a sequence. Return as a list
        return list(value)
