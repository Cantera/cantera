from glob import glob
import os
import shutil
import sys
from os.path import join as pjoin
import textwrap
import re
import subprocess
import difflib
import time
import types

class DefineDict(object):
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


def regression_test(target, source, env):
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

    if diff or code:
        print 'FAILED'

        if os.path.exists(target[0].abspath):
            os.path.unlink(target[0].abspath)
        return -1
    else:
        print 'PASSED'
        open(target[0].path, 'w').write(time.asctime()+'\n')


def compareFiles(env, file1, file2):
    if file1.endswith('csv') and file2.endswith('csv'):
        return compareCsvFiles(env, file1, file2)
    else:
        return compareTextFiles(env, file1, file2)


def compareTextFiles(env, file1, file2):
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
    except IOError as e:
        print e
        return 1

    relerror = np.abs(data2-data1) / (np.maximum(np.abs(data2), np.abs(data1)) + env['test_csv_threshold'])
    maxerror = np.nanmax(relerror.flat)
    tol = env['test_csv_tolerance']
    if maxerror > tol: # Threshold based on printing 6 digits in the CSV file
        print "Files differ. %i / %i elements above specified tolerance" % (np.sum(relerror > tol), relerror.size)
        return 1
    else:
        return 0


def regression_test_message(target, source, env):
    return """* Running test '%s'...""" % env['active_test_name']


def add_RegressionTest(env):
    env['BUILDERS']['RegressionTest'] = env.Builder(
        action=env.Action(regression_test, regression_test_message))


class CopyNoPrefix(object):
    """
    Copy a file, ignoring leading directories that are part
    of 'prefix' (e.g. the variant directory)
    """
    def __init__(self, prefix):
        self.prefix = prefix

    def __call__(self, source, target, env):
        sourcepath = psplit(str(source[0]))
        targetpath = psplit(str(target[0]))

        depth = 0
        for a,b in zip(targetpath, psplit(self.prefix)):
            if a == b:
                depth += 1
            else:
                break
        print str(source[0]), pjoin(*targetpath[depth:])
        shutil.copyfile(str(source[0]), pjoin(*targetpath[depth:]))


class BuildOpts(object):
    def __init__(self, subdir, name, exts=('cpp',), **kwargs):
        self.subdir = subdir
        self.name = name
        self.extensions = exts
        self.linklibs = kwargs.get('libs', [])


def quoted(s):
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
    head, tail = os.path.split(s)
    path = [tail]
    while head:
        head, tail = os.path.split(head)
        path.append(tail)

    path.reverse()
    return path


def which(program):
    """ Replicates the functionality of the 'which' shell command """
    import os
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
    Convert an option specified as a string to a list.
    Allow both comma and space as delimiters. Passes
    lists transparently.
    """
    if isinstance(value, types.StringTypes):
        return value.replace(',', ' ').split()
    else:
        # Already a sequence. Return as a list
        return list(value)

# This tool adds the builder:
#
#   env.RecursiveInstall(target, path)
#
# This is useful for doing:
#
#   k = env.RecursiveInstall(dir_target, dir_source)
#
# and if any thing in dir_source is updated the install is rerun
#
# It behaves similar to the env.Install builtin. However it expects
# two directories and correctly sets up the dependencies between each
# sub file instead of just between the two directories.
#
# Note in also traverses the in memory node tree for the source
# directory and can detect things that are not built yet. Internally
# we use the env.Glob function for this support.
#
# You can see the effect of this function by doing:
#
#   scons --tree=all,prune
#
# and see the one to one correspondence between source and target
# files within each directory.

def RecursiveInstall(env, target, dir):
    nodes = _recursive_install(env, dir)

    dir = env.Dir(dir).abspath
    target = env.Dir(target).abspath

    l = len(dir) + 1

    relnodes = [n.abspath[l:] for n in nodes]

    out = []
    for n in relnodes:
        t = os.path.join(target, n)
        s = os.path.join(dir, n)
        out.extend(env.InstallAs(env.File(t), env.File(s)))

    return out

def _recursive_install(env, path):
    nodes = env.Glob(os.path.join(path, '*'), strings=False)
    nodes.extend(env.Glob(os.path.join(path, '*.*'), strings=False))
    out = []
    for n in nodes:
        if n.isdir():
            out.extend(_recursive_install(env, n.abspath))
        else:
            out.append(n)

    return out
