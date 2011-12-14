from glob import glob
import os
import shutil
import sys
from os.path import join as pjoin

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
