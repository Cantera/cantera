from glob import glob
import os
import shutil
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
    """ each arg in args is assumed to be file extension """
    return sum((env.Glob('%s/*.%s' % (subdir, ext)) for ext in args), [])


def psplit(s):
    head, tail = os.path.split(s)
    path = [tail]
    while head:
        head, tail = os.path.split(head)
        path.append(tail)

    path.reverse()
    return path
