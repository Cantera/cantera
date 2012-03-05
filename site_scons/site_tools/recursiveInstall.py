import os
import re

def RecursiveInstall(env, target, dir, exclude=None):
    """
    This tool adds the builder:

        env.RecursiveInstall(target, path)

    This is useful for doing:

        k = env.RecursiveInstall(dir_target, dir_source)

    and if any thing in dir_source is updated the install is rerun

    'exclude' is a list of regular expression patterns for files
    to skip, e.g. ['\\.o$', '^~']

    It behaves similar to the env.Install builtin. However it expects
    two directories and correctly sets up the dependencies between each
    sub file instead of just between the two directories.

    Note in also traverses the in memory node tree for the source
    directory and can detect things that are not built yet. Internally
    we use the env.Glob function for this support.

    You can see the effect of this function by doing:

        scons --tree=all,prune

    and see the one to one correspondence between source and target
    files within each directory.
    """
    if exclude:
        excludePatterns = [re.compile(e) for e in exclude]
    else:
        excludePatterns = []

    nodes = _recursive_install(env, dir, excludePatterns)

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


def _recursive_install(env, path, exclude):
    """ Helper function for RecursiveInstall """
    nodes = env.Glob(os.path.join(path, '*'), strings=False)
    nodes.extend(env.Glob(os.path.join(path, '*.*'), strings=False))
    out = []
    for n in nodes:
        skip = False
        for e in exclude:
            if e.search(n.name):
                skip = True
                break
        if skip:
            continue

        if n.isdir():
            out.extend(_recursive_install(env, n.abspath, exclude))
        else:
            out.append(n)

    return out


def generate(env, **kw):
    env.AddMethod(RecursiveInstall)


def exists(env):
    return True
