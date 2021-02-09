#!/usr/bin/env python3

"""
Collect test coverage data and generate an html report.
"""

import os
import subprocess
import shutil

def getDirectories():
    """
    Return a list of all directories containing coverage data.
    """
    sourcedirs = set()
    rootdir = os.getcwd()
    for dirpath, dirnames, filenames in os.walk(rootdir):
        if 'test_problems' in dirpath or '/ext/' in dirpath:
            continue
        for fname in filenames:
            if fname.endswith('.gcda'):
                dirpath.replace(rootdir, '', 1)
                sourcedirs.add(dirpath)

    return sourcedirs


def clean():
    """
    Remove all coverage data.
    """
    sourcedirs = getDirectories()
    if not sourcedirs:
        return

    dirflags = []
    for d in sourcedirs:
        dirflags.append('-d')
        dirflags.append(d)
    subprocess.call(['lcov', '--zerocounters'] + dirflags)
    subprocess.call(['lcov', '-c', '-i', '-b', '.'] + dirflags + ['-o', 'coverage-base.info'])


def test():
    """
    Run the full test suite.
    """
    subprocess.call(['scons', 'test-reset'])
    subprocess.call(['scons', 'test'])


def collect():
    """
    Collect the generated coverage data into 'coverage.info'
    """
    sourcedirs = getDirectories()
    if not sourcedirs:
        print("Warning! Didn't find any coverage data.")
        return

    command = ['lcov', '-c',
               '-b', '.',
               '-o', 'coverage-raw.info']
    for d in sourcedirs:
        command.append('-d')
        command.append(d)
    subprocess.call(command)

    # Add the baseline
    subprocess.call(['lcov',
                    '-a', 'coverage-raw.info',
                    '-a', 'coverage-base.info',
                    '-o', 'coverage-tmp.info'])

    # Filter to remove non-Cantera code
    subprocess.call(['lcov',
                     '-o', 'coverage.info',
                     '-e', 'coverage-tmp.info',
                     os.getcwd() + '/include/*',
                     os.getcwd() + '/src/*'])
    os.remove('coverage-raw.info')
    os.remove('coverage-base.info')
    os.remove('coverage-tmp.info')

def genhtml():
    """
    Produce an html report from the collected data.
    """
    if os.path.exists('coverage'):
        shutil.rmtree('coverage')

    os.mkdir('coverage')
    subprocess.call(['genhtml', 'coverage.info',
                     '-o', 'coverage',
                     '-p', os.getcwd()])


if __name__ == '__main__':
    clean()
    test()
    collect()
    genhtml()
