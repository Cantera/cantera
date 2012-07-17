#!/usr/bin/python

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

    command = ['lcov', '--zerocounters']
    for d in sourcedirs:
        command.append('-d')
        command.append(d)
    subprocess.call(command)


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
        print "Warning! Didn't find any coverage data."
        return

    command = ['lcov', '-c',
               '-b', '.',
               '-o', 'coverage.info']
    for d in sourcedirs:
        command.append('-d')
        command.append(d)
    subprocess.call(command)


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
