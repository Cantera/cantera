#
# SCons builder for gcc's precompiled headers
# Copyright (C) 2006  Tim Blechmann (C) 2011 Pedro Larroy (C) 2015 Carl Cerecke
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.

# 1.3
#
# 09-11-2011 Pedro Larroy: Fixed dependency emitter not working with variant dir
# 07-09-2012 Pedro Larroy: Create the resulting pch on the variant directory,
#   the variant dir should be added to the includes, before any of the other
#   local includes for the compiler to search the pch file on the build directory
#   that we just created. In doubt execute 'strace -f g++ ...' to check that it
#   opens the correct pch.
# 2015-08-07 Carl Cerecke: Hack to work around corruption of .sconsign.dblite
#   caused by SCons.Scanner.C.CScanner() when using both variant dirs, and also
#   separate non-gch static libraries which depend on generated c++ source files.
#   Instead of using SCons.Scanner.C.CScanner, use our own hack.
#
# FIXME: the original precompiled header has to be set
#   in the environment as
#   env['precompiled_header'], this should be fixed
#
#
#  env['precompiled_header'] = File('src/common/includes/all.h')
#  env['Gch'] = env.Gch(target='common/includes/all.h.gch', source=env['precompiled_header'])


from __future__ import print_function
import SCons.Action
import SCons.Builder
import SCons.Scanner.C
import SCons.Util
import SCons.Script
import os
import functools
import re
import subprocess

SCons.Script.EnsureSConsVersion(0,96,92)

GchAction = SCons.Action.Action('$GCHCOM', '$GCHCOMSTR')
GchShAction = SCons.Action.Action('$GCHSHCOM', '$GCHSHCOMSTR')

def gen_suffix(env, sources):
    return sources[0].get_suffix() + env['GCHSUFFIX']

GchShBuilder = SCons.Builder.Builder(action = GchShAction,
                                     source_scanner = SCons.Scanner.C.CScanner(),
                                     suffix = gen_suffix)

GchBuilder = SCons.Builder.Builder(action = GchAction,
                                   source_scanner = SCons.Scanner.C.CScanner(),
                                   suffix = gen_suffix)

def header_path(node):
    h_path = node.abspath
    idx = h_path.rfind('.gch')
    if idx != -1:
        h_path = h_path[0:idx]
        if not os.path.isfile(h_path):
            raise SCons.Errors.StopError("can't find header file: {0}".format(h_path))
        return h_path

    else:
        raise SCons.Errors.StopError("{0} file doesn't have .gch extension".format(h_path))

def dump_node(node):
    print('~Attrs for '+str(node))
    for attr_name in sorted(dir(node)):
        if attr_name[0] == '_':
            continue
        attr = getattr(node, attr_name)
        if type(attr) == list:
            attr = [str(x) for x in attr]
        if not callable(attr):
            print('~',attr_name+":", attr)


inc_re = re.compile(r'#include\s*"([a-zA-Z0-9._]+)"')
def _directly_includes_header(node, header):
    '''
    Simple file search.
    Using SCons.Scanner.C.CScanner with variant_dir can cause corruption of .SConsign.dblite
    and you will get messages during build like:
    "scons: Cannot explain why `output/libGch1.o' is being rebuilt: No previous build information found"
    '''
    for line in open(str(node)):
        if line.startswith('#include'):
            match = inc_re.match(line)
            if match:
                fname = match.group(1)
                if fname == header:
                    return True
            #return False # Should first #include should be precompiled header?
    return False

def pch_emitter(pch_env_key, target, source, env):
    src = source[0].srcnode()

    if not os.path.exists(str(src)):
        return target, source

    if env.get(pch_env_key):
        if _directly_includes_header(src, os.path.basename(env['precompiled_header'].path)):
            if 'explain' in env.GetOption('debug'):
                print('Found dep. on pch: ', source[0], ' -> ', env[pch_env_key])
            env.Depends(target, env[pch_env_key])
    return (target, source)

def generate(env):
    """
    Add builders and construction variables for the Gch builder.
    """
    env.Append(BUILDERS = {
        'gch': env.Builder(
        action = GchAction,
        suffix = 'gch',
        target_factory = env.fs.File,
        ),
        'gchsh': env.Builder(
        action = GchShAction,
        target_factory = env.fs.File,
        ),
        })

    try:
        bld = env['BUILDERS']['Gch']
        bldsh = env['BUILDERS']['GchSh']
    except KeyError:
        bld = GchBuilder
        bldsh = GchShBuilder
        env['BUILDERS']['Gch'] = bld
        env['BUILDERS']['GchSh'] = bldsh

    env['GCHCOM']     = '$CXX -Wall -o $TARGET -x c++-header -c $CXXFLAGS $CCFLAGS $_CCCOMCOM $SOURCE'
    env['GCHSHCOM']   = '$CXX -o $TARGET -x c++-header -c $SHCXXFLAGS $SHCCFLAGS $_CCCOMCOM $SOURCE'
    env['GCHSUFFIX']  = '.gch'


def exists(env):
    return env.Detect('g++')
