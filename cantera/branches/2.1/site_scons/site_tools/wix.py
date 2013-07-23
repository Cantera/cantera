"""
Tool to support WiX (Windows Installer XML toolset)
http://blogs.msdn.com/robmen/
http://sourceforge.net/projects/wix
http://www.scons.org/wiki/WiX_Tool
"""
__revision__ = "Revision: 1.1"
__date__ = "Date: 2004/05/21 20:44:46"
__author__ = "elliot.murphy@veritas.com"
__credits__ = ""

import os

import SCons.Defaults
import SCons.Util
import SCons.Scanner

def generate(env):
    """Add Builders and construction variables for WiX to an Environment."""
    if not exists(env):
        return

    env['WIXCANDLE'] = '"%sbin/candle.exe"' % os.environ['WIX']
    env['WIXCANDLEFLAGS'] = ['-nologo']
    env['WIXCANDLEINCLUDE'] = []
    env['WIXCANDLECOM'] = '$WIXCANDLE $WIXCANDLEFLAGS -I $WIXCANDLEINCLUDE -o ${TARGET} ${SOURCE}'

    env['WIXLIGHT'] = '"%sbin/light.exe"' % os.environ['WIX']
    env['WIXLIGHTFLAGS'] = ['-nologo', '-spdb']
    env['WIXLIGHTCOM'] = "$WIXLIGHT $WIXLIGHTFLAGS -out ${TARGET} ${SOURCES}"

    object_builder = SCons.Builder.Builder(
        action = '$WIXCANDLECOM',
        suffix = '.wxiobj',
        src_suffix = '.wxs')

    linker_builder = SCons.Builder.Builder(
        action = '$WIXLIGHTCOM',
        src_suffix = '.wxiobj',
        src_builder = object_builder)

    env['BUILDERS']['WiX'] = linker_builder

def exists(env):
    return 'WIX' in os.environ
