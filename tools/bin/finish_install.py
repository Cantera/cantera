import sys, os
prefix = sys.argv[1]
pycmd = sys.argv[2]
if prefix == '-': prefix = '/usr/local'

bindir = prefix+'/bin'
libdir = prefix+'/lib/cantera'
hdrdir = prefix+'/include/cantera'

f = open(prefix+'/cantera/cantera_init','w')
f.write('#!/bin/sh\n')
libpath = os.getenv('LD_LIBRARY_PATH')
if libpath:
    f.write('LD_LIBRARY_PATH='+libdir+':'+libpath+'\n')
else:
    f.write('LD_LIBRARY_PATH='+libdir+'\n')    
f.write('PATH='+bindir+':$PATH\n')
f.write('PYTHON_CMD='+pycmd+'\n')
if pycmd <> 'python':
    f.write('alias ctpython='+pycmd+'\n')
f.close()

f = open(bindir+'/mixmaster.py','w')
f.write("""from MixMaster import MixMaster
MixMaster()
""")
f.close()

try:
    import Cantera
    ctpath = Cantera.__path__[0]
except:
    ctpath = "-"
    
print """

Cantera has been successfully installed.

File locations:

    applications      """+bindir+"""
    library files     """+libdir+"""
    C++ headers       """+hdrdir+"""
    demos             """+prefix+"""/cantera/demos
    data files        """+prefix+"""/cantera/data
    
    Matlab toolbox    """+prefix+"""/matlab/toolbox/cantera/cantera
    Matlab demos      """+prefix+"""/matlab/toolbox/cantera/cantera-demos
    Matlab tutorials  """+prefix+"""/matlab/toolbox/cantera/cantera-tutorials"""
if ctpath <> "-":
    print """
    Python package    """+ctpath
print """
    A shell script 'cantera_init' has been written that configures the
    environment for Cantera. It may be found in
    """+prefix+"""/cantera. It is recommended that you run this script
    before using Cantera, or include its contents in your shell login
    script.
    """

    
