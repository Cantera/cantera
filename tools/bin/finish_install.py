import sys, os, string
prefix = '/usr/local' 
pycmd = '/usr/bin/python' # sys.argv[2]
localinst = 0

icyg = prefix.find('/cygdrive/')
if icyg == 0:
    prefix = prefix[10]+':'+prefix[11:]

bindir = '/usr/local/bin'
libdir = '/usr/local/lib/cantera-1.5.3' 
hdrdir = '/usr/local/include/cantera'
demodir = '/usr/local/share/cantera/1.5.3/demos'
datadir = '/usr/local/share/cantera/1.5.3/data'
ctdir = '/usr/local/share/cantera/1.5.3'

f = open(ctdir+'/setup_cantera','w')
f.write('#!/bin/sh\n')
f.write('LD_LIBRARY_PATH='+libdir+':$LD_LIBRARY_PATH\nexport LD_LIBRARY_PATH\n')
f.write('PATH='+bindir+':$PATH\nexport PATH\n')
f.write('PYTHON_CMD='+pycmd+'\nexport PYTHON_CMD\n')
if pycmd <> 'python':
    f.write('alias ctpython='+pycmd+'\n')

ctloc = '-'
warn = ''
warn2 = ''
if localinst:
    try:
        v = sys.version_info
        ctloc = prefix+'/lib/python'+`v[0]`+'.'+`v[1]`+'/site-packages'
        try:
            import Cantera
            ctpath = Cantera.__path__[0]
            if ctpath <> ctloc:
                warn = """
 ######################################################################
    Warning: the Cantera Python package is already installed at
    """+ctpath+""". The newly-installed package at
    """+ctloc+"""/Cantera
    cannot be accessed until the existing one is removed.
 ######################################################################

"""
        except:
            pass
                
        sys.path.append(ctloc)
        f.write('PYTHONPATH='+ctloc+':$PYTHONPATH\nexport PYTHONPATH\n')
    except:
        print 'error'
f.close()

# write the script to run MixMaster
f = open(bindir+'/mixmaster','w')
f.write('#!/bin/sh\n'+pycmd+"""w -c 'from MixMaster import MixMaster; MixMaster()
""")
f.close()

try:
    import Cantera
    ctpath = Cantera.__path__[0]
except:
    ctpath = "-"

fm = open(prefix+"/cantera/ctpath.m","w")
fm.write("""path(path,'"""+prefix+"""/matlab/toolbox/cantera/cantera')\n""")
fm.write("""path(path,'"""+prefix+"""/matlab/toolbox/cantera/cantera/1D')\n""")
fm.close()

print """

Cantera has been successfully installed.

File locations:

    applications      """+bindir+"""
    library files     """+libdir+"""
    C++ headers       """+hdrdir+"""
    demos             """+demodir+"""
    data files        """+datadir+"""
    
    Matlab toolbox    """+prefix+"""/matlab/toolbox/cantera/cantera
    Matlab demos      """+prefix+"""/matlab/toolbox/cantera/cantera-demos
    Matlab tutorials  """+prefix+"""/matlab/toolbox/cantera/cantera-tutorials

    An m-file to set the correct matlab path for Cantera
    is at """+prefix+"""/cantera/ctpath.m"""
    
if ctpath <> "-":
    print """
    Python package    """+ctpath
    if warn <> '':
        print warn
else:
    print """
 ######################################################################
    Warning: the Cantera Python package is not installed. If you
    intentionally skipped it, ignore this message. Otherwise, type
    'make python' and/or 'make python-install' and look for error messages.
    Note that you must first install the 'Numeric' package before installing
    the Cantera package.
 ######################################################################            
"""

print """

    setup script      """+prefix+"""/cantera/setup_cantera
    
    The setup script configures the environment for Cantera. It is
    recommended that you run this script by typing

      source """+prefix+"""/cantera/setup_cantera
    
    before using Cantera, or else
    include its contents in your shell login script.
    """

    
