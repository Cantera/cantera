########################################################################
#
#   Cantera installation / upgrade utility.
#
#   Download this file, save it as ctupdate.py, and process it with python
#   to install Cantera or update an existing installation.
#
#########################################################################


"""
    usage:  python ctupdate.py <options>

    valid options are:
    --help          print this message
    --force         update everything, whether it needs it or not
    --src           compile/install from source code
    --bin           do a binary installation
    --[no]python    [don't] install the Cantera Python interface
    --[no]matlab    [don't] install the Cantera Matlab toolbox
    
    This script installs Cantera or updates an existing installation.
    Both Windows and linux/unix installations are supported.  It
    installs from binaries on Windows, and builds everything from the
    source code on all non-Windows platforms.
    
    If you are updating and have more than one Python interpreter
    installed, be sure to run this script using the interpreter
    you use with Cantera. On linux systems, this may be 'python2.'

"""

import urllib
import os
import sys
import time

v = float(sys.version[:3])
if v < 2.0:
    print """

### ERROR ###

This script must be run with python version 2.0 or greater.
This is version """+sys.version+""".
"""
    if sys.platform[:5] == 'linux':
        print """On many linux systems, python 2.x is installed as 'python2'.
"""

    sys.exit(-1)



#
# definitions
#
_CTINFO_FILE = 'http://www.cantera.org/cantera_dist/etc/_ctinfo.py'


#
# set these to zero to skip installation of Python and Matlab components
# 
install_python_module = 1
install_matlab_toolbox = 1

    
#
#  global variables
#
_instdir = '--'    # installation directory
_srcdir = '--'     # directory where downloaded files go
force_update = 0   # force updating if nonzero

_updates = []
_reason = []       # 
_num = 0
_bininstall = -1
_new_python = 0

_options = {}
_local = {}
_srv = {}


def read_options():
    """Read local configuration information."""
    optionfile = get_options_fname()
    if not os.access(optionfile,os.F_OK):
        return
    f = open(optionfile)
    lines = f.readlines()
    for line in lines:
        if line:
            w = line.split()+['null','null']
            if w[1][-1] == '\n':
                v = w[1][:-1]
            else:
                v = w[1]
            _local[w[0]] = v
            _options[w[0]] = v


def get_options_fname():
    
    # if the OS is Windows, store the config file in
    # C:\Documents and Settings\<username>\Application Data
    if sys.platform == 'win32':
        appdata = os.getenv('APPDATA')
        if not appdata:
            print """
            Application Data directory cannot be found, because
            environment variable APPDATA is not set.
            """
            sys.exit(-1)
            
        home = appdata+os.sep+'cantera'
        if not os.access(home,os.F_OK):
            os.makedirs(home)

    # otherwise, look for environment variable 'HOME' to determine
    # the home directory
    else:
        home = os.getenv('HOME')

    # if HOME is not set, then look in the current working directory.
    if not home:
        home = os.getcwd()

    optionfile = home+os.sep+'.ctupdate'
    return optionfile


    
def write_options():
    """
    Save configuration options.
    """
    optionfile = get_options_fname()
    f = open(optionfile,'w')
    lines = []
    for k in _options.keys():
        if k:
            if type(_options[k]) == types.StringType:
                lines.append(k+' '+_options[k]+'\n')
            else:
                lines.append(k+' '+`_options[k]`+'\n')
    f.writelines(lines)
    f.close()

            
def ask(msg, options, deflt, helpmsg=''):
    """
    Prompt for input.
    options -- string listing the possible answers
    deflt   -- default response
    helpmsg -- message printed if 'h' is input.
    """
    print msg,
    if options:
        print '['+options+'] ',
    if deflt:
        print '('+deflt+') ',
    ans = sys.stdin.readline()[:-1]
    if ans == 'q' or ans == 'quit':
        sys.exit(0)
    elif ans == 'h':
        print helpmsg
        ans = ask(msg, options, deflt)
        
    if not ans: return deflt
    else: return ans



def report(prog, blocks, size):
    """Print dots and percent completed during file download."""
    global _num
    sys.stdout.write('.')
    _num = _num + 1
    if 20*(prog/20) == prog: sys.stdout.write('%5.1f' %
                                              ((100.0)*prog*blocks/size,)+'%')


def check_write():
    """
    Check to see if the Python library directory is writable.
    """
    prefix = sys.prefix
    return os.access(prefix+'/lib',os.W_OK)


def pycmd():
    """
    Return the location of the Python interpreter used to run Cantera,
    or if this cannot be found, then return the string 'python'.
    """
    if _options.has_key('pycmd'):
        return _options['pycmd']
    else:
        return 'python'


def get_srcdir():

    """ Enter the directory where the installation files downloaded
    from the web should be put. Once the installation is complete, any
    files in this directory may be deleted. (But if you are building
    from the source code, you might want to keep them, since the
    source code is not copied to the installation directory.)  """
    
    global _srcdir
    if _srcdir <> '--':
        return _srcdir

    if _local.has_key('download_dir'):
        d = _local['download_dir']
    else:
        try:
            d = os.getenv('HOME')
        except:
            d = os.getcwd()
            
    dir = ask('Download directory:','<directory>|q|h',d,
              get_srcdir.__doc__)

    dir = os.path.expanduser(dir)
    _srcdir = os.path.abspath(dir)
    if not os.access(_srcdir,os.F_OK):
        a = ask('Directory '+_srcdir+' does not exist. Create it?','y/n/q/h','y')
        if a == 'y':
            os.makedirs(_srcdir)
    os.chdir(_srcdir)
    _options['download_dir'] = _srcdir
    return _srcdir


def get_instdir():

    """ Enter the root directory under which binary executables,
    library files, etc., should be installed. If you are doing a
    system-wide install and have write-access to the /usr partition,
    just press <Return>. Otherwise, enter the name of an existing
    directory you have write access to. In this case, subdirectories
    bin, lib, cantera, etc. will be created in this directory if they
    don't already exist.  """

    global _instdir
    if _instdir <> '--':
        return _instdir
    if _local.has_key('install_dir'):
        d = _local['install_dir']
    else:
        d = ''
    dir = ask('Installation directory:','<directory>|q|h',d,get_instdir.__doc__)
    ndir = dir
    try:
        if dir <> '' and dir[0] == '~':
            ndir = os.getenv('HOME')+dir[1:]
    except:
        ndir = dir
    _instdir = os.path.abspath(ndir)
    _options['install_dir'] = _instdir
    return _instdir



def get_python():
    """
    Download the version of Python Cantera requires.
    """
    global _new_python
    _new_python = 1
    server = _ctinfo_srv.server
    
    if platform == 'win32':
        path = _ctinfo_srv.WinPythonPath
        urllib.urlretrieve(server+path,'installPython.exe',report)
        _install_script.write('installPython.exe\n')
        _options['pycmd'] = 'python'

    elif _bininstall == 0:
        _prefix = get_instdir()
        path = _ctinfo_srv.SrcPythonPath
        urllib.urlretrieve(server+path,'python.tar.gz',report)
        
        _install_script.write('gunzip python.tar.gz\ntar xvf python.tar\n')
        _install_script.write('cd '+_ctinfo_srv.SrcPythonDir+'\n')
        if _prefix:
            _install_script.write('./configure --prefix='+_prefix
                                  +' --exec-prefix='+_prefix+'\n')
        else:
            _install_script.write('./configure\n')
        _install_script.write('make\nmake install\ncd ..\n\n')
        if _prefix:
            _options['pycmd'] = _prefix+'/bin/python'
        else:
            _options['pycmd'] = 'python'
        _options['python_version'] = _srv['python_version']



def get_graphviz():
    """Download graphviz."""
    server = _ctinfo_srv.server

    _options['dot_version'] = _srv['dot_version']
            
    if platform == 'win32' and _bininstall == 1:
        path = _ctinfo_srv.WinGraphvizPath
        urllib.urlretrieve(server+path,'installGraphviz.exe',report)
        _install_script.write('installGraphviz.exe\n')

    elif _bininstall == 0:
        _prefix = get_instdir()
        path = _ctinfo_srv.SrcGraphvizPath
        urllib.urlretrieve(server+path,'graphviz.tar.gz',report)
        
        _install_script.write('gunzip graphviz.tar.gz\ntar xvf graphviz.tar\n')
        _install_script.write('cd '+_ctinfo_srv.SrcGraphvizDir+'\n')
        if _prefix:
            _install_script.write('./configure --prefix='+_prefix
                                  +' --exec-prefix='+_prefix+'\n')
        else:
            _install_script.write('./configure\n')
        _install_script.write('make\nmake install\ncd ..\n\n')


def get_numeric():
    """Download the version of Numeric Cantera requires."""
    server = _ctinfo_srv.server
    
    if platform == 'win32' and _bininstall == 1:
        path = _ctinfo_srv.WinNumericPath
        urllib.urlretrieve(server+path,'installNumeric.exe',report)
        _install_script.write('installNumeric.exe\n')
    else:
        path = _ctinfo_srv.SrcNumericPath
        urllib.urlretrieve(server+path,'numeric.tar.gz',report)
        
        _install_script.write('gunzip numeric.tar.gz\ntar xvf numeric.tar\n')
        _install_script.write('cd '+_ctinfo_srv.SrcNumericDir+'\n')
        _install_script.write(pycmd()+' setup.py install\ncd ..\n\n')


def get_ct():
    """Download cantera."""
    server = _ctinfo_srv.server
    server = 'http://prdownloads.sourceforge.net/cantera/'
    ctfile = 'cantera-1.5.4.tar.gz'
    mirror = 'easynews'
    url = server+ctfile+'?use_mirror='+mirror
    
    _options['cantera'] = _srv['cantera']
    
    if platform == 'win32' and _bininstall == 1:
        path = _ctinfo_srv.WinCanteraPath
        urllib.urlretrieve(server+path,url,report)
        _install_script.write('Cantera13.msi\n')
        _install_script.write('Cantera13.msi\n')        
        path = _ctinfo_srv.WinCanteraPyPath
        urllib.urlretrieve(server+path,'installCanteraPy.exe',report)
        _install_script.write('installCanteraPy.exe\n')        

    elif _bininstall == 0:
        _prefix = get_instdir()
        path = _ctinfo_srv.SrcCanteraPath
        print urllib.urlopen(url) #,'get.html',report)
        
        _install_script.write('gunzip cantera.tar.gz\ntar xvf cantera-1.3.tar\n')
        _install_script.write('cd '+_ctinfo_srv.SrcCanteraDir+'\n')
        _install_script.write('PYTHON_CMD='+pycmd()+'; export PYTHON_CMD\n')
        if install_matlab_toolbox == 0:
            _install_script.write('BUILD_MATLAB_TOOLBOX="n"; export BUILD_MATLAB_TOOLBOX\n')
        if install_python_module == 0:
            _install_script.write('BUILD_PYTHON_INTERFACE="n"; export BUILD_PYTHON_INTERFACE\n')                        
        if _prefix:
            _install_script.write('configure --prefix='+_prefix
                                  +' --exec-prefix='+_prefix+'\n')
        else:
            _install_script.write('configure\n')
        _install_script.write('make\nmake install\ncd ..\n\n\n')



def check_python():
    py_update_reason = ''
    v = sys.version_info
    vsrv =_ctinfo_srv.PythonVersion
    _srv['python_version'] = vsrv
    
    need_upgrade = 0
    if _bininstall and (v[0] <> vsrv[0] or v[1] <> vsrv[1]):
        py_update_reason = """
        Binary installs of Cantera requires Python version
        """+`vsrv[0]`+'.'+`vsrv[1]`+""", but this is version """+`v[0]`+'.'+`v[1]`+""".
        While any Python version 2.x can be used if you build Cantera
        from the source, binary installs require that the Python
        version matches that used when the binary distribution was
        created.  """
        need_upgrade = 1
        
    elif not check_write():
        
        py_update_reason = """
        You do not have write access to the library directory
        associated with this Python executable
        ("""+sys.prefix+"""/lib), and so the Python packages Cantera,
        Numeric, and MixMaster can't be installed.

        Options:
        
        a) quit and re-run this script as super-user (unix)
           or Administrator (Windows); or
           
        b) use this script to install a local version of Python
           in a directory where you have write access
           (won't work on Windows); or

        c) install only the Matlab and/or C++ components.
        
        """
        need_upgrade = 1


    elif force_update:
        need_upgrade = 1
        py_update_reason = '\nforced update\n'
        
    if need_upgrade:
        _updates.append(('Python '+`vsrv[0]`+'.'+`vsrv[1]`,get_python))
        _reason.append(py_update_reason)
    return need_upgrade


def check_numeric():
    if platform == 'win32' and _bininstall == 1:
        nv = _ctinfo_srv.WinNumericVersion
    else:
        nv = _ctinfo_srv.SrcNumericVersion
        
    num_update_reason = '--'
    
    try:
        import Numeric
        v = Numeric.__version__
        if _bininstall and v <> nv:
            num_update_reason = """
            
        Binary installs of Cantera requires Numeric (NumPy) version
        """+nv+""", but this system has version """+v+""".  While any
        version can be used if you build Cantera from the source,
        binary installs require that the version matches that used
        when the binary distribution was created.
        """

    except:
        num_update_reason = """
        The Cantera Python package requires the 'Numeric Extensions for
        Python' package ('Numeric'), but it is not installed on this
        system.
        """ 

    if force_update and num_update_reason == '--':
        num_update_reason = '\nForced update.\n'

    if num_update_reason <> '--':
        _updates.append(('Numeric Extensions for Python '+nv,get_numeric))
        _reason.append(num_update_reason)
        return 1
    else:
        return 0






print """
           Cantera Installation / Upgrade Utility
                     version 1.0

               type 'h' for help, or 'q' to quit
"""


platform = sys.platform
read_options()

#
# process command-line arguments
# 
args = sys.argv
nargs = len(args)

if nargs >= 2:
    for n in range(nargs):
        if  args[n] == '--force':
            force_update = 1
        elif args[n] == '--help':
            print __doc__
            sys.exit(0)
        elif args[n] == '--src':
            _bininstall = 0
        elif args[n] == '--python':
            install_python_module = 1
        elif args[n] == '--nopython':
            install_python_module = 0
        elif args[n] == '--matlab':
            install_matlab_toolbox = 1            
        elif args[n] == '--nomatlab':
            install_matlab_toolbox = 0            

_options['build_python'] = install_python_module
_options['build_matlab'] = install_matlab_toolbox


# open the install script files.
if sys.platform == 'win32' and _bininstall < 0:
    _bininstall = 1
    _install_script = open('install.bat','w')
    _install_script.write("""@echo off
REM run this script to install or update Cantera.
""")
else:
    _bininstall = 0
    _install_script = open('install.sh','w')
    _install_script.write("""#!/bin/sh
# run this script to install or update Cantera.
""")
    
print 'checking for updates...'
    
# get the information file from the server. This contains information about
# the versions on the server.
urllib.urlretrieve(_CTINFO_FILE, '_ctinfo_srv.py')
import _ctinfo_srv


#########################################################################


_thisdir = os.getcwd()


_files = []


##################################################################
#
#   check the Python version
#
##################################################################

if install_python_module:
    need_update = check_python()
    if need_update:
        install_python_module = 0

    
##################################################################
#
#   install/update Numeric
#
##################################################################

if install_python_module:
    check_numeric()


##################################################################
#
#   install/update Graphviz
#
##################################################################


if install_python_module:
    if platform == 'win32' and _bininstall == 1:
        gv = _ctinfo_srv.WinGraphvizVersion
    else:
        gv = _ctinfo_srv.SrcGraphvizVersion
    _srv['dot_version'] = gv
    graphviz_update_reason = '--'
    need_dot = 0
    dot_update_reason = '--'
    v = ''
    if _local.has_key('dot_version'):
        v = _local['dot_version']
    else:
        err = os.system('dot -V')
        if err == 0:
            print """
        Dot is installed, but the local configuration file
        does not contain version information.
        """
            v = ask('Enter dot version:','','','')
            _options['dot_version'] = v
            
    if v <> '':
        if v <> gv:
            need_dot = 1
            dot_update_reason = """
        A newer version of 'dot' is available.
        Installed dot version: """+v+""",
        Latest dot version:    """+gv+"""
            """
    else:
        need_dot = 1
        dot_update_reason = """
        Dot is not installed.
        """

    if force_update and dot_update_reason == '--':
        need_dot = 1        
        dot_update_reason = '\nforced update\n'
        
    if need_dot == 1:
        _updates.append(('Graphviz version '+gv, get_graphviz))
        _reason.append(dot_update_reason)
        


##################################################################
#
#   install/update Cantera
#
##################################################################
import types
ct_update_reason = 'none'

if _local.has_key('cantera'):
    v = float(_local['cantera'])
else:
    v = 0.0
srvv = float(_ctinfo_srv.CanteraDate)
_srv['cantera'] = `srvv`
if v <> 0.0 and (abs(v - srvv) > 1.0):
    ct_update_reason = """
    A newer version of Cantera is available.
    Installed version date: """+time.ctime(v)+"""
    Latest version date:    """+time.ctime(srvv)
    _updates.append(('Cantera',get_ct))
    _reason.append(ct_update_reason)
elif v == 0:
    ct_update_reason = """
    Cantera is not installed, or version information can't be found.
    """
    _updates.append(('Cantera',get_ct))
    _reason.append(ct_update_reason)        
elif force_update:
    ct_update_reason = '\nforced update\n'
    _updates.append(('Cantera',get_ct))
    _reason.append(ct_update_reason)            
#except:
#    _updates.append(('Cantera',get_ct)) 


        
urllib.urlcleanup()
os.remove('_ctinfo_srv.pyc')

ninst = len(_updates)

if ninst > 0:
    print '\n'+`ninst`+'  update(s) found.\n'

    i = 0
    for u in _updates:
        i = i + 1
        print '--------------------------------------------------------\n'
        print ' ['+`i`+'] '+u[0] + '\nReason for update: '+_reason[i-1]


    print '----------------------------------------------------------\n'
    print ' The following packages will be downloaded and installed:'
     
    a = ''
    print
    _selected = [' ']*ninst
    while a <> 'q':
        print
        i = 0
        for u in _updates:
            i = i + 1
            print _selected[i-1]+' ['+`i`+'] '+u[0]

        a = ''
        print
        print """
        Enter 'a' to select all packages, or enter a package number.
        Press <Enter> to begin download.
        """
        a = ask('Package:','1-'+`i`+'/a/h/q','')
        if a == '':
            break
        if a == 'a':
            for n in range(ninst):
                _selected[n] = '*'
        else:
            #try:
                nt = int(a)
                if nt > 0 and nt <= ninst:
                    if _selected[nt-1] == '*':
                        _selected[nt-1] = ' '
                    else:
                        _selected[nt-1] = '*'
            #except:
            #    pass

    if _srcdir == '--':
        print '\n\n'
        srcdir = get_srcdir()
        if srcdir: _install_script.write('cd '+srcdir+'\n')
        if _bininstall == 0:
            instdir = get_instdir()

        
    n = 0
    for u in _updates:
        n = n + 1            
        if _selected[n-1] == '*':
            print 'downloading '+u[0]
            _num = 0
            u[1]()
            print 'done.\n\n'


    _install_script.close()
    write_options()

    hlp = """
    The shell script 'install.sh' must be run to install the packages
    downloaded. If you answer 'y' or hit <Enter>, this script will
    be run now. If you enter 'n' or 'q', nothing will be installed, but you
    can later run 'install.sh' yourself to install the packages.
    """
    a = ask('Install the downloaded packages?','y/n/h/q','y',hlp)
    if a == 'y':
        if _new_python > 0:
            print """

************************ NOTE ***************************************

You are installing a new Python interpreter. After the installation
finishes, you need to run this script again with the new interpreter
in order to set up the Cantera Python interface.

*********************************************************************

"""
        os.chdir(_thisdir)
        if sys.platform == 'win32':
            os.execl('install.bat','install.bat')
        else:
            os.execl('/bin/sh','/bin/sh','install.sh')
else:
    write_options()    
    print 'Cantera is up to date.'




