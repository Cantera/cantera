
import sys

bindir = 'C:/cygwin/usr/local/cantera/bin'
libdir = 'c:/vc_env/cant17/build/lib/i686-pc-win32'
incdir = 'c:/vc_env/cant17/build/include'
dflibdir = ''

bllibstr = "-lctlapack -lctblas"
bllibs = bllibstr.replace('-l',' ')
bllist = bllibs.split()

bldir = "c:/vc_env/cant17/build/lib/i686-pc-win32"

libs = ['clib', 'oneD', 'zeroD', 'transport', 'equil', 'kinetics', 'thermo', 'numerics',  'converters', 'base',
        'ctcxx', 'tpx', 'NVEC_SER', 'CVODES', 'SUNDIALS_SHARED',  'ctmath', 'ctlapack', 'ctblas', 'ctf2c']

f = open('setup.m','w')
f.write('cd cantera\nbuild_cantera\nexit\n')
f.close()

fb = open('cantera/build_cantera.m','w')
fb.write("""
disp('building Cantera..');
mex -v  -I"""+incdir+""" private/ctmethods.cpp private/ctfunctions.cpp ...
    private/xmlmethods.cpp private/phasemethods.cpp  ...
    private/thermomethods.cpp private/kineticsmethods.cpp ...
    private/transportmethods.cpp private/reactormethods.cpp ...
    private/reactornetmethods.cpp ...
    private/wallmethods.cpp private/flowdevicemethods.cpp ...
    private/funcmethods.cpp ...
    private/mixturemethods.cpp ...
    private/onedimmethods.cpp private/surfmethods.cpp  ...
""")
s = ''
for lib in libs:
    s += '    '+libdir+'/'+lib+'.lib ...\n'
fb.write(s)
if bllist:
    s = ''
    for lib in bllist:
        s += '    '+bldir+'/'+lib+'.lib ...\n'
    fb.write(s)
# fb.write('    "'+dflibdir+'/dformd.lib" ...\n')
# fb.write('    "'+dflibdir+'/dfconsol.lib" ...\n')
# fb.write('    "'+dflibdir+'/dfport.lib" \n')
fb.write('   \n')

fb.close()

fp = open('cantera/ctbin.m','w')
fp.write("""function path = ctbin
path = '"""+bindir+"""';
""")
fp.close()

