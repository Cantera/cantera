
import sys

bindir = '/home/goodwin/ct154/bin'
libdir = '/home/goodwin/dv/sf/cantera/build/lib/i686-pc-linux-gnu'
incdir = '/home/goodwin/dv/sf/cantera/build/include'
dflibdir = ''

libs = ['clib', 'oneD', 'zeroD', 'transport', 'cantera', 'recipes',
        'cvode', 'ctlapack', 'ctmath', 'ctblas', 'tpx']

f = open('setup.m','w')
f.write('cd cantera\nbuild_cantera\nexit\n')
f.close()

fb = open('cantera/build_cantera.m','w')
fb.write("""
disp('building Cantera..');
mex -I"""+incdir+""" private/ctmethods.cpp private/ctfunctions.cpp ...
    private/xmlmethods.cpp private/phasemethods.cpp  ...
    private/thermomethods.cpp private/kineticsmethods.cpp ...
    private/transportmethods.cpp private/reactormethods.cpp ...
    private/wallmethods.cpp private/flowdevicemethods.cpp ...
    private/onedimmethods.cpp private/surfmethods.cpp private/write.cpp ...
""")
s = ''
for lib in libs:
    s += '    '+libdir+'/'+lib+'.lib ...\n'
fb.write(s)
fb.write('    "'+dflibdir+'/dformd.lib" ...\n')
fb.write('    "'+dflibdir+'/dfconsol.lib" ...\n')
fb.write('    "'+dflibdir+'/dfport.lib" \n')
fb.close()

fp = open('cantera/ctbin.m','w')
fp.write("""function path = ctbin
path = '"""+bindir+"""';
""")
fp.close()

