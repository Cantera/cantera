
import sys

bindir = 'c:/cantera/bin'
libdir = 'd:/dgg/dv/sf/cantera/build/lib/i686-pc-win32'
incdir = 'd:/dgg/dv/sf/cantera/build/include'
libs   = ' -loneD -lzeroD -ltransport -lcantera -lrecipes -lcvode -lctlapack -lctmath -lctblas -ltpx  '

f = open('setup.m','w')
f.write('cd cantera\nbuildux\nexit\n')
f.close()

fb = open('cantera/buildux.m','w')
fb.write("""
disp('building Cantera..');
mex private/ctmethods.cpp private/ctfunctions.cpp ...
    private/xmlmethods.cpp private/phasemethods.cpp  ...
    private/thermomethods.cpp private/kineticsmethods.cpp ...
    private/transportmethods.cpp private/reactormethods.cpp ...
    private/wallmethods.cpp private/flowdevicemethods.cpp ...
    private/onedimmethods.cpp private/surfmethods.cpp private/write.cpp ...
"""+'-I'+incdir+'     -L'+libdir+' '+libs+'\n'+"""disp('done.');
""")
fb.close()

fp = open('cantera/ctbin.m','w')
fp.write("""function path = ctbin
path = '"""+bindir+"""';
""")
fp.close()

