
import sys

bindir = '/usr/local/bin'
libdir = '/Users/dgg/dv/sf/cantera/build/lib/powerpc-apple-darwin7.2.0'
incdir = '/Users/dgg/dv/sf/cantera/build/include'
libs   = '-lclib  -loneD -lzeroD -ltransport -lcantera -lrecipes -lcvode -lctlapack -lctmath -lctblas -ltpx  -lg2c -lgcc'

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
    private/funcmethods.cpp ... 
    private/onedimmethods.cpp private/surfmethods.cpp private/write.cpp ...
"""+'-I'+incdir+'     -L'+libdir+' '+libs+'\n'+"""disp('done.');
""")
fb.close()

fp = open('cantera/ctbin.m','w')
fp.write("""function path = ctbin
path = '"""+bindir+"""';
""")
fp.close()

