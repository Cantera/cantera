
import sys

if len(sys.argv) >= 3:
    libdir = sys.argv[1]
    libs = '-l'+sys.argv[2]+' '+sys.argv[3]
else:
    print 'usage: python setup_matlab.py <libdir> <lib>'
    sys.exit(0)
    
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
    private/onedimmethods.cpp private/write.cpp ...
"""+'     -L'+libdir+' '+libs+'\n'+"""disp('done.');
""")
fb.close()
