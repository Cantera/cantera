#
# python script to create a few key files for the windows build
#
# $Id: setup_winmatlab.py,v 1.29 2009/07/20 19:38:53 hkmoffa Exp $
#
import sys
import os
import os.path

#
# We are currently in the CanteraRoot/Cantera/matlab directory
CurrDir=os.getcwd()
#
# Get a list of path + 'matlab'
#
CantDir=os.path.split(CurrDir)
#
# Get a list of CanteraRoot + 'Cantera'
#
CantRootList=os.path.split(CantDir[0])
#
# CanteraRoot will contain the absolute path name
#             of the development root directory
#
CanteraRoot=CantRootList[0]

bindir = CanteraRoot + '/bin'
libdir = CanteraRoot + '/build/lib/i686-pc-win32'
incdir = CanteraRoot + '/build/include'
dflibdir = ''

bllibstr = "-lctlapack -lctblas"
bllibs = bllibstr.replace('-l',' ')
bllist = bllibs.split()

bldir = libdir

libs = ['clib']
#
#  setup.m  file
#     This is a utility matlab file that tells matlab to
#     jump down a directory and execute build_cantera.m
#
f = open('setup.m','w')
f.write('cd cantera\nbuild_cantera\nexit\n')
f.close()
#
#  cantera/build_Cantera.m
#
#    Here we create the file cantera/build_Cantera.m
#    This file contains the command which is executed from within matlab
#    to build the cantera extension
#
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
#
#  Here we create the file ctbin.m.
#  This is a matlab command which specified the ctbin directory
#  within matlab. However we have a problem. This refers to the
#  development computer, while we may be on the target computer.
#  Don't know how to resolve this atm.
#
fp = open('cantera/ctbin.m','w')
fp.write("""function path = ctbin
path = '"""+bindir+"""';
""")
fp.close()

