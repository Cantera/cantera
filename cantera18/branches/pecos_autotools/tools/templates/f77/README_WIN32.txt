CVF/VC++ 6.0 Build Procedure
============================

This demo program can be built on a PC using Compaq Visual Fortran 6.0
and Microsoft Visual C++ 6.0, which share a common development
environment (Visual Studio). To build it, open project file
f77demo.dsp in Visual Studio by double-clicking its icon. Build the
project, and then run it.

To modify this project file to build your own Fortran 77 application,
make a copy with a different name, open the copy in Visual Studio,
delete file demo.f and replace it with your Fortran program. If you
need to change demo_ftnlib.cpp, then also replace it with your
modified version.

Note that the include and library directories are specified relative
to the location of the project file. If you move it to another
directory, you will need to update these locations under 'Settings'.
