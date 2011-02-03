# Microsoft Developer Studio Project File - Name="oneD" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=oneD - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "oneD.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "oneD.mak" CFG="oneD - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "oneD - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "oneD - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "oneD - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /assume:underscore /compile_only /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /nologo /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "../../Cantera/src" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\oneD.lib"

!ELSEIF  "$(CFG)" == "oneD - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /debug:full /libs:dll /nologo /threads /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /Gm /GX /ZI /Od /I "../../Cantera/src" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\oneD_d.lib"

!ENDIF 

# Begin Target

# Name "oneD - Win32 Release"
# Name "oneD - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\boundaries1D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\MultiJac.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\MultiNewton.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\newton_utils.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\OneDim.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\refine.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\Sim1D.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\StFlow.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\Inlet1D.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\Jac1D.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\MultiJac.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\MultiNewton.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\Newton1D.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\OneDim.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\Resid1D.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\oneD\StFlow.h
# End Source File
# End Group
# End Target
# End Project
