# Microsoft Developer Studio Project File - Name="zeroD" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=zeroD - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "zeroD.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "zeroD.mak" CFG="zeroD - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "zeroD - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "zeroD - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "zeroD - Win32 Release"

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
# ADD CPP /nologo /MD /W3 /GX /O2 /Ob2 /I "../../Cantera/src" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\zeroD.lib"

!ELSEIF  "$(CFG)" == "zeroD - Win32 Debug"

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
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\zeroD_d.lib"

!ENDIF 

# Begin Target

# Name "zeroD - Win32 Release"
# Name "zeroD - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\FlowDevice.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\Reactor.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\ReactorBase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\ReactorNet.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\Wall.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\flowControllers.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\FlowDevice.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\PID_Controller.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\Reactor.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\ReactorBase.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\ReactorNet.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\Reservoir.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\zeroD\Wall.h
# End Source File
# End Group
# End Target
# End Project
