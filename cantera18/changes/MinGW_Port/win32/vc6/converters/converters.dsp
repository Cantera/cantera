# Microsoft Developer Studio Project File - Name="converters" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=converters - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "converters.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "converters.mak" CFG="converters - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "converters - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "converters - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "converters - Win32 Release"

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
# ADD F90 /assume:underscore /compile_only /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /nologo /reentrancy:threaded /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "../../Cantera/src" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"../../build/lib/i686-pc-win32/converters.lib"

!ELSEIF  "$(CFG)" == "converters - Win32 Debug"

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
# ADD LIB32 /nologo /out:"../../build/lib/i686-pc-win32/converters_d.lib"

!ENDIF 

# Begin Target

# Name "converters - Win32 Release"
# Name "converters - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\..\Cantera\src\converters\atomicWeightDB.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\ck2ct.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\CKParser.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\ckr_utils.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\CKReader.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\Reaction.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\thermoFunctions.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\writelog.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=..\..\Cantera\src\converters\ck2ct.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\CKParser.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\ckr_defs.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\ckr_utils.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\CKReader.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\Reaction.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\thermoFunctions.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\converters\writelog.h
# End Source File
# End Group
# End Target
# End Project
