# Microsoft Developer Studio Project File - Name="ckreader" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ckreader - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ckreader.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ckreader.mak" CFG="ckreader - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ckreader - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "ckreader - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "ckreader - Win32 Release"

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
# ADD F90 /assume:underscore /iface:nomixed_str_len_arg /iface:cref /libs:static /math_library:fast /names:lowercase
# SUBTRACT F90 /threads
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /I "../../include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\lib\ckreader.lib"

!ELSEIF  "$(CFG)" == "ckreader - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD F90 /assume:underscore /dbglibs /iface:nomixed_str_len_arg /iface:cref /names:lowercase
# SUBTRACT F90 /threads
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "../../include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\lib\ckreader_d.lib"

!ENDIF 

# Begin Target

# Name "ckreader - Win32 Release"
# Name "ckreader - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\CKReader\src\atomicWeightDB.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\CKParser.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\ckr_utils.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\CKReader.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\filter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\Reaction.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\thermoFunctions.cpp
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\writelog.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\CKReader\src\CKParser.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\CKParser.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\ckr_defs.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\ckr_defs.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\ckr_utils.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\ckr_utils.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\CKReader.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\CKReader.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\config.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\config.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\Constituent.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\Constituent.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\Element.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\Element.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\Reaction.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\Reaction.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\Species.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\Species.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\thermoFunctions.h
# End Source File
# Begin Source File

SOURCE="..\..\dv\cantera-1.1b\CKReader\src\thermoFunctions.h"
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\writelog.h
# End Source File
# End Group
# End Target
# End Project
