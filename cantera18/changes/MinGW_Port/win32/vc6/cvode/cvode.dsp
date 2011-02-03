# Microsoft Developer Studio Project File - Name="cvode" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=cvode - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cvode.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cvode.mak" CFG="cvode - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cvode - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cvode - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "cvode - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD F90 /assume:underscore /iface:nomixed_str_len_arg /iface:cref /libs:dll /math_library:fast /names:lowercase /threads
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /Ob2 /I "..\..\ext\cvode\include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\cvode.lib"

!ELSEIF  "$(CFG)" == "cvode - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "../../lib"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD F90 /assume:underscore /dbglibs /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /threads
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /Gm /GX /ZI /Od /I "..\..\ext\cvode\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\cvode_d.lib"

!ENDIF 

# Begin Target

# Name "cvode - Win32 Release"
# Name "cvode - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=..\..\ext\cvode\source\band.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\cvband.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\cvbandpre.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\cvdense.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\cvdiag.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\cvode.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\cvspgmr.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\dense.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\iterativ.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\llnlmath.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\nvector.c
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\source\spgmr.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\ext\cvode\include\band.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvband.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvbandpre.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvdense.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvdiag.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvode.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvspgmr.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\dense.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\iterativ.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\llnlmath.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\llnltyps.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\nvector.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\spgmr.h
# End Source File
# End Group
# End Target
# End Project
