# Microsoft Developer Studio Project File - Name="lapack" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=lapack - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "lapack.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "lapack.mak" CFG="lapack - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "lapack - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "lapack - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "lapack - Win32 Release"

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
# ADD F90 /assume:underscore /compile_only /iface:nomixed_str_len_arg /iface:cref /libs:dll /math_library:fast /names:lowercase /nologo /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /Ob2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\ctlapack.lib"

!ELSEIF  "$(CFG)" == "lapack - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /assume:underscore /check:bounds /compile_only /dbglibs /debug:full /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /nologo /threads /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\ctlapack_d.lib"

!ENDIF 

# Begin Target

# Name "lapack - Win32 Release"
# Name "lapack - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\..\ext\lapack\dbdsqr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgbsv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgbtf2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgbtrf.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgbtrs.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgebd2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgebrd.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgelq2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgelqf.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgelss.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgeqr2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgeqrf.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgetf2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgetrf.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgetri.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dgetrs.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlabad.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlabrd.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlacpy.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlamch.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlange.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlapy2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlarf.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlarfb.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlarfg.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlarft.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlartg.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlas2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlascl.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlaset.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlasq1.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlasq2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlasq3.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlasq4.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlasr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlasrt.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlassq.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlasv2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dlaswp.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dorg2r.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dorgbr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dorgl2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dorglq.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dorgqr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dorm2r.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dormbr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dorml2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dormlq.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\dormqr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\drscl.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\ilaenv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\lapack\lsame.f
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# End Target
# End Project
