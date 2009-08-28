# Microsoft Developer Studio Project File - Name="ctmath" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ctmath - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ctmath.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ctmath.mak" CFG="ctmath - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ctmath - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "ctmath - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "ctmath - Win32 Release"

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
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\ctmath.lib"

!ELSEIF  "$(CFG)" == "ctmath - Win32 Debug"

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
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /assume:underscore /check:bounds /compile_only /dbglibs /debug:none /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /nologo /optimize:0 /threads /traceback /warn:argument_checking /warn:nofileopt
# SUBTRACT F90 /automatic /fast
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\ctmath_d.lib"

!ENDIF 

# Begin Target

# Name "ctmath - Win32 Release"
# Name "ctmath - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\..\ext\math\daux.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\ddaspk.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\dgbefa.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\dgbsl.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\dgefa.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\dgesl.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\dp1vlu.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\dpcoef.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\dpolft.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\fdump.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\idamax.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\j4save.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\mach.cpp
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\pcoef.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\polfit.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\pvalue.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\xercnt.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\xerhlt.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\xermsg.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\xerprn.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\xersve.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\xgetua.f
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=..\..\ext\math\cblas.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\math\gmres.h
# End Source File
# End Group
# End Target
# End Project
