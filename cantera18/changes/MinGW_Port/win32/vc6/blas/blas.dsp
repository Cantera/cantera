# Microsoft Developer Studio Project File - Name="blas" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=blas - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "blas.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "blas.mak" CFG="blas - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "blas - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "blas - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "blas - Win32 Release"

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
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\ctblas.lib"

!ELSEIF  "$(CFG)" == "blas - Win32 Debug"

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
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\ctblas_d.lib"

!ENDIF 

# Begin Target

# Name "blas - Win32 Release"
# Name "blas - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\..\ext\blas\dasum.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\daxpy.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dcabs1.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dcopy.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\ddot.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dgbmv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dgemm.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dgemv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dger.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dnrm2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\drot.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\drotg.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\drotm.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\drotmg.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsbmv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dscal.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsdot.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dspmv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dspr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dspr2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dswap.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsymm.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsymv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsyr.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsyr2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsyr2k.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dsyrk.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtbmv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtbsv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtpmv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtpsv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtrmm.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtrmv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtrsm.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dtrsv.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dzasum.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\dznrm2.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\icamax.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\idamax.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\isamax.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\izamax.f
# End Source File
# Begin Source File

SOURCE=..\..\ext\blas\xerbla.f
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# End Target
# End Project
