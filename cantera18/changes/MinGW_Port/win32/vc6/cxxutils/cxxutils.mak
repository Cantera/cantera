# Microsoft Developer Studio Generated NMAKE File, Based on cxxutils.dsp
!IF "$(CFG)" == ""
CFG=cxxutils - Win32 Debug
!MESSAGE No configuration specified. Defaulting to cxxutils - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "cxxutils - Win32 Release" && "$(CFG)" != "cxxutils - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cxxutils.mak" CFG="cxxutils - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cxxutils - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cxxutils - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "cxxutils - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : "..\..\lib\cxxutils.lib"


CLEAN :
	-@erase "$(INTDIR)\cxxutils.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\writelog.obj"
	-@erase "..\..\lib\cxxutils.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/compile_only /libs:dll /nologo /threads /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/
CPP_PROJ=/nologo /MD /W3 /GX /O2 /I "../../Cantera/src" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\cxxutils.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\cxxutils.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"..\..\lib\cxxutils.lib" 
LIB32_OBJS= \
	"$(INTDIR)\cxxutils.obj" \
	"$(INTDIR)\writelog.obj"

"..\..\lib\cxxutils.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "cxxutils - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "..\..\lib\cxxutils_d.lib"


CLEAN :
	-@erase "$(INTDIR)\cxxutils.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\writelog.obj"
	-@erase "..\..\lib\cxxutils_d.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/check:bounds /compile_only /debug:full /libs:dll /nologo /threads /traceback /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/
CPP_PROJ=/nologo /MD /W3 /Gm /GX /ZI /Od /I "../../Cantera/src" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\cxxutils.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\cxxutils.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"..\..\lib\cxxutils_d.lib" 
LIB32_OBJS= \
	"$(INTDIR)\cxxutils.obj" \
	"$(INTDIR)\writelog.obj"

"..\..\lib\cxxutils_d.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("cxxutils.dep")
!INCLUDE "cxxutils.dep"
!ELSE 
!MESSAGE Warning: cannot find "cxxutils.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "cxxutils - Win32 Release" || "$(CFG)" == "cxxutils - Win32 Debug"
SOURCE=..\..\Cantera\cxx\cxxutils.cpp

"$(INTDIR)\cxxutils.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\Cantera\cxx\writelog.cpp

"$(INTDIR)\writelog.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

