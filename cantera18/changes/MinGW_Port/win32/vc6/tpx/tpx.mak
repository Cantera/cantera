# Microsoft Developer Studio Generated NMAKE File, Based on tpx.dsp
!IF "$(CFG)" == ""
CFG=tpx - Win32 Debug
!MESSAGE No configuration specified. Defaulting to tpx - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "tpx - Win32 Release" && "$(CFG)" != "tpx - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "tpx.mak" CFG="tpx - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "tpx - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "tpx - Win32 Debug" (based on "Win32 (x86) Static Library")
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

!IF  "$(CFG)" == "tpx - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release

ALL : "..\..\lib\tpx.lib"


CLEAN :
	-@erase "$(INTDIR)\Hydrogen.obj"
	-@erase "$(INTDIR)\Methane.obj"
	-@erase "$(INTDIR)\Nitrogen.obj"
	-@erase "$(INTDIR)\Oxygen.obj"
	-@erase "$(INTDIR)\Sub.obj"
	-@erase "$(INTDIR)\utils.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\Water.obj"
	-@erase "..\..\lib\tpx.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/assume:underscore /compile_only /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /nologo /reentrancy:threaded /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/
CPP_PROJ=/nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\tpx.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\tpx.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"..\..\lib\tpx.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Methane.obj" \
	"$(INTDIR)\Nitrogen.obj" \
	"$(INTDIR)\Oxygen.obj" \
	"$(INTDIR)\Sub.obj" \
	"$(INTDIR)\utils.obj" \
	"$(INTDIR)\Water.obj" \
	"$(INTDIR)\Hydrogen.obj"

"..\..\lib\tpx.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "tpx - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "..\..\lib\tpx_d.lib"


CLEAN :
	-@erase "$(INTDIR)\Hydrogen.obj"
	-@erase "$(INTDIR)\Methane.obj"
	-@erase "$(INTDIR)\Nitrogen.obj"
	-@erase "$(INTDIR)\Oxygen.obj"
	-@erase "$(INTDIR)\Sub.obj"
	-@erase "$(INTDIR)\utils.obj"
	-@erase "$(INTDIR)\vc60.idb"
	-@erase "$(INTDIR)\vc60.pdb"
	-@erase "$(INTDIR)\Water.obj"
	-@erase "..\..\lib\tpx_d.lib"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/assume:underscore /check:bounds /compile_only /debug:full /iface:nomixed_str_len_arg /iface:cref /libs:dll /names:lowercase /nologo /reentrancy:threaded /traceback /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/
CPP_PROJ=/nologo /MD /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /Fp"$(INTDIR)\tpx.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ  /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\tpx.bsc" 
BSC32_SBRS= \
	
LIB32=link.exe -lib
LIB32_FLAGS=/nologo /out:"..\..\lib\tpx_d.lib" 
LIB32_OBJS= \
	"$(INTDIR)\Methane.obj" \
	"$(INTDIR)\Nitrogen.obj" \
	"$(INTDIR)\Oxygen.obj" \
	"$(INTDIR)\Sub.obj" \
	"$(INTDIR)\utils.obj" \
	"$(INTDIR)\Water.obj" \
	"$(INTDIR)\Hydrogen.obj"

"..\..\lib\tpx_d.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
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
!IF EXISTS("tpx.dep")
!INCLUDE "tpx.dep"
!ELSE 
!MESSAGE Warning: cannot find "tpx.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "tpx - Win32 Release" || "$(CFG)" == "tpx - Win32 Debug"
SOURCE=..\..\ext\tpx\Hydrogen.cpp

"$(INTDIR)\Hydrogen.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\ext\tpx\Methane.cpp

"$(INTDIR)\Methane.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\ext\tpx\Nitrogen.cpp

"$(INTDIR)\Nitrogen.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\ext\tpx\Oxygen.cpp

"$(INTDIR)\Oxygen.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\ext\tpx\Sub.cpp

"$(INTDIR)\Sub.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\ext\tpx\utils.cpp

"$(INTDIR)\utils.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)


SOURCE=..\..\ext\tpx\Water.cpp

"$(INTDIR)\Water.obj" : $(SOURCE) "$(INTDIR)"
	$(CPP) $(CPP_PROJ) $(SOURCE)



!ENDIF 

