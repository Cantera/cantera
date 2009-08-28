# Microsoft Developer Studio Project File - Name="cantera" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=cantera - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "cantera.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "cantera.mak" CFG="cantera - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "cantera - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "cantera - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "cantera - Win32 Release"

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
# ADD CPP /nologo /MD /W3 /GX /O2 /Ob2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /D "USE_MKL" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\cantera.lib"

!ELSEIF  "$(CFG)" == "cantera - Win32 Debug"

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
# ADD LIB32 /nologo /out:"..\..\build\lib\i686-pc-win32\cantera_d.lib"

!ENDIF 

# Begin Target

# Name "cantera - Win32 Release"
# Name "cantera - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\..\Cantera\src\BandMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ChemEquil.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ConstDensityThermo.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Constituents.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ct2ctml.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ctml.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ctvector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\CVode.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\DenseMatrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\EdgeKinetics.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Elements.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\FalloffFactory.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\funcs.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\GasKinetics.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\GasKineticsWriter.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\GRI_30_Kinetics.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Group.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\IdealGasPhase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ImplicitSurfChem.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\importCTML.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\InterfaceKinetics.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\KineticsFactory.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\misc.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Phase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\phasereport.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\plots.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\PureFluidPhase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ReactionPath.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ReactionStoichMgr.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\sort.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\SpeciesThermoFactory.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\State.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\StoichSubstance.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\stringUtils.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\SurfPhase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ThermoFactory.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ThermoPhase.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\xml.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=..\..\Cantera\src\Array.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\BandMatrix.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ChemEquil.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ck2ctml.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\config.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ConstDensityThermo.h
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\Constituent.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Constituents.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ct_defs.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ctDebug.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ctexceptions.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ctlapack.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ctml.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ctvector.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvdense.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvdiag.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\CVode.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvode.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\cvspgmr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\DASPK.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\dense.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\DenseMatrix.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\EdgeKinetics.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\EdgePhase.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Elements.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Enhanced3BConc.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\exceptions.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Falloff.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\FalloffFactory.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\FalloffMgr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\FuncEval.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\GasKinetics.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\GasKineticsWriter.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\global.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\GRI_30_Kinetics.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Group.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\IdealGasPhase.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ImplicitSurfChem.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\importCK.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\importCTML.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\importXML.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Integrator.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\InterfaceKinetics.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\iterativ.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Jac.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Jac1D.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Jac2.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Kinetics.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\KineticsFactory.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\lapack.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\llnlmath.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\llnltyps.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\mix_defs.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\MixFactory.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\MixTransport.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\MMCollisionInt.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\MultiJac.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\MultiResid.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\MultiTransport.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\NasaPoly1.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\NasaThermo.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Newton.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Newton1D.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\nvector.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Phase.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\plots.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\polyfit.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\PropertyCalculator.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\PropertyUpdater.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\pureSubstances.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\RateCoeffMgr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\reaction_defs.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ReactionData.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ReactionMechanism.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ReactionPath.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ReactionStoichMgr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Reactor.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ReactorBase.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\recipes.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Resid1D.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\RxnRates.h
# End Source File
# Begin Source File

SOURCE=..\..\CKReader\src\RxnSpecies.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\sort.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\SpeciesThermo.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\SpeciesThermoFactory.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\SpeciesThermoMgr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\speciesThermoTypes.h
# End Source File
# Begin Source File

SOURCE=..\..\ext\cvode\include\spgmr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\State.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\StFlow.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\StoichManager.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\stringUtils.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\surfKinetics.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\SurfPhase.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\TempCacher.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\TempFuncMgr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ThermoFactory.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ThermoPhase.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\ThirdBodyMgr.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\TigerPolynomial.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\TransportFactory.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\units.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\utilities.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\vec_functions.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\Wall.h
# End Source File
# Begin Source File

SOURCE=..\..\Cantera\src\xml.h
# End Source File
# End Group
# End Target
# End Project
