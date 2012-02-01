h_sources  = State.h Elements.h Constituents.h Phase.h mix_defs.h  \
             ThermoPhase.h IdealGasPhase.h ConstDensityThermo.h    \
             SpeciesThermoFactory.h ThermoFactory.h                \
             NasaPoly1.h NasaPoly2.h NasaThermo.h Nasa9Poly1.h     \
	     Nasa9PolyMultiTempRegion.h WaterPropsIAPWS.h          \
             ShomateThermo.h ShomatePoly.h ConstCpPoly.h           \
             SimpleThermo.h SpeciesThermoMgr.h                     \
             SpeciesThermoInterpType.h VPSSMgr_Water_ConstVol.h    \
             GeneralSpeciesThermo.h Mu0Poly.h                      \
             speciesThermoTypes.h SpeciesThermo.h SurfPhase.h      \
             EdgePhase.h IdealSolnGasVPSS.h MolalityVPSSTP.h       \
             VPSSMgr.h  VPSSMgrFactory.h  VPSSMgr_General.h        \
             VPSSMgr_IdealGas.h VPSSMgr_ConstVol.h PDSS_ConstVol.h \
	     PDSS_IdealGas.h PDSS_SSVol.h DebyeHuckel.h            \
	     VPStandardStateTP.h PDSS_Water.h PDSS.h WaterProps.h  \
	     WaterPropsIAPWSphi.h VPSSMgr_Water_HKFT.h PDSS_HKFT.h \
             PDSS_IonsFromNeutral.h IonsFromNeutralVPSSTP.h        \
             GibbsExcessVPSSTP.h

cc_sources = State.cpp Elements.cpp Constituents.cpp Phase.cpp             \
             ThermoPhase.cpp IdealGasPhase.cpp ConstDensityThermo.cpp      \
             SpeciesThermoFactory.cpp ConstCpPoly.cpp Nasa9Poly1.cpp       \
	     Nasa9PolyMultiTempRegion.cpp PDSS_Water.cpp PDSS_HKFT.cpp     \
             Mu0Poly.cpp GeneralSpeciesThermo.cpp SurfPhase.cpp            \
             ThermoFactory.cpp phasereport.cpp SpeciesThermoInterpType.cpp \
             VPSSMgr.cpp  VPSSMgrFactory.cpp  VPSSMgr_General.cpp          \
	     IdealSolnGasVPSS.cpp MolalityVPSSTP.cpp VPStandardStateTP.cpp \
             VPSSMgr_IdealGas.cpp VPSSMgr_ConstVol.cpp PDSS_ConstVol.cpp   \
	     PDSS_IdealGas.cpp PDSS_SSVol.cpp DebyeHuckel.cpp PDSS.cpp     \
	     WaterProps.cpp WaterPropsIAPWS.cpp	WaterPropsIAPWSphi.cpp     \
             VPSSMgr_Water_HKFT.cpp VPSSMgr_Water_ConstVol.cpp             \
             PDSS_IonsFromNeutral.cpp IonsFromNeutralVPSSTP.cpp            \
             GibbsExcessVPSSTP.cpp  

# semiconductorphase
# stoich substance
# pure fluid phase
# lattice solide phase
# lattice phase
#
# g++ -c VPStandardStateTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c IdealMolalSoln.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c WaterPropsIAPWSphi.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c WaterPropsIAPWS.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c WaterProps.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c PDSS.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c PDSS_Water.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c PDSS_HKFT.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c HMWSoln.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c HMWSoln_input.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c DebyeHuckel.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c WaterSSTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c MetalSHEelectrons.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c VPSSMgr_Water_ConstVol.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c VPSSMgr_Water_HKFT.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c IdealSolidSolnPhase.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c StoichSubstanceSSTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c SingleSpeciesTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c MineralEQ3.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c GibbsExcessVPSSTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c PseudoBinaryVPSSTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c MargulesVPSSTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c IonsFromNeutralVPSSTP.cpp -O3 -Wall     -fPIC  -I../base 
# g++ -c PDSS_IonsFromNeutral.cpp -O3 -Wall     -fPIC  -I../base 



	     PDSS_Water.cpp PDSS_HKFT.cpp     \
	     DebyeHuckel.cpp PDSS.cpp     \
	     WaterProps.cpp WaterPropsIAPWS.cpp	WaterPropsIAPWSphi.cpp     \
             VPSSMgr_Water_HKFT.cpp VPSSMgr_Water_ConstVol.cpp             \
             PDSS_IonsFromNeutral.cpp

AM_CPPFLAGS             = -I../base -I../../../
AM_CXXFLAGS             = $(AM_CPPFLAGS)

lib_LTLIBRARIES         = $(top_builddir)/build/lib/libthermo.la
library_includedir      = $(includedir)

#-----------------------
# Cantera thermo C/C++ library
#-----------------------

__top_builddir__build_lib_libthermo_la_LDFLAGS      = $(all_libraries) -release $(GENERIC_RELEASE)
__top_builddir__build_lib_libthermo_la_SOURCES      = $(cc_sources) $(h_sources)
