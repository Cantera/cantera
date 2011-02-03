SET (THERMO_SRCS State.cpp Elements.cpp Constituents.cpp Phase.cpp 
             ThermoPhase.cpp IdealGasPhase.cpp ConstDensityThermo.cpp 
             SpeciesThermoFactory.cpp ConstCpPoly.cpp Nasa9Poly1.cpp 
             Nasa9PolyMultiTempRegion.cpp 
             Mu0Poly.cpp GeneralSpeciesThermo.cpp SurfPhase.cpp 
             ThermoFactory.cpp phasereport.cpp  SemiconductorPhase.cpp 
             StoichSubstance.cpp PureFluidPhase.cpp LatticeSolidPhase.cpp 
             LatticePhase.cpp)

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)

ADD_LIBRARY(thermo ${THERMO_SRCS})

TARGET_LINK_LIBRARIES (thermo ctbase)

SET(THERMO_H   State.h Elements.h Constituents.h Phase.h mix_defs.h 
             ThermoPhase.h IdealGasPhase.h ConstDensityThermo.h 
             SpeciesThermoFactory.h ThermoFactory.h 
             NasaPoly1.h NasaPoly2.h NasaThermo.h Nasa9Poly1.h 
             Nasa9PolyMultiTempRegion.h 
             ShomateThermo.h ShomatePoly.h ConstCpPoly.h 
             SimpleThermo.h SpeciesThermoMgr.h 
             SpeciesThermoInterpType.h 
             GeneralSpeciesThermo.h Mu0Poly.h 
             speciesThermoTypes.h SpeciesThermo.h SurfPhase.h 
             EdgePhase.h )
#             @phase_header_files@)

INSTALL_FILES(/include/cantera/kernel FILES ${THERMO_H})
