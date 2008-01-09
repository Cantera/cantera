SET (THERMO_SRCS State.cpp Elements.cpp Constituents.cpp Phase.cpp 
             ThermoPhase.cpp IdealGasPhase.cpp ConstDensityThermo.cpp 
             SpeciesThermoFactory.cpp ConstCpPoly.cpp Nasa9Poly1.cpp 
             Nasa9PolyMultiTempRegion.cpp 
             Mu0Poly.cpp GeneralSpeciesThermo.cpp SurfPhase.cpp 
             ThermoFactory.cpp phasereport.cpp  SemiconductorPhase.cpp 
             StoichSubstance.cpp PureFluidPhase.cpp LatticeSolidPhase.cpp 
             LatticePhase.cpp)
ADD_LIBRARY(thermo ${THERMO_SRCS})
