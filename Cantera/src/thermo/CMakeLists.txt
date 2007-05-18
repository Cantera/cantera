SET (CTTHERMO_SRCS State.cpp Elements.cpp Constituents.cpp Phase.cpp 
             ThermoPhase.cpp IdealGasPhase.cpp ConstDensityThermo.cpp 
             SpeciesThermoFactory.cpp ConstCpPoly.cpp 
             Mu0Poly.cpp GeneralSpeciesThermo.cpp SurfPhase.cpp 
             ThermoFactory.cpp phasereport.cpp  StoichSubstance.cpp 
             PureFluidPhase.cpp LatticeSolidPhase.cpp LatticePhase.cpp)

INCLUDE_DIRECTORIES (../base)

ADD_LIBRARY(ctthermo ${CTTHERMO_SRCS})
