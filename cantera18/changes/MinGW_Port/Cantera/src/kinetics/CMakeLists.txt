SET (KINETICS_SRCS importKinetics.cpp GRI_30_Kinetics.cpp KineticsFactory.cpp
            GasKinetics.cpp FalloffFactory.cpp ReactionStoichMgr.cpp Kinetics.cpp 
            solveSP.cpp InterfaceKinetics.cpp ImplicitSurfChem.cpp Group.cpp 
            ReactionPath.cpp)

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/base)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/thermo)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/Cantera/src/numerics)

ADD_LIBRARY(kinetics ${KINETICS_SRCS})

TARGET_LINK_LIBRARIES (kinetics numerics thermo ctbase)

SET(KINETICS_H  importKinetics.h GRI_30_Kinetics.h KineticsFactory.h 
             Kinetics.h GasKinetics.h 
             FalloffFactory.h  ReactionStoichMgr.h reaction_defs.h 
             FalloffMgr.h ThirdBodyMgr.h RateCoeffMgr.h ReactionData.h 
             RxnRates.h Enhanced3BConc.h StoichManager.h solveSP.h InterfaceKinetics.h
             ImplicitSurfChem.h EdgeKinetics.h Group.h ReactionPath.h)

INSTALL_FILES(/cantera/kernel FILES ${KINETICS_H})
