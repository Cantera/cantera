
CONFIGURE_FILE (${PROJECT_SOURCE_DIR}/BackEnd/src/algorithms/DAEs/IDA_Solver.hpp.in  ${PROJECT_BINARY_DIR}/Headers/IDA_Solver.hpp)

INCLUDE_DIRECTORIES (${PROJECT_BINARY_DIR}/Headers)        
ADD_LIBRARY (cantera-dae IDA_Solver)
TARGET_LINK_LIBRARY (sundials)=

INSTALL_TARGETS(/lib)

