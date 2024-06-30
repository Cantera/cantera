
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/_cantera.cpp"
  COMMENT
    "Making ${CMAKE_CURRENT_BINARY_DIR}/_cantera.cpp from ${CMAKE_CURRENT_SOURCE_DIR}/_cantera.pyx"
  COMMAND Python::Interpreter -m cython
          "${CMAKE_CURRENT_SOURCE_DIR}/_cantera.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/_onedim.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/_utils.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/constants.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/delegator.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/func1.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/kinetics.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/mixture.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/preconditioners.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/reaction.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/reactionpath.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/reactor.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/solutionbase.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/speciesthermo.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/thermo.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/transport.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/units.pyx"
          "${CMAKE_CURRENT_SOURCE_DIR}/yamlwriter.pyx"
          --output "${CMAKE_CURRENT_BINARY_DIR}/_cantera.cpp" --cplus --directive binding=True
  DEPENDS _cantera.pyx _onedim.pyx _utils.pyx constants.pyx delegator.pyx func1.pyx kinetics.pyx mixture.pyx preconditioners.pyx reaction.pyx
    reactionpath.pyx reactor.pyx solutionbase.pyx speciesthermo.pyx thermo.pyx transport.pyx units.pyx yamlwriter.pyx
  VERBATIM)

python_add_library(cantera MODULE _cantera.cpp WITH_SOABI)
target_include_directories(cantera PRIVATE "${Python_NumPy_INCLUDE_DIRS}" "../include" "${CMAKE_CURRENT_BINARY_DIR}/../include" "${CMAKE_CURRENT_BINARY_DIR}/../include/cantera/base")
target_compile_definitions(cantera PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
target_link_libraries(cantera PRIVATE cantera_lib)
