# TODO: Follow the progress of https://github.com/scikit-build/cython-cmake,
# when it seems stable switch to that rather than custom-compiling
file(GLOB_RECURSE CY_SOURCES "*.pyx")
set(CY_OUTPUTS)
foreach(CY_SOURCE IN LISTS CY_SOURCES)
  cmake_path(GET CY_SOURCE FILENAME CY_PATH)
  cmake_path(REPLACE_EXTENSION CY_PATH ".cpp" OUTPUT_VARIABLE CPP_OUTPUT)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CPP_OUTPUT} OUTPUT_VARIABLE CPP_OUTPUT)
  list(APPEND CY_OUTPUTS ${CPP_OUTPUT})
  if(${CY_SOURCE} MATCHES "delegator.pyx$")
    cmake_path(REPLACE_EXTENSION CPP_OUTPUT ".h" OUTPUT_VARIABLE DELEGATOR_H)
    list(APPEND CPP_OUTPUT ${DELEGATOR_H})
    add_custom_target(DELEGATOR_H_EXT ALL DEPENDS ${DELEGATOR_H})
  endif()

  add_custom_command(
    OUTPUT ${CPP_OUTPUT}
    COMMENT "Making ${CPP_OUTPUT} from ${CY_SOURCE}"
    COMMAND Python::Interpreter -m cython ${CY_SOURCE} --output ${CMAKE_CURRENT_BINARY_DIR} --cplus --directive binding=True
    DEPENDS ${CY_SOURCE}
    VERBATIM)
endforeach()

# The extension manager is compiled here rather than with the rest of the
# Cantera source code because it needs delegator.h, built by the Cython step
# above.
add_library(ext_manager ../src/extensions/PythonExtensionManager.cpp)
target_include_directories(ext_manager PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}"
  "../include"
  "${Python_INCLUDE_DIRS}")
target_link_libraries(ext_manager PRIVATE fmt::fmt)
add_dependencies(ext_manager DELEGATOR_H_EXT)

python_add_library(_cantera MODULE ${CY_OUTPUTS} WITH_SOABI)
target_include_directories(_cantera PRIVATE "${Python_NumPy_INCLUDE_DIRS}" "../include")
target_compile_definitions(_cantera PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
target_link_libraries(_cantera PRIVATE ext_manager cantera_lib)
