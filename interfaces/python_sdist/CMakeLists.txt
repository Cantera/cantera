# Set minimum version 3.27 so that CMP0144, search for packages by upper-cased
# environment variables will use upper-cased variables. OLD behavior is to ignore
# upper-cased variables
cmake_minimum_required(VERSION 3.27...4.4)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
option(
  CANTERA_PYODIDE
  "Build Cantera Python module with a feature set suitable for Pyodide/emscripten."
  OFF)

# 1. Locate the target WebAssembly Python headers
find_package(Python COMPONENTS Development.Module REQUIRED)

# 2. Extract NumPy include path directly from the host isolated environment
execute_process(
    COMMAND python -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)

# 3. Create the standard FindPython target and variable manually
set(Python_NumPy_INCLUDE_DIRS "${NumPy_INCLUDE_DIRS}")
if(NOT TARGET Python::NumPy)
    add_library(Python::NumPy INTERFACE IMPORTED)
    set_target_properties(Python::NumPy PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${NumPy_INCLUDE_DIRS}"
    )
endif()

add_subdirectory("src")
add_subdirectory("cantera")

install(TARGETS _cantera DESTINATION ${SKBUILD_PLATLIB_DIR}/cantera)
