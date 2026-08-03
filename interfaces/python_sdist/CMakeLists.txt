cmake_minimum_required(VERSION 3.27...4.4)

project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
option(
  CANTERA_PYODIDE
  "Build Cantera Python module with a feature set suitable for Pyodide/emscripten."
  OFF
)

# Emscripten's dependency scanner cannot parse all compiler flags
# supplied by the Pyodide toolchain.
if(EMSCRIPTEN AND CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
      set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
endif()

# CMake 4.1 changes how cross compiling works and especially how host vs. target
# dependencies like Python are discovered. This means we need to manually disover
# NumPy. I'm not sure what needs to change in CMake, scikit-build-core, or pyodide
# to fix this behavior.
if (EMSCRIPTEN AND CMAKE_VERSION VERSION_GREATER_EQUAL 4.1)
  # Explicitly use modern cross-compilation behavior.
  cmake_policy(SET CMP0190 NEW)

  # scikit-build-core provides Python_EXECUTABLE pointing to the
  # isolated build environment's Python interpreter.
  set(_host_python "${Python_EXECUTABLE}")

  # Find target Python development headers from Pyodide.
  find_package(Python COMPONENTS Development.Module REQUIRED)

  # Discover NumPy using the build environment's Python.
  execute_process(
      COMMAND "${_host_python}" -c
          "import numpy; print(numpy.get_include())"
      RESULT_VARIABLE _numpy_result
      OUTPUT_VARIABLE Python_NumPy_INCLUDE_DIRS
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_VARIABLE _numpy_error
  )

  if(NOT _numpy_result EQUAL 0)
      message(FATAL_ERROR
          "Failed to find NumPy headers: ${_numpy_error}")
  endif()

  set(Python_NumPy_INCLUDE_DIR "${Python_NumPy_INCLUDE_DIRS}")
  set(Python_NumPy_FOUND TRUE)

  add_library(Python::NumPy INTERFACE IMPORTED)
  set_target_properties(Python::NumPy PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Python_NumPy_INCLUDE_DIRS}"
  )

  message(STATUS "Target Python: ${Python_INCLUDE_DIRS}")
  message(STATUS "Host Python: ${_host_python}")
  message(STATUS "NumPy includes: ${Python_NumPy_INCLUDE_DIRS}")
else()
    find_package(Python COMPONENTS Development.Module NumPy REQUIRED)
endif()

find_package(Cython MODULE REQUIRED VERSION 3)
include(UseCython)
add_subdirectory("src")
add_subdirectory("cantera")

install(TARGETS _cantera DESTINATION ${SKBUILD_PLATLIB_DIR}/cantera)
