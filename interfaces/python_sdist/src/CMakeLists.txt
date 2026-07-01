include(FetchContent)

if(CANTERA_PYODIDE)
  set(CANTERA_ENABLE_LAPACK OFF)
  set(CT_USE_HDF5 0)
  set(CT_USE_SYSTEM_HIGHFIVE 0)
else()
  set(CANTERA_ENABLE_LAPACK ON)
  set(CT_USE_HDF5 1)
  set(CT_USE_SYSTEM_HIGHFIVE 1)
endif()

set(EIGEN_VERSION 3.4.0)
FetchContent_Declare(
  eigen
  URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.bz2
  URL_HASH SHA256=b4c198460eba6f28d34894e3a5710998818515104d6e74e5cc331ce31e46e626
  # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
  # everything after it to pass to `find_package()`
  FIND_PACKAGE_ARGS NAMES Eigen3
)

set(YAML_CPP_VERSION 0.8.0)
FetchContent_Declare(
  yaml-cpp
  URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/${YAML_CPP_VERSION}.tar.gz
  URL_HASH SHA256=fbe74bbdcee21d656715688706da3c8becfd946d92cd44705cc6098bb23b3a16
  # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
  # everything after it to pass to `find_package()`
  FIND_PACKAGE_ARGS NAMES yaml-cpp
)

set(SUNDIALS_VERSION 7.4.0)
FetchContent_Declare(
  sundials
  URL https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz
  URL_HASH SHA256=679ddacdd77610110e613164e8297d6d0cd35bae8e9c3afc8e8ff6f99a1c2a7b
  # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
  # everything after it to pass to `find_package()`
  FIND_PACKAGE_ARGS NAMES SUNDIALS
)

set(FMT_VERSION 11.2.0)
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/releases/download/${FMT_VERSION}/fmt-${FMT_VERSION}.zip
  URL_HASH SHA256=203eb4e8aa0d746c62d8f903df58e0419e3751591bb53ff971096eaa0ebd4ec3
  # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
  # everything after it to pass to `find_package()`
  FIND_PACKAGE_ARGS NAMES fmt GLOBAL
)

if(CT_USE_HDF5)
  set(HIGHFIVE_VERSION .3.1.1)
  FetchContent_Declare(
    HighFive
    URL https://github.com/highfive-devs/HighFive/archive/refs/tags/v${HIGHFIVE_VERSION}.tar.gz
    URL_HASH SHA256=c0bc76b01868a133bf4550a07321923c5582fe3967edc102864b082e40799914
    # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
    # everything after it to pass to `find_package()`
    FIND_PACKAGE_ARGS NAMES HighFive
  )
endif()

if(NOT DEFINED Boost_INCLUDE_DIRS)
  find_package(Boost 1.70 CONFIG REQUIRED)
endif()
message(STATUS "Using Boost at ${Boost_INCLUDE_DIRS}")

# Common options
set(BUILD_TESTING OFF)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_BUILD_TYPE "Release")
# We don't care about Fortran support here so just disable this check
set(CMAKE_Fortran_COMPILER NOTFOUND)
# This is set to OFF to resolve linker errors on Windows and macOS. It doesn't
# seem to affect Linux builds.
set(BUILD_SHARED_LIBS OFF)

# HighFive build options
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)

# SUNDIALS build options
set(EXAMPLES_INSTALL OFF)
set(EXAMPLES_ENABLE_C OFF)
set(BUILD_ARKODE OFF)
set(BUILD_CVODE OFF)
set(BUILD_IDA OFF)
set(BUILD_KINSOL OFF)
set(BUILD_CPODES OFF)
set(BUILD_FORTRAN_MODULE_INTERFACE OFF)
if(CANTERA_ENABLE_LAPACK AND "${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    set(ENABLE_LAPACK ON)
    set(BLA_VENDOR Apple)
    set(SUNDIALS_LAPACK_CASE "LOWER")
    set(SUNDIALS_LAPACK_UNDERSCORES "NONE")
elseif(CANTERA_ENABLE_LAPACK AND "${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    set(ENABLE_LAPACK ON)
    set(BLA_VENDOR OpenBLAS)
else()
    set(ENABLE_LAPACK OFF)
endif()

# fmt build options
set(FMT_INSTALL OFF)
set(FMT_DOC OFF)
set(FMT_TEST OFF)
set(FMT_MASTER_PROJECT OFF)

# Eigen options
set(EIGEN_BUILD_DOC OFF)

# yaml-cpp options
set(YAML_CPP_FORMAT_SOURCE OFF)
set(YAML_CPP_INSTALL OFF)

set(CT_FETCHCONTENT_TARGETS eigen yaml-cpp sundials fmt)
if(CT_USE_HDF5)
  list(APPEND CT_FETCHCONTENT_TARGETS HighFive)
endif()
FetchContent_MakeAvailable(${CT_FETCHCONTENT_TARGETS})

file(GLOB_RECURSE CT_LIB_SOURCES "*.cpp")
# The Python extension handler is built with the rest of the Python library because it
# needs to be compiled after the Cython step runs.
list(FILTER CT_LIB_SOURCES EXCLUDE REGEX "extensions/.*")
list(APPEND CT_LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/extensions/canteraShared.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/extensions/pythonShim.cpp")

if(CANTERA_ENABLE_LAPACK AND NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)
    set(CT_USE_LAPACK 1)
    set(CT_SUNDIALS_USE_LAPACK 1)
    set (SUNDIALS_LIBRARIES BLAS::BLAS LAPACK::LAPACK SUNDIALS::sunlinsollapackdense SUNDIALS::sunlinsollapackband)
else()
    set (SUNDIALS_LIBRARIES SUNDIALS::sunlinsoldense SUNDIALS::sunlinsolband)
    set(CT_USE_LAPACK 0)
    set(CT_SUNDIALS_USE_LAPACK 0)
endif()

configure_file("../include/cantera/base/config.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/../include/cantera/base/config.h")

# On native platforms (Linux/macOS/Windows) libcantera is a SHARED library, so the
# Cython pure-Python-syntax modules can be standalone extensions that link it and ship
# their .py source beside the compiled .so. On Pyodide this split layout is not viable:
# each extension and libcantera would be an emscripten "side module", and emscripten
# cannot resolve C++ vague-linkage *data* symbols (RTTI typeinfo / vtables, which
# exception handling depends on) across the side-module boundary. So on Pyodide we fall
# back to a STATIC libcantera merged into a single _cantera extension (see
# cantera/CMakeLists.txt). Revisit if emscripten's dynamic linking of data symbols
# improves.
if(CANTERA_PYODIDE)
  add_library(cantera_lib STATIC ${CT_LIB_SOURCES})
else()
  add_library(cantera_lib SHARED ${CT_LIB_SOURCES})
endif()
target_include_directories(cantera_lib PRIVATE "../include")
target_include_directories(cantera_lib SYSTEM PRIVATE "${Python_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}")
set_target_properties(cantera_lib PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    # Cantera's headers don't annotate the public API with dllexport/dllimport, so
    # have CMake generate a .def exporting all symbols from the shared library
    WINDOWS_EXPORT_ALL_SYMBOLS ON)
# pythonShim.cpp (embedded interpreter support) calls into the Python C API, so the
# shared library must resolve those symbols at link time on Windows. Python::Module
# links the Python import library on Windows and adds nothing on POSIX, where the
# symbols are provided by the interpreter hosting the extension.
target_link_libraries(cantera_lib PRIVATE Python::Module)

set(CT_DEPS
    yaml-cpp::yaml-cpp
    Eigen3::Eigen
    SUNDIALS::core
    SUNDIALS::nvecserial
    SUNDIALS::cvodes
    SUNDIALS::idas
    SUNDIALS::sunmatrixband
    SUNDIALS::sunmatrixdense
    SUNDIALS::sunmatrixsparse
    SUNDIALS::sunlinsolspgmr
    SUNDIALS::sunnonlinsolnewton
    ${SUNDIALS_LIBRARIES})
if(CT_USE_HDF5)
  list(APPEND CT_DEPS HighFive)
endif()

if(CANTERA_PYODIDE)
  # Static libcantera is merged into the single _cantera extension, which therefore links
  # the third-party deps transitively. Link them (and fmt) PUBLIC so the merged extension
  # picks them up.
  target_link_libraries(cantera_lib PUBLIC ${CT_DEPS} fmt::fmt)
else()
  # Shared libcantera: bake the deps in (PRIVATE) so the consumer extensions link only
  # cantera_lib rather than carrying a duplicate copy of each static archive (e.g. the
  # whole SUNDIALS stack) -- this mirrors the SCons build and matters as more modules are
  # split out (thermo, kinetics, oneD, ...). The deps' headers are still part of
  # Cantera's public API, so propagate their *compile* usage requirements (include
  # directories and compile definitions such as SUNDIALS_STATIC_DEFINE) to consumers via
  # $<COMPILE_ONLY>, which evaluates to nothing on any link line (CMake >= 3.27).
  target_link_libraries(cantera_lib PRIVATE ${CT_DEPS})
  target_link_libraries(cantera_lib PUBLIC "$<COMPILE_ONLY:${CT_DEPS}>")
  # fmt is the exception: Cantera's public headers contain templated/inlined fmt calls, so
  # consumers reference fmt's compiled symbols and must link their own copy, as the SCons
  # build does. fmt is small and self-contained, so the duplicate is cheap.
  target_link_libraries(cantera_lib PUBLIC fmt::fmt)

  # Ship the shared libcantera inside the package, beside the extension modules that link
  # against it. On Windows the loader searches the directory containing the .pyd for
  # dependent DLLs; on Linux/macOS the extensions carry an RPATH of $ORIGIN/@loader_path
  # (set where the extensions are defined).
  install(TARGETS cantera_lib
      LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/cantera
      RUNTIME DESTINATION ${SKBUILD_PLATLIB_DIR}/cantera)
endif()
