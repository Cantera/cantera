include(FetchContent)

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

set(SUNDIALS_VERSION 7.0.0)
FetchContent_Declare(
  sundials
  URL https://github.com/LLNL/sundials/releases/download/v${SUNDIALS_VERSION}/sundials-${SUNDIALS_VERSION}.tar.gz
  URL_HASH SHA256=d762a7950ef4097fbe9d289f67a8fb717a0b9f90f87ed82170eb5c36c0a07989
  # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
  # everything after it to pass to `find_package()`
  FIND_PACKAGE_ARGS 7.0.0 EXACT NAMES SUNDIALS
)

set(FMT_VERSION 10.2.1)
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/releases/download/${FMT_VERSION}/fmt-${FMT_VERSION}.zip
  URL_HASH SHA256=312151a2d13c8327f5c9c586ac6cf7cddc1658e8f53edae0ec56509c8fa516c9
  # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
  # everything after it to pass to `find_package()`
  FIND_PACKAGE_ARGS NAMES fmt GLOBAL
)

set(HIGHFIVE_VERSION 2.9.0)
FetchContent_Declare(
  HighFive
  URL https://github.com/BlueBrain/HighFive/archive/refs/tags/v${HIGHFIVE_VERSION}.tar.gz
  URL_HASH SHA256=6301def8ceb9f4d7a595988612db288b448a3c0546f6c83417dab38c64994d7e
  # FIND_PACKAGE_ARGS has to be the last thing in this call because it greedily takes
  # everything after it to pass to `find_package()`
  FIND_PACKAGE_ARGS NAMES HighFive
)

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
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    set(ENABLE_LAPACK ON)
    set(BLA_VENDOR Apple)
    set(SUNDIALS_LAPACK_CASE "LOWER")
    set(SUNDIALS_LAPACK_UNDERSCORES "NONE")
elseif ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
    set(ENABLE_LAPACK ON)
    set(BLA_VENDOR OpenBLAS)
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

FetchContent_MakeAvailable(eigen yaml-cpp sundials fmt HighFive)

file(GLOB_RECURSE CT_LIB_SOURCES "*.cpp")
# The Python extension handler is built with the rest of the Python library because it
# needs to be compiled after the Cython step runs.
list(FILTER CT_LIB_SOURCES EXCLUDE REGEX "extensions/.*")
list(APPEND CT_LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/extensions/canteraShared.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/extensions/pythonShim.cpp")

if (NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
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

add_library(cantera_lib STATIC ${CT_LIB_SOURCES})
target_include_directories(cantera_lib PRIVATE "../include")
target_include_directories(cantera_lib SYSTEM PRIVATE "${Python_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}")
set_target_properties(cantera_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(cantera_lib PUBLIC
    yaml-cpp::yaml-cpp
    fmt::fmt
    Eigen3::Eigen
    HighFive
    SUNDIALS::core
    SUNDIALS::nvecserial
    SUNDIALS::cvodes
    SUNDIALS::idas
    SUNDIALS::sunmatrixband
    SUNDIALS::sunmatrixdense
    SUNDIALS::sunmatrixsparse
    SUNDIALS::sunlinsolspgmr
    SUNDIALS::sunnonlinsolnewton
    ${SUNDIALS_LIBRARIES}
)
