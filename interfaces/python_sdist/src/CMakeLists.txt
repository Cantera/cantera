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
  FIND_PACKAGE_ARGS NAMES SUNDIALS
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

# HighFive build options
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)

# SUNDIALS build options
set(EXAMPLES_INSTALL OFF)
set(EXAMPLES_ENABLE_C OFF)

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

execute_process(
  COMMAND ${Python_EXECUTABLE} src/configure.py "include/cantera/base/config.h.in"
  ECHO_OUTPUT_VARIABLE
  COMMAND_ECHO STDOUT
)

file(GLOB_RECURSE CT_LIB_SOURCES "*.cpp")
list(FILTER CT_LIB_SOURCES EXCLUDE REGEX "extensions/.*")
list(APPEND CT_LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/extensions/canteraShared.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/extensions/pythonShim.cpp")

add_library(cantera_lib STATIC ${CT_LIB_SOURCES})
target_include_directories(cantera_lib PRIVATE "../include" "${Python_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}")
set_target_properties(cantera_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(cantera_lib PUBLIC
    yaml-cpp::yaml-cpp
    fmt::fmt
    Eigen3::Eigen
    HighFive
    SUNDIALS::core_static
    SUNDIALS::nvecserial_static
    SUNDIALS::cvodes_static
    SUNDIALS::idas_static
    SUNDIALS::sunmatrixband_static
    SUNDIALS::sunmatrixdense_static
    SUNDIALS::sunmatrixsparse_static
    SUNDIALS::sunlinsoldense_static
    SUNDIALS::sunlinsolband_static
    SUNDIALS::sunlinsolspgmr_static
    SUNDIALS::sunnonlinsolnewton_static)
