include(FetchContent)

FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen
  GIT_TAG 3.4.0
)
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0
)
FetchContent_Declare(
  sundials
  GIT_REPOSITORY https://github.com/LLNL/sundials
  GIT_TAG v7.0.0
)
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt
  GIT_TAG 10.2.1
)
FetchContent_Declare(
  HighFive
  GIT_REPOSITORY https://github.com/BlueBrain/HighFive
  GIT_TAG v2.9.0
)

find_package(Boost 1.70 CONFIG REQUIRED)

set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(EXAMPLES_INSTALL OFF)
set(FMT_INSTALL OFF)
set(FMT_MASTER_PROJECT OFF)
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
target_include_directories(cantera_lib PRIVATE "../include" "${Python_INCLUDE_DIRS}")
set_target_properties(cantera_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(cantera_lib PUBLIC
    yaml-cpp::yaml-cpp
    fmt::fmt
    Eigen3::Eigen
    Boost::boost
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
