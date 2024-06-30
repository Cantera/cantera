include(FetchContent)

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0
)
FetchContent_Declare(
  sundials
  GIT_REPOSITORY https://github.com/LLNL/sundials
  GIT_TAG v7.0.0
  INSTALL_DIR /tmp
)
FetchContent_MakeAvailable(yaml-cpp sundials)
set(CANTERA_VERSION "3.1.0a2")
set(CANTERA_SHORT_VERSION "3.1")
set(CT_SUNDIALS_VERSION 70)
set(CT_USE_SYSTEM_YAMLCPP 1)
set(FTN_TRAILING_UNDERSCORE 1)
set(LAPACK_FTN_STRING_LEN_AT_END 1)
set(LAPACK_FTN_TRAILING_UNDERSCORE 1)
set(CT_SUNDIALS_USE_LAPACK 0)
# TODO: Replace with a Python script to avoid changing the template file
configure_file("../include/cantera/base/config.h.in" "../include/cantera/base/config.h" @ONLY)

file(GLOB_RECURSE CT_LIB_SOURCES "*.cpp")
list(FILTER CT_LIB_SOURCES EXCLUDE REGEX "${CMAKE_CURRENT_SRC_DIR}/extensions/.*")

add_library(cantera_lib ${CT_LIB_SOURCES})
target_include_directories(cantera_lib PRIVATE "${CMAKE_CURRENT_SRC_DIR}" "../include" "${CMAKE_CURRENT_BINARY_DIR}/../include" "${CMAKE_CURRENT_BINARY_DIR}/../include/cantera/base")
target_link_libraries(cantera_lib PUBLIC
    yaml-cpp::yaml-cpp
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
