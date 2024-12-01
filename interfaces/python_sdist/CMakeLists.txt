# Set minimum version 3.27 so that CMP0144, search for packages by upper-cased
# environment variables will use upper-cased variables. OLD behavior is to ignore
# upper-cased variables
cmake_minimum_required(VERSION 3.27...3.31)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

add_subdirectory("src")
add_subdirectory("cantera")

install(TARGETS _cantera DESTINATION ${SKBUILD_PLATLIB_DIR}/cantera)
