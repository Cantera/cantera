cmake_minimum_required(VERSION 3.20...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

add_subdirectory("src")
add_subdirectory("cantera")

# TODO: Disable install of Eigen includes and libsundials
install(TARGETS _cantera DESTINATION ${SKBUILD_PLATLIB_DIR}/cantera)
