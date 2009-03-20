cmake_minimum_required (VERSION 2.6)

PROJECT(Cantera)

INCLUDE( ${CMAKE_ROOT}/Modules/FindPythonInterp.cmake)
INCLUDE( ${CMAKE_ROOT}/Modules/FindPythonLibs.cmake)

INCLUDE( ${CMAKE_ROOT}/Modules/FindMatlab.cmake)


ADD_SUBDIRECTORY (src)

