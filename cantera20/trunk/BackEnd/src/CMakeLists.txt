cmake_minimum_required (VERSION 2.6)

ADD_SUBDIRECTORY(base)

ADD_SUBDIRECTORY (models)

ADD_SUBDIRECTORY (algorithms)

ADD_SUBDIRECTORY (processes)

# ADD_EXECUTABLE (testit testit.cpp )

SET (CMAKE_INSTALL_PREFIX /usr/local)

