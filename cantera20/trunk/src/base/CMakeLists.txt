cmake_minimum_required (VERSION 2.6.)

PROJECT (CanteraBase CXX)

ADD_LIBRARY (ctbase  cttypes.cpp)
ADD_EXECUTABLE(test_base test_types.cpp )
TARGET_LINK_LIBRARIES(test_base  ctbase )
