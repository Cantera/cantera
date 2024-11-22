cmake_minimum_required(VERSION 3.1)
project (Counterflow)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

find_package(Threads REQUIRED)

include_directories("/usr/local/include")
link_directories("/usr/local/lib")

add_executable(Counterflow main.cpp)
target_link_libraries(Counterflow cantera_shared Threads::Threads)
