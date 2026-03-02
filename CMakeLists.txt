cmake_minimum_required(VERSION 3.14)
project(lammpsParser VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# include folder for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SRC_CPP
    src_cpp/lammps_parser.cpp
    src_cpp/bindings.cpp
)

pybind11_add_module(lammpsParser ${SRC_CPP})

# link stdc++fs for older compilers if needed
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(lammpsParser PRIVATE stdc++fs)
  endif()
endif()
