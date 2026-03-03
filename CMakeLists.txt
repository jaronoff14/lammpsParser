cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

#pybind11_add_module(lammpsParser example.cpp)


set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/lammpsParser")
set(INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(SRCS
    "${SRC_DIR}/lammps_parser.cpp"
    "${SRC_DIR}/bindings.cpp"
)

pybind11_add_module(lammpsParser ${SRCS})
target_include_directories(lammpsParser PRIVATE ${INC_DIR})

install(TARGETS lammpsParser LIBRARY DESTINATION .)