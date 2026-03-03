cmake_minimum_required(VERSION 3.15)
project(lammpsParser VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(SRCS
    "${SRC_DIR}/lammps_parser.cpp"
    "${SRC_DIR}/bindings.cpp"
)

foreach(_f IN LISTS SRCS)
    if(NOT EXISTS "${_f}")
        message(FATAL_ERROR "Required source file not found: ${_f}")
    endif()
endforeach()

pybind11_add_module(lammpsParser ${SRCS})
target_include_directories(lammpsParser PRIVATE ${INC_DIR})
target_compile_features(lammpsParser PRIVATE cxx_std_17)

# Compatibility for older libstdc++
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(lammpsParser PRIVATE stdc++fs)
  endif()
endif()

# --- SIMPLE, RELATIVE INSTALL (robust for scikit-build) ---
install(TARGETS lammpsParser
        LIBRARY DESTINATION lammpsParser)

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser
            DESTINATION .)
endif()
