cmake_minimum_required(VERSION 3.15)
project(lammpsParser VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# include folder for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SRC
    src/lammps_parser.cpp
    src/bindings.cpp
)

pybind11_add_module(lammpsParser ${SRC_CPP})
# --- compute the platlib (site-packages) directory using the configured Python ---
# Python_EXECUTABLE is provided by find_package(Python COMPONENTS Interpreter ...)
execute_process(
  COMMAND ${Python_EXECUTABLE} -c
    "import sysconfig, json; print(sysconfig.get_path('platlib'))"
  RESULT_VARIABLE _py_exec_res
  OUTPUT_VARIABLE PYTHON_PLATLIB
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT _py_exec_res EQUAL 0)
  message(FATAL_ERROR "Failed to query Python platlib via ${Python_EXECUTABLE}")
endif()

message(STATUS "Python platlib (absolute): ${PYTHON_PLATLIB}")

# Make a path relative to the install prefix so CMake will install under the staging dir.
# Remove leading '/' if present to create a relative DESTINATION.
string(REGEX REPLACE "^/" "" PYTHON_PLATLIB_REL "${PYTHON_PLATLIB}")

message(STATUS "Python platlib (relative DESTINATION): ${PYTHON_PLATLIB_REL}")

# Install the compiled extension target into site-packages/lammpsParser (relative to staging dir)
install(TARGETS lammpsParser
        LIBRARY DESTINATION ${PYTHON_PLATLIB_REL}/lammpsParser)
# Optionally install the python package folder (if present in source tree) to ensure __init__.py is included
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser
            DESTINATION ${PYTHON_PLATLIB_REL})
endif()
#install(TARGETS lammpsParser
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/lammpsParser)

# link stdc++fs for older compilers if needed
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(lammpsParser PRIVATE stdc++fs)
  endif()
endif()
