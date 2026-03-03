cmake_minimum_required(VERSION 3.15)
project(lammpsParser VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

message(STATUS "Using Python executable: ${Python_EXECUTABLE}")
message(STATUS "pybind11 version: ${pybind11_VERSION}")

# Source layout (use CMAKE_CURRENT_SOURCE_DIR so paths are correct inside scikit-build's sandbox)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Required source files (explicit)
set(SRCS
    ${SRC_DIR}/lammps_parser.cpp
    ${SRC_DIR}/bindings.cpp
)

# Sanity: check the files exist and show errors early if not
foreach(f ${SRCS})
    if(NOT EXISTS "${f}")
        message(FATAL_ERROR "Required source file not found: ${f} (current source dir: ${CMAKE_CURRENT_SOURCE_DIR})")
    endif()
endforeach()

# Include headers
target_include_directories(${PROJECT_NAME} INTERFACE ${INC_DIR}) # placeholder until module target exists

# Create the pybind11 module
# Note: pybind11_add_module expects source files list; pass our SRCS
pybind11_add_module(lammpsParser ${SRCS})

# Now ensure include dir is attached to target
target_include_directories(lammpsParser PRIVATE ${INC_DIR})

# Platform compatibility (older libstdc++)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(lammpsParser PRIVATE stdc++fs)
  endif()
endif()

# Compute Python platlib and make a relative destination for staged install
execute_process(
  COMMAND ${Python_EXECUTABLE} -c
    "import sysconfig; print(sysconfig.get_path('platlib'))"
  OUTPUT_VARIABLE PYTHON_PLATLIB
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (NOT PYTHON_PLATLIB)
  message(FATAL_ERROR "Failed to query Python platlib path via ${Python_EXECUTABLE}")
endif()
string(REGEX REPLACE "^/" "" PYTHON_PLATLIB_REL "${PYTHON_PLATLIB}")
message(STATUS "Python platlib (relative): ${PYTHON_PLATLIB_REL}")

# Install compiled shared object into site-packages/lammpsParser inside staging dir
install(TARGETS lammpsParser
        LIBRARY DESTINATION ${PYTHON_PLATLIB_REL}/lammpsParser)

# Also ensure python package files are installed if present in source tree
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser
            DESTINATION ${PYTHON_PLATLIB_REL})
endif()
