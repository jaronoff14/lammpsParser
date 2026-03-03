cmake_minimum_required(VERSION 3.15)
project(lammpsParser VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

message(STATUS "Using Python executable: ${Python_EXECUTABLE}")
message(STATUS "pybind11 version: ${pybind11_VERSION}")
message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")

# Resolve source and include dirs explicitly (absolute paths)
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
message(STATUS "SRC_DIR = ${SRC_DIR}")
message(STATUS "INC_DIR = ${INC_DIR}")

# Candidate source files (adjust/add if you name differently)
set(SRCS
    "${SRC_DIR}/lammps_parser.cpp"
    "${SRC_DIR}/bindings.cpp"
)

# Check for existence of sources and fail with clear message if missing
foreach(_f IN LISTS SRCS)
    if(EXISTS "${_f}")
        message(STATUS "  FOUND: ${_f}")
    else()
        message(FATAL_ERROR "Required source file not found: ${_f}\nMake sure the file exists and is committed to the repo.")
    endif()
endforeach()

# Create the pybind11 extension target (this must come before target_* calls)
pybind11_add_module(lammpsParser ${SRCS})
message(STATUS "Created pybind11 module target 'lammpsParser'")

# Now attach include dirs and basic compile features
target_include_directories(lammpsParser PRIVATE ${INC_DIR})
target_compile_features(lammpsParser PRIVATE cxx_std_17)

# Compatibility for older libstdc++
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(lammpsParser PRIVATE stdc++fs)
  endif()
endif()

# Compute Python platlib path and strip leading slash so install DESTINATION is relative
execute_process(
  COMMAND ${Python_EXECUTABLE} -c
    "import sysconfig; print(sysconfig.get_path('platlib'))"
  OUTPUT_VARIABLE PYTHON_PLATLIB
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (NOT PYTHON_PLATLIB)
  message(FATAL_ERROR "Failed to query Python platlib via ${Python_EXECUTABLE}")
endif()
string(REGEX REPLACE "^/" "" PYTHON_PLATLIB_REL "${PYTHON_PLATLIB}")
message(STATUS "Python platlib (relative): ${PYTHON_PLATLIB_REL}")

# Install the compiled module into staging_dir/<platlib>/lammpsParser so it will be
# included inside the wheel as lammpsParser/<so>.
install(TARGETS lammpsParser
        LIBRARY DESTINATION ${PYTHON_PLATLIB_REL}/lammpsParser)

# Also ensure python package files (shim __init__.py) are included in the wheel
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lammpsParser
            DESTINATION ${PYTHON_PLATLIB_REL})
endif()
