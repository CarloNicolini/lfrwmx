set(CMAKE_CXX_COMPILER "g++-4.8")
set(CMAKE_C_COMPILER "gcc-4.8")
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(LFRW)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(PROJECT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build/")
set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/" CACHE PATH "Single directory for all libraries")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/" CACHE PATH "Single directory for all executables")

# Check for the support of C++11 compliant compilers
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# Specify if matlab support and python wrappers are to be compiled. Default OFF
# MATLAB support needs MEX (Tested on Ubuntu 14.04, Matlab R2013a)
option(MATLAB_SUPPORT "Mex with MATLAB" OFF)
option(OCTAVE_SUPPORT "Oct with Octave mkoctfile" OFF)

if(MATLAB_SUPPORT)
    find_package(MyMatlab)
    message(STATUS "Matlab include directories: ${MATLAB_INCLUDE_DIRS}")
    include_directories(${MATLAB_INCLUDE_DIRS})
    if (MATLAB_NOT_FOUND)
        message(FATAL "Error matlab not found")
    endif()
    add_definitions("-DMATLAB_SUPPORT")
endif()

if(OCTAVE_SUPPORT)
    find_package(OCTAVE)
    include_directories(${OCTAVE_INCLUDE_DIRS})
    add_definitions("-DOCTAVE_SUPPORT")
endif(OCTAVE_SUPPORT)


# PYTHON_SUPPORT needs Cython (install it via "pip install cython")
option(PYTHON_SUPPORT "Python support through Cython extensions" OFF)
if(PYTHON_SUPPORT)
    set(Python_ADDITIONAL_VERSIONS "2.7.12")
    find_package(PythonInterp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    find_package(PythonLibs 2.7.12 EXACT REQUIRED)
    find_package(Numpy)
    include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
        # Print some informations on python libraries linked and included
        message(STATUS "Python libraries are ${PYTHON_LIBRARIES}")
        message(STATUS "Python libraries are included in ${PYTHON_INCLUDE_DIRS}")
        message(STATUS "Python versions found is ${PYTHONLIBS_VERSION_STRING}")
        include_directories("${PYTHON_INCLUDE_DIR}")
    include(UseCython)
    include_directories("${PYTHON_INCLUDE_DIR}/numpy")
    include_directories("src/lfrwmx") #necessary to include set_parameters.h in Cython
    message(STATUS "${PYTHON_INCLUDE_DIR}/numpy")
endif(PYTHON_SUPPORT)

# Add the subdirectory with the source code
add_subdirectory(src)
