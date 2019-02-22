# Find the Python NumPy package
# Python3_NUMPY_INCLUDE_DIR
# Python3_NUMPY_FOUND
# will be set by this script

cmake_minimum_required(VERSION 2.6)

if (Python3_EXECUTABLE)
  # Find out the include path
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
            "from __future__ import print_function\ntry: import numpy; print(numpy.get_include(), end='')\nexcept:pass\n"
            OUTPUT_VARIABLE __numpy_path)
  # And the version
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
            "from __future__ import print_function\ntry: import numpy; print(numpy.__version__, end='')\nexcept:pass\n"
    OUTPUT_VARIABLE __numpy_version)
elseif(__numpy_out)
  message(STATUS "Python3 executable not found.")
endif(Python3_EXECUTABLE)

find_path(Python3_NUMPY_INCLUDE_DIR numpy/arrayobject.h HINTS "${__numpy_path}" "${PYTHON_INCLUDE_PATH}" NO_DEFAULT_PATH)

if(Python3_NUMPY_INCLUDE_DIR)
  set(Python3_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy found")
  set(Python3_NUMPY_VERSION ${__numpy_version})
  include_directories(${Python3_NUMPY_INCLUDE_DIR})
endif(Python3_NUMPY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumPy REQUIRED_VARS Python3_NUMPY_INCLUDE_DIR VERSION_VAR __numpy_version)
