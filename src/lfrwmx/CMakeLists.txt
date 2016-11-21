cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(LFRW)

set(LFRW_SRCS
benchm.cpp
binary_benchm.cpp
cast.cpp
cc.cpp
combinatorics.cpp
histograms.cpp
print.cpp
random.cpp
set_parameters.cpp
)

set(LFRW_HDRS
benchm.h
binary_benchm.h
cast.h
cc.h
combinatorics.h
FileLogger.h
histograms.h
print.h
random.h
set_parameters.h
standard_include.h
)


add_library(LFRW SHARED ${LFRW_SRCS} ${LFRW_HDRS})

add_executable(run_benchm run_benchm.cpp)
add_executable(test_benchm test_benchm.cpp)
target_link_libraries(run_benchm LFRW)
target_link_libraries(test_benchm LFRW)

add_executable(test_powerlaw test_powerlaw.cpp)
target_link_libraries(test_powerlaw LFRW)

if(MATLAB_SUPPORT)
    add_mex(lfrw_mx lfr.cpp)
    target_link_libraries(lfrw_mx LFRW ${MATLAB_LIBRARIES})
    add_mex(powerlaw_comm_sizes powerlaw_comm_sizes.cpp)
    target_link_libraries(powerlaw_comm_sizes LFRW ${MATLAB_LIBRARIES})
endif(MATLAB_SUPPORT)

if(OCTAVE_SUPPORT)
    octave_add_oct(lfrw_oct SOURCES lfr.cpp LINK_LIBRARIES LFRW EXTENSION "mex")
    #target_link_libraries(lfrw LFRW ${MATLAB_LIBRARIES})
endif(OCTAVE_SUPPORT)


if(PYTHON_SUPPORT)
    set_source_files_properties(pylfrw.pyx PROPERTIES CYTHON_IS_CXX 1 )
    cython_add_module(pylfrw pylfrw.pyx)
    target_link_libraries(pylfrw LFRW)
endif(PYTHON_SUPPORT)
