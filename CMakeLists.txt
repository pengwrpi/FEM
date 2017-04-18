# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.0)

project(a1) 

set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS}")

# Let CMake know where to find custom FindFoo.cmake files
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

get_filename_component(COMPILER_DIR "${CMAKE_CXX_COMPILER}" PATH)
find_program(MPIRUN NAMES mpirun PATHS "${COMPILER_DIR}")
set(MPIRUN_PROCFLAG "-np" CACHE string
  "the command line flag to give process count to MPIRUN")
message(STATUS "MPIRUN: ${MPIRUN} ${MPIRUN_PROCFLAG}")
set(VALGRIND "" CACHE string
  "the valgrind executable")
set(VALGRIND_ARGS "" CACHE string
  "the command line arguments to VALGRIND")

set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")
find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)

#set(ENV{EIGEN3_INCLUDE_DIR} "/usr/local/include/Eigen")
include_directories(${SCOREC_PREFIX}/include)
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCES
    FEA.cpp
    shapefunction.cpp
    Analysis.cpp
    Assemble.cpp
    Ele_Load.cpp
    Ele_Stiff.cpp
    ReadEle.cpp
    Reorder.cpp
    SolveEq.cpp
    Ele_Stress.cpp
    )
add_library(fep ${SOURCES})

# this INTERFACE target bundles all the enabled libraries together
#add_library(core INTERFACE)
#target_link_libraries(core INTERFACE ${SCOREC_EXPORTED_TARGETS})

add_executable(FEA_main FEA_main.cpp) 
target_link_libraries(FEA_main PUBLIC SCOREC::core fep)
