project(rnp_core)
cmake_minimum_required(VERSION 2.4)


if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)


set(COMPILE_GUROBI false)
set(COMPILE_CPLEX true)

set(COMPILE_32 FALSE)
set(COMPILE_64 TRUE)

set(CPLEX_CPP false)


find_package(GLPK)

if( GLPK_FOUND )
  include_directories(${GLPK_INCLUDE_DIRS})
  ADD_DEFINITIONS(-D__WITH_GLPK=1)
else (GLPK_FOUND )
  message( "GLPK Support Disabled")
endif( GLPK_FOUND )


IF(CMAKE_BUILD_TYPE MATCHES Debug)
  MESSAGE("DEBUG ENABLED")
  #SET(CMAKE_CXX_COMPILER '/home/eduardo/.vim/bin/cc_args.py clang++')
  #SET(CMAKE_CXX_COMPILER "/home/eduardo/.vim/bin/cc_args_clang.sh")
ELSE()
  MESSAGE("RELEASE MODE")
  #SET(CMAKE_CXX_COMPILER "clang++")
ENDIF()


if(COMPILE_GUROBI)
  add_subdirectory(wrapper )
  find_package(GUROBI REQUIRED)
  if(GUROBI_FOUND)
    include_directories(${GUROBI_INCLUDE_DIRS})
  ENDIF()
  find_package(GUROBICPP REQUIRED)
  IF(GUROBICPP_FOUND)
    include_directories(${GUROBICPP_INCLUDE_DIRS})
  ENDIF()
endif()

if(COMPILE_CPLEX)
  find_package(Cplex REQUIRED)
  IF(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIRS})
    set(CORELIBS ${CORELIBS} ${CPLEX_LIBRARIES})
  ENDIF()
endif()


find_package(CAIROMM)
IF(CAIROMM_FOUND)
  include_directories(${CAIROMM_INCLUDE_DIRS})
ENDIF()


SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin CACHE PATH

  "Single Directory for all executables"
  )

SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/libs CACHE PATH

  "Single Directory for all dynamic Libraries"
  )


#include_directories(${CMAKE_CURRENT_SRC_DIR}/generator)
add_subdirectory( solver )
add_subdirectory( cputime )
add_subdirectory( socket_lib )
add_subdirectory( generator )
add_subdirectory( multimth )





