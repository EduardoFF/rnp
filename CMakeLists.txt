PROJECT(thesis)

cmake_minimum_required(VERSION 2.4)
set(COMPILE_GUROBI false)
set(COMPILE_CPLEX true)

set(COMPILE_32 FALSE)
set(COMPILE_64 TRUE)
set(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH})

set(CPLEX_CPP false)

find_package(GLPK REQUIRED)
include_directories(${GLPK_INCLUDE_DIRS})


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
  find_package(CPLEX REQUIRED)
  IF(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIRS})
    set(CORELIBS ${CORELIBS} ${CPLEX_LIBRARIES})
  ENDIF()
endif()

ADD_DEFINITIONS(-Wno-deprecated -Wno-write-strings -O3)

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

add_subdirectory( generator )
add_subdirectory( multimth )





