
if(COMPILE_GUROBI)
	find_package(GUROBICPP REQUIRED)
	IF(GUROBICPP_FOUND)
		include_directories(${GUROBICPP_INCLUDE_DIRS})
	ENDIF()
ENDIF()

find_package(GLPK)

if( GLPK_FOUND )
  include_directories(${GLPK_INCLUDE_DIRS})
  ADD_DEFINITIONS(-D__WITH_GLPK=1)
  set(CORELIBS ${CORELIBS} ${GLPK_LIBRARIES})
endif( GLPK_FOUND )
  
include_directories(${CMAKE_SOURCE_DIR}/cputime)
include_directories(${CMAKE_SOURCE_DIR}/socket_lib)

set(CORELIBS ${CORELIBS} ${PTHREADS_LIBRARY} "-lm")

if(COMPILE_CPLEX)
	find_package(Cplex REQUIRED)
	IF(CPLEX_FOUND)
		include_directories(${CPLEX_INCLUDE_DIRS})
	ENDIF()
endif()

ADD_DEFINITIONS(-Wno-deprecated -Wno-write-strings -O2 -g)

IF(CPLEX_FOUND)
	set(CORELIBS ${CORELIBS} ${CPLEX_LIBRARIES})
ENDIF()

IF(GUROBICPP_FOUND)
		set(CORELIBS ${CORELIBS} ${GUROBICPP_LIBRARIES})
ENDIF()


if(COMPILE_CPLEX)
	add_definitions("-D_COMPILE_CPLEX=1")
ENDIF()

if(COMPILE_GUROBI)
	add_definitions("-D_COMPILE_GUROBI=1")
ENDIF()

if(COMPILE_64)
	add_definitions("-m64")
endif()
if(COMPILE_32)	
	add_definitions("-m32")
endif()

find_package(Pthreads REQUIRED)

find_package(Boost COMPONENTS serialization iostreams REQUIRED)

find_package(ZLIB REQUIRED)
set( CORELIBS 
	${CORELIBS} 
	${ZLIB_LIBRARIES}
	${Boost_SERIALIZATION_LIBRARY}
	${Boost_IOSTREAMS_LIBRARY}
	)


set(CORELIBS ${CORELIBS} ${PTHREADS_LIBRARIES})


#ADD_LIBRARY(solver SHARED solver.cpp)
ADD_LIBRARY(solver solver.cpp)

add_executable(tests tests.cpp)

Target_Link_Libraries(tests CpuTime)
TARGET_LINK_LIBRARIES(solver CpuTime  ${CORELIBS} ${GLOG_LIBRARIES})

#TARGET_LINK_LIBRARIES(solver CpuTime ${CPLEX_LIBRARY} ${GLPK_LIBRARY})


install (TARGETS solver DESTINATION lib)
install (FILES solver.h DESTINATION include)
install (FILES solverparameter.h DESTINATION include)


