
add_executable(main main.cpp)

include_directories(../../utils/boost-geometry)


#IF(GUROBICPP_FOUND)
#	include_directories(${GUROBICPP_INCLUDE_DIRS})
#ENDIF()
#IF(CPLEX_FOUND)
#	include_directories(${CPLEX_INCLUDE_DIRS})
#ENDIF()


link_directories(../../utils/boost-geometry/libs)

enable_testing()



add_library(metric metric.cpp)
add_library(graph graph.cpp)

#add_library(cputime CpuTime.cc)

add_library(mst mst.cpp)

add_library(network network.cpp)

#add_library(geometry Geometry.cpp)

add_library(rnpsolution rnpsolution.cc)
add_library(lpsolver lpsolver.cc)

find_package(GLPK REQUIRED)
set(CORELIBS ${CORELIBS} ${GLPK_LIBRARIES})

find_package(SOLVER REQUIRED)
set(CORELIBS ${CORELIBS} ${SOLVER_LIBRARIES})

find_package(CpuTime REQUIRED)
set(CORELIBS ${CORELIBS} ${CpuTime_LIBRARIES})

find_package(ZLIB REQUIRED)

set(CORELIBS ${CORELIBS} ${ZLIB_LIBRARIES})
find_package(SOCKETLIB REQUIRED)

target_link_libraries(lpsolver ${SOLVER_LIBRARIES}) 

add_library(matheuristic matheuristic.cc)
target_link_libraries(matheuristic ${SOCKETLIB_LIBRARIES})
add_library(random Random.cpp)
add_library(heuristic heuristic.cpp)

add_library(calc_metric calc_metric.cc)

add_library(GetPot GetPot.cc)

add_library(ga_solver ga_solver.cc)

add_library(cplex_solver cplex_solver.cc)


add_library(generator generator.cpp)
add_library(Plotter Plotter.cpp)

find_package(Pthreads REQUIRED)
if(CPLEX_CPP)
	find_package(CPLEX REQUIRED)

endif(CPLEX_CPP)

find_package(Boost COMPONENTS serialization iostreams REQUIRED)
MESSAGE("Boost LIBRARIES " ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CPLEX_INCLUDE_DIRS})
#set(CORELIBS ${GLPK_LIBRARY} ${PTHREADS_LIBRARY} "-lm")
#IF(GUROBI_FOUND)
#set(CORELIBS ${CORELIBS} ${GUROBI_LIBRARY})
#ENDIF()

#IF(GUROBICPP_FOUND)
#set(CORELIBS ${CORELIBS} ${GUROBICPP_LIBRARY})
#ENDIF()

#IF(CPLEX_FOUND)
#	set(CORELIBS ${CORELIBS} ${CPLEX_LIBRARY})
#ENDIF()

#find_package(Boost COMPONENTS filesystem)
#If(Boost_FOUND)
#	set(CORELIBS ${CORELIBS} ${Boost_LIBRARIES})
#ENDIF()

IF(CAIROMM_FOUND)
	set(CORELIBS ${CORELIBS} ${Cairomm_LIBRARY})
	ADD_DEFINITIONS("`pkg-config --libs --cflags cairomm-1.0`")
ENDIF()


find_package(GALib)

#add_definitions(-DGALIB_USE_NO_TEMPLATES=1)

#if(COMPILE_GUROBI)
#	add_definitions("-D_COMPILE_GUROBI=1")
#	set(CORELIBS ${CORELIBS} wrapper)
#endif()
#if(COMPILE_CPLEX)
#	add_definitions("-D_COMPILE_CPLEX=1")
#ENDIF()

IF( COMPILE_32)
	add_definitions("-m32")
ENDIF()

IF(COMPILE_64)
	message("compiling 64bits")
	add_definitions("-m64 -g ")
endif()

target_link_libraries(ga_solver ${GALIB_LIBRARY})

target_link_libraries(cplex_solver network ${CPLEX_LIBRARIES})

set( CORELIBS 
	${CORELIBS} 
	${Boost_LIBRARIES})
	#${Boost_IOSTREAMS_LIBRARY}
	#${Boost_SERIALIZATION_LIBRARY})
set( COMMON_LIBRARIES
	rnpsolution 
	lpsolver 
	mst metric 
	calc_metric 
	network 
	GetPot 
	graph 
	generator 
	random
# geometry depends on old boost geometry lib 
#	geometry 
	heuristic 
	Plotter 
	ga_solver 
	matheuristic 
	cplex_solver)

IF( CPLEX_CPP)
	add_library(cplex_solver_v2 cplex_solver_v2.cc)
	SET( COMMON_LIBRARIES ${COMMON_LIBRARIES} cplex_solver_v2)
	add_definitions(-DCPLEX_CPP)
ENDIF( CPLEX_CPP)

#add_subdirectory(tests)
target_link_libraries(main ${COMMON_LIBRARIES} ${CORELIBS} -lpthread)
add_executable(mytest mytest.cpp)

#link_directories(${Boost_LIBRARY_DIR})
set( TESTLIBS  
  ${Boost_LIBRARIES}
  ${SOLVER_LIBRARIES})
target_link_libraries(mytest ${TESTLIBS} )
#target_link_libraries(main ${Cairomm_LIBRARIES} solver metric calc_metric network GetPot graph wrapper generator random geometry heuristic cputime Plotter ${CORELIBS} ${GUROBI_LIBRARIES} ${GUROBICPP_LIBRARIES} -lpthread)
