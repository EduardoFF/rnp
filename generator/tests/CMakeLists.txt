
add_executable(test test.cpp)
include_directories (. ../)

link_directories(. ../)

MESSAGE(STATUS "TEST ${CORELIBS}")
target_link_libraries(test ${CORELIBS} 
  ${GUROBI_LIBRARIES} 
  ${Cairomm_LIBRARIES} 
  generator solver network GetPot graph random geometry Plotter) 

ADD_TEST(test ${EXECUTABLE_OUTPUT_PATH}/test)

