add_library(wrapper wrapper.c)

ADD_DEFINITIONS(-std=c99)

include_directories(${GUROBI_INCLUDE_DIRS})

#include_directories(/opt/gurobi/gurobi300/linux32/include)
#target_link_libraries(wrapper ${CORELIBS})

