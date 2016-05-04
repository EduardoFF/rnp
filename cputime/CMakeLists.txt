

ADD_LIBRARY(CpuTime SHARED CpuTime.cpp)
ADD_LIBRARY(CpuTime-static CpuTime.cpp)

# hack to create both shared and static with same name
SET_TARGET_PROPERTIES(CpuTime-static PROPERTIES 
  OUTPUT_NAME CpuTime CLEAN_DIRECT_OUTPUT 1)
#ADD_LIBRARY(solver SHARED solver.cpp)


install(TARGETS CpuTime
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )

install(TARGETS CpuTime-static
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
install (FILES CpuTime.h DESTINATION include)


