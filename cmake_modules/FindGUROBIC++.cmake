
SET(GUROBI_ROOT_DIR "$ENV{GUROBI_HOME}")


MESSAGE(STATUS "Searching gurobi cpp in: ${GUROBI_ROOT_DIR}")

FIND_PATH(GRBCPP_INCLUDE_DIR
  gurobi_c++.h
  PATHS ${GUROBI_ROOT_DIR}/include
)
MESSAGE(STATUS "GUROBI CPP: INCLUDE DIR FOUND ${GRB_INCLUDE_DIR}")
FIND_LIBRARY(GRBCPP_LIBRARY
   gurobi_c++
  PATHS ${GUROBI_ROOT_DIR}/lib
)

IF (GRBCPP_INCLUDE_DIR AND GRBCPP_LIBRARY)
	SET(GUROBICPP_INCLUDE_DIRS ${GRBC++_INCLUDE_DIR})
	SET(GUROBICPP_LIBRARIES ${GRBC++_LIBRARY})
	SET(GUROBICPP_FOUND TRUE)
ENDIF (GRBCPP_INCLUDE_DIR AND GRBCPP_LIBRARY)


IF (GUROBICPP_FOUND)
	IF (NOT GUROBICPP_FIND_QUIETLY)
		MESSAGE(STATUS "Found Gurobi C++: ${GRBPP_LIBRARY}")
	ENDIF (NOT GUROBICPP_FIND_QUIETLY)
ELSE (GUROBICPP_FOUND)
	IF (GUROBICPP_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could not find Gurobi C++")
	ENDIF (GUROBICPP_FIND_REQUIRED)
ENDIF (GUROBICPP_FOUND)
