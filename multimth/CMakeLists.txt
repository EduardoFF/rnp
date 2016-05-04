include_directories(${CMAKE_SOURCE_DIR}/socket_lib)
#include_directories(../generator/)
#include_directories(../../utils/boost-geometry)
add_executable(mth_manager mth_manager.cc)



#find_package(SOCKETLIB REQUIRED)
target_link_libraries(mth_manager socket_lib)
