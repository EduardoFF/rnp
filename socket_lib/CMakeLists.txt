ADD_LIBRARY(socket_lib socket_lib.cc)

install (TARGETS socket_lib DESTINATION lib)
install (FILES socket_lib.h DESTINATION include)



