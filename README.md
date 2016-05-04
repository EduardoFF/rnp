# rnp


Dependencies:

++++++++++ CPLEX ++++++++++

tested with cplex studio 12.6

++++++++++ GAlib ++++++++++


apt-get install libga-dev

alternatively, you can install it manually:

extract dependencies/galib247.zip and go to that folder from terminal.

edit the file 'makevars'. This file contains the compiler flags that
are unique for each cpu/os/compiler, as well as the destination
directory for doing 'make install'.  add -fpermissive flag to avoid
compilation errors

To build the library and the non-graphic examples:
% make
To install GAlib on your system:
% make install


++++++++++ Other ++++++++++
sudo apt-get install cmake build-essential g++-4.8 libstdc++-4.8-dev libboost-dev libboost-iostreams-dev libboost-serialization-dev zlib1g-dev libgoogle-glog-dev libglpk-dev graphviz


++++++++++ Building ++++++++++


cd rnp_core
mkdir build
cd build
cmake -DCPLEX_ROOT_DIR=<your CPLEX path> ../