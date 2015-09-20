

CFLAGS =  -Wno-deprecated -O2 

GUROBI_PATH = /opt/gurobi/gurobi300

INCLUDE =  -I$(GUROBI_PATH)/linux32/include 

LIB = -L$(GUROBI_PATH)/linux32/lib -lpthread -lgurobi30 -lglpk

CPPFILE = generator.cpp

OUT = generator

OBJ = solver.o network.o graph.o GetPot.o 

WRAPPER = ../wrapper/wrapper.o


all: generator



GetPot.o: GetPot.cc
	g++ -Wall -O2 -c GetPot.cc -o GetPot.o

solver.o: solver.cpp solver.h
	g++ $(CFLAGS) -c  $(INCLUDE) $(LIB) solver.cpp -o solver.o

network.o: network.cpp network.h
	g++ $(CFLAGS) -c $(INCLUDE) $(LIB) network.cpp -o network.o

graph.o: graph.cpp graph.h
	g++ $(CFLAGS) -c $(INCLUDE) $(LIB) graph.cpp -o graph.o

#wrapper: ../wrapper/wrapper.cc
#	gcc ../wrapper/wrapper.c -lglpk -L/opt/gurobi/gurobi300/linux32/lib -I/opt/gurobi/gurobi300/linux32/include -lgurobi30 -lpthread -std=c99 -o ../wrapper/wrapper


generator: generator.cpp network.o GetPot.o solver.o graph.o
	g++ $(CFLAGS) $(INCLUDE) $(LIB) $(CPPFILE) -o $(OUT) $(OBJ) $(WRAPPER)

clean: 
	rm -f $(OBJ) $(OUT)

