program = main.app
OBJS    =  ClusterProducer.o

ROOTCFLAGS    	= $(shell root-config --cflags)
ROOTLIBS      	= $(shell root-config --libs)
includes       = -I./interface

CC = g++
DEBUG = -g
CFLAGS = -c $(includes)

$(program) : $(OBJS) 
	$(CC) $(ROOTLIBS) $< -o $(program)

$(OBJS) : interface/ClusterProducer.h interface/EventGenerator.h src/ClusterProducer.cpp
	$(CC) $(CFLAGS) src/ClusterProducer.cpp $(ROOTCFLAGS)

clean:
	\rm *.o *~ p1