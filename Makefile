CXX = g++
CFLAGS = -Wall -O1 -g
LDLIBS = -lpthread
RM = rm -f

.SUFFIXES : .cpp

all : P1 P2 P3

P1 : Matrix.h Matrix.cpp P1.cpp
	$(CXX) $(CFLAGS) Matrix.cpp P1.cpp -o P1 $(LDLIBS)

P2 : Matrix.h Matrix.cpp P2.cpp
	$(CXX) $(CFLAGS) Matrix.cpp P2.cpp -o P2 $(LDLIBS)

P3 : mmreader.hpp mmreader.cpp P3.cpp
	$(CXX) $(CFLAGS) mmreader.cpp P3.cpp -o P3 $(LDLIBS) -std=c++11

test3 :
	./P3 ./matrix/2cubes_sphere.mtx 2048


.clean : 
	$(RM) P1 P2 P3
