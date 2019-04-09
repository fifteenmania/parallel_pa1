CXX = g++
CFLAGS = -Wall -O3 -g
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

test : test1 test2 test3 test4 test5 test6

test1 :
	./P3 ./matrix/2cubes_sphere.mtx 2048
test2 : 
	./P3 ./matrix/cage12.mtx 1024
test3 : 
	./P3 ./matrix/consph.mtx 2048
test4 :
	./P3 ./matrix/cop20k_A.mtx 2048
test5 : 
	./P3 ./matrix/filter3D.mtx 2048
test6 :
	./P3 ./matrix/hood.mtx 1024
test7 :  
	./sparse ./matrix/m133-b3.mtx 1024
test12 :
	./sparse ./matrix/mac_econ_fwd500.mtx 1024
test13 :
	./sparse ./matrix/majorbasis.mtx 1024
test14 :
	./sparse ./matrix/mario002.mtx 512
test15 :
	./sparse ./matrix/mc2depi.mtx 512
test16 :
	./sparse ./matrix/offshore.mtx 1024
test17 :
	./sparse ./matrix/patents_main.mtx 1024
test18 :
	./sparse ./matrix/pdb1HYS.mtx 4096
test19 :
	./sparse ./matrix/poisson3Da.mtx 16384
test20 :
	./sparse ./matrix/pwtk.mtx 1024
test21 :
	./sparse ./matrix/rma10.mtx 4096
test22 :
	./sparse ./matrix/scircuit.mtx 1024
test23 :
	./sparse ./matrix/shipsec1.mtx 1024
test24 : 
	./sparse ./matrix/webbase-1M.mtx 256


clean : 
	$(RM) P1 P2 P3
