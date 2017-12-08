
main: test
	mpirun -np 8 ./test 1
test: main.o Mymat.o FMymat.o 
	mpic++ -o test main.o FMymat.o Mymat.o Mymat.h -g -Wall -std=c++14
main.o: main.cc
	mpic++ -c main.cc -std=c++14
Mymat.o: Mymat.cc
	mpic++ -c Mymat.cc -std=c++14
FMymat.o: FMymat.cc
	mpic++ -c FMymat.cc -std=c++14

clean:
	rm *.o test
gdb: main.cc Mymat.cc FMymat.cc Mymat.h
	mpic++ -o test main.cc FMymat.cc Mymat.cc Mymat.h -g -Wall -std=c++14
run: 
	mpirun -np 8 ./test 1


