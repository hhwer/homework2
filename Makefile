
test: main.cc Mymat.cc Mymat.h
	mpic++ -o test main.cc Mymat.h Mymat.cc -g -Wall -std=c++14


clean:
	rm *.o test

run: 
	mpirun -np ./test 1


