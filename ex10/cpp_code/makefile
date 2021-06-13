CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: ex1.exe ex2.exe

ex1.exe: ex1.o random.o salesman.o
	$(CC) ex1.o random.o salesman.o -o ex1.exe
ex1.o: ex1.cpp
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
ex2.exe: ex2.o random.o salesman.o
	mpicxx ex2.o random.o salesman.o -o ex2.exe
ex2.o: ex2.cpp
	mpicxx -c ex2.cpp -o ex2.o $(CFLAGS)
random.o: ../../tools/random.cpp ../../tools/random.h
	$(CC) -c ../../tools/random.cpp -o random.o $(CFLAGS)
salesman.o: ../../tools/salesman.cpp ../../tools/salesman.h
	$(CC) -c ../../tools/salesman.cpp -o salesman.o $(CFLAGS)

clean:
	rm *.o *.exe
