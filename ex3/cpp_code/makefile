CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex1.exe: ex1.o random.o blocking.o
	$(CC) ex1.o random.o blocking.o -o ex1.exe
ex1.o: ex1.cpp
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
random.o: ../../tools/random.cpp ../../tools/random.h
	$(CC) -c ../../tools/random.cpp -o random.o $(CFLAGS)
blocking.o: ../../tools/blocking.cpp ../../tools/blocking.h
	$(CC) -c ../../tools/blocking.cpp -o blocking.o $(CFLAGS)

clean:
	rm *.o *.exe seed.out
