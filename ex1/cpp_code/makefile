CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex1.exe: main_ex1.o random.o tools.o
	$(CC) main_ex1.o random.o tools.o -o ex1.exe

main_ex1.o: main_ex1.cpp
	$(CC) -c main_ex1.cpp -o main_ex1.o $(CFLAGS)

random.o: random/random.cpp
	$(CC) -c random/random.cpp -o random.o $(CFLAGS)

tools.o: tools/tools.cpp
	$(CC) -c tools/tools.cpp -o tools.o $(CFLAGS)

clean:
	rm *.o *.exe seed.out
