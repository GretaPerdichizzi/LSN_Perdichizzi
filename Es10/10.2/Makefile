CC=mpicxx

CFLAGS= -Wall --std=c++11

main.exe: main.o
	$(CC) random.cpp Population.cpp Chromosome.cpp main.o  -o main.exe

main.o: main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

Population.o : Population.cpp Population.h
	$(CC) -c Population.cpp -o Population.o $(CFLAGS)

Chromosome.o : Chromosome.cpp Chromosome.h City.h
	$(CC) -c Chromosome.cpp -o Chromosome.o $(CFLAGS)

clean:
	rm -rf *o main.exe seed.out
	


