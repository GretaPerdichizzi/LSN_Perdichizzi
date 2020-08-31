CC=g++

CFLAGS= -Wall --std=c++11

main.exe: main.o random.o pricing.o
	$(CC) random.cpp pricing.cpp main.o -o main.exe

main.o: main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

pricing.o : pricing.cpp pricing.h
	$(CC) -c pricing.cpp -o pricing.o $(CFLAGS)

clean:
	rm -rf *o main.exe seed.out c_dir c_discr p_dir p_discr
	
