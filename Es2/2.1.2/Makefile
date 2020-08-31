CC=g++

CFLAGS= -Wall --std=c++11

main.exe: main_px.o
	$(CC) random.cpp main_px.o  -o main_px.exe

main.o: main_px.cpp
	$(CC) -c main_px.cpp -o main_px.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean:
	rm -rf *o main_px.exe seed.out
	
