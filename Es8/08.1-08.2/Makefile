CC = g++
CFLAGS = -Wall -O3 --std=c++11

QMC_1d_free_particle.exe : QMC_1d_free_particle.o random.o
	$(CC) random.o QMC_1d_free_particle.o -o QMC_1d_free_particle.exe
QMC_1d_free_particle.o : QMC_1d_free_particle.cpp
	$(CC) -c QMC_1d_free_particle.cpp -o QMC_1d_free_particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o QMC_1d_free_particle.exe
