#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include "random.h"
#include "Population.h"
#include "Chromosome.h"
#include "main.h"
#include "mpi.h"
/***************************/

using namespace std;

int main (int argc, char* argv[]){
	int Ngen=200;
	int Nmigr=100;
	int parent1, parent2;

	int size=4;
	int rank=4;
	int count = 0;

	Input(size);

	Population tt;
	tt.GetPop(); 

	
	MPI_Init(&argc, &argv);

	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	//Assigning different Primes to each rank	to init rnd
	for (int i=0; i<size; i++){
		if(rank==i){
			tt.prime[0]=prime[count];
			tt.prime[1]=prime[count+1];
		}
		count+=2;
	}

	tt.SetRnd(rnd);	
	parent1 = tt.Search(rnd);

	
	if(!tt.CheckFunct()){	
		cerr << "The simulation was interrupted due to an error." << endl;
		return -1;
	}
	
	tt.Sorting();
	for (int j=0; j<Ngen; j++){
		for (int i=0; i<Nmigr; i++){	
			//cout<<"Ciao"<<endl;		
			parent1 = tt.Search(rnd);
			parent2 = tt.Search(rnd);
			tt.Crossover(rnd,parent1,parent2);
			tt.Mutate(rnd, 98);
			tt.Mutate(rnd, 99);	
			if(!tt.CheckFunct()){	
				cerr << "The simulation was interrupted due to an error." << endl;
				return -1;
			}
			tt.Sorting();
			tt.PrintFitness(i, size, rank);
		}
		//Exchanging best path random between processes
		for (int i=0; i<size; i++){
			if(rank==i){
				int rnd_rank=int(rnd.Rannyu(0,4));
				for(int k=0; k<=31; k++){
					MPI_Gather(&tt.pop[0].chromo[k], 2, MPI_INTEGER, &tt.pop[0].chromo[k], 2, MPI_INTEGER, rnd_rank, comm);
				}
			}
		}
	}
	PrintConfig(tt.pop[0].chromo, size, rank);	
	rnd.SaveSeed();
	MPI_Finalize();
return 0;
}

/***************************/
//Initialize
void Input(int size){
	//cout << "Traveling Salesman Problem        " << endl;
	//cout << "Genetic Algorithm             " << endl << endl;
	int count = 0;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		for (int i=0; i<size; i++){
			Primes >> prime[count] >> prime[count+1];
			count+=2;
		} 
	}
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	
	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
	input >> property;
	if( property == "RANDOMSEED" ){
		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	//So I saved 4 couples of Primes and the 4 seeds (equal for all ranks)
}

void PrintConfig(vector<City> c, int size, int rank){
	for (int j=0; j<size; j++){
		if(rank==j){  
			ofstream WriteConfig;
			WriteConfig.open("node"+j+"config.final", ofstream::app);
			int wd = 12;
			for(int i=0; i<32; i++)  WriteConfig << c[i].x << setw(wd)<< c[i].y <<endl;  
			WriteConfig << c[0].x << setw(wd)<< c[0].y <<endl;
			WriteConfig.close();
		}
	}
}

