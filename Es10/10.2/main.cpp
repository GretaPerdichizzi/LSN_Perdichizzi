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
	
	int Ngen=1000;
	int Nmigr=200;
	int parent1, parent2;

	int size=4;
	int rank=4;
	int count = 0;

	int itag=1; int itag2=2;

	Input(size);

	Population tt;
	tt.GetPop(); 

	MPI_Init(&argc, &argv);

	MPI_Comm comm;
	comm=MPI_COMM_WORLD;
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);
	MPI_Status stat1, stat2;
	MPI_Request req;
	
	//Assigning different Primes to each rank to init rnd
	for (int i=0; i<size; i++){
		if(rank==i){
			tt.prime[0]=prime[rank];
			tt.prime[1]=prime[rank+1];
		}
	}
	tt.SetRnd(rnd);	
	if(!tt.CheckFunct()){	
		cerr << "The simulation was interrupted due to an error." << endl;
		return -1;
	}
	tt.Sorting();
	
	for (int j=0; j<Nmigr; j++){
		for (int i=0; i<Ngen; i++){			
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
			tt.PrintFitness(i+j*Ngen, size, rank);
			if(rank==0){cout<<"migr = "<<j<<" gen = "<<i<<endl;}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		int migr1=0; int migr2=0; int migr3=0; int migr4=0;
		//generating the index of the ones which will be switched
		if(rank==0){
			migr1=int(rnd.Rannyu(0,4));
			while(migr2==migr1){migr2=int(rnd.Rannyu(0,4));}
			while(migr3==migr1 || migr3==migr2){migr3=int(rnd.Rannyu(0,4));}
			while(migr4==migr1 || migr4==migr2 || migr4==migr3){migr4=int(rnd.Rannyu(0,4));}
		}
		//sending the index to all
		for(int i=1; i<size; i++){
			if(rank==0){
			MPI_Send(&migr1,1,MPI_INTEGER,i,itag, MPI_COMM_WORLD);
			MPI_Send(&migr2,1,MPI_INTEGER,i,itag, MPI_COMM_WORLD);
			MPI_Send(&migr3,1,MPI_INTEGER,i,itag, MPI_COMM_WORLD);
			MPI_Send(&migr4,1,MPI_INTEGER,i,itag, MPI_COMM_WORLD);
			}
			if(rank==i){
			MPI_Recv(&migr1,1,MPI_INTEGER,0,itag, MPI_COMM_WORLD,&stat2);
			MPI_Recv(&migr2,1,MPI_INTEGER,0,itag, MPI_COMM_WORLD,&stat2);
			MPI_Recv(&migr3,1,MPI_INTEGER,0,itag, MPI_COMM_WORLD,&stat2);
			MPI_Recv(&migr4,1,MPI_INTEGER,0,itag, MPI_COMM_WORLD,&stat2);
			}
		}

		tt.GetBest();

		//Switching best elements between couples migr1/migr2 and migr3/migr4
		if(rank==migr1){
			MPI_Isend(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr2,itag, MPI_COMM_WORLD,&req);
			MPI_Recv(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr2,itag2, MPI_COMM_WORLD,&stat2);
			
			MPI_Isend(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr2,itag, MPI_COMM_WORLD,&req);
			MPI_Recv(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr2,itag2,MPI_COMM_WORLD,&stat2);

		}
		
		if(rank==migr2){
			MPI_Send(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr1,itag2,MPI_COMM_WORLD);
			MPI_Recv(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr1,itag,MPI_COMM_WORLD, &stat1);

			MPI_Send(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr1,itag2,MPI_COMM_WORLD);
			MPI_Recv(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr1,itag,MPI_COMM_WORLD, &stat1);
		}
/*
		if(rank==migr3){
			MPI_Isend(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr4,itag, MPI_COMM_WORLD,&req);
			MPI_Recv(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr4,itag2, MPI_COMM_WORLD,&stat2);
			
			MPI_Isend(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr4,itag, MPI_COMM_WORLD,&req);
			MPI_Recv(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr4,itag2,MPI_COMM_WORLD,&stat2);

		}
		
		if(rank==migr4){
			MPI_Send(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr3,itag2,MPI_COMM_WORLD);
			MPI_Recv(&tt.best_path_x[0],32,MPI_DOUBLE_PRECISION,migr3,itag,MPI_COMM_WORLD, &stat1);

			MPI_Send(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr3,itag2,MPI_COMM_WORLD);
			MPI_Recv(&tt.best_path_y[0],32,MPI_DOUBLE_PRECISION,migr3,itag,MPI_COMM_WORLD, &stat1);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		tt.CopyBest();
		tt.Sorting();*/
	}
		
		
	PrintConfig(tt.pop[0].chromo, size, rank);	
	//rnd.SaveSeed();
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
		if(rank==0){  
			ofstream WriteConfig;
			WriteConfig.open("node0_config.final", ofstream::app);
			int wd = 12;
			for(int i=0; i<32; i++)  WriteConfig << c[i].x << setw(wd)<< c[i].y <<endl;  
			WriteConfig << c[0].x << setw(wd)<< c[0].y <<endl;
			WriteConfig.close();
		}
		if(rank==1){  
			ofstream WriteConfig;
			WriteConfig.open("node1_config.final", ofstream::app);
			int wd = 12;
			for(int i=0; i<32; i++)  WriteConfig << c[i].x << setw(wd)<< c[i].y <<endl;  
			WriteConfig << c[0].x << setw(wd)<< c[0].y <<endl;
			WriteConfig.close();
		}
		if(rank==2){  
			ofstream WriteConfig;
			WriteConfig.open("node2_config.final", ofstream::app);
			int wd = 12;
			for(int i=0; i<32; i++)  WriteConfig << c[i].x << setw(wd)<< c[i].y <<endl;  
			WriteConfig << c[0].x << setw(wd)<< c[0].y <<endl;
			WriteConfig.close();
		}
		if(rank==3){  
			ofstream WriteConfig;
			WriteConfig.open("node3_config.final", ofstream::app);
			int wd = 12;
			for(int i=0; i<32; i++)  WriteConfig << c[i].x << setw(wd)<< c[i].y <<endl;  
			WriteConfig << c[0].x << setw(wd)<< c[0].y <<endl;
			WriteConfig.close();
		}
}

