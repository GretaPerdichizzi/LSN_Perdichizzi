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
#include"main.h"
/***************************/

using namespace std;

int main (void){
	int parent1, parent2;
  Input(); 

	Population tt;
	tt.SetRnd(rnd);	
	tt.GetPop();
	if(!tt.CheckFunct()){	
		cerr << "The simulation was interrupted due to an error." << endl;
		return -1;
	}
	tt.Sorting();
	
	double T = 25;
	int nsteps = 150;
	int rejected = 0;

	for (int k=0; k<25; k++){	
		T/=1.1;
		nsteps += 150;
		for (int i=0; i<nsteps; i++){
			vector <City> temp1, temp2;
			double L1, L2;

			L1=tt.EvalFitness();
			parent1 = tt.Search(rnd);
			parent2 = tt.Search(rnd);	
			temp1 = tt.pop[parent1].chromo;
			temp2 = tt.pop[parent2].chromo;
			tt.Crossover(rnd,parent1,parent2);
			tt.Mutate(rnd, n_chromo-2);
			tt.Mutate(rnd, n_chromo-1);
			L2=tt.EvalFitness();
		
			if(exp(1/T*(L1-L2))<rnd.Rannyu()){
				tt.pop[n_chromo-2].chromo=temp1;
				tt.pop[n_chromo-1].chromo=temp2;
				rejected+=1;
			}
			else{}

			if(!tt.CheckFunct()){	
				cerr << "The simulation was interrupted due to an error." << endl;
				return -1;
			}
			
			tt.Sorting();
			tt.PrintFitness(k*nsteps+i);
		}
	}

	PrintConfig(tt.pop[0].chromo);	
	rnd.SaveSeed();
	cout<<"Rejected = "<<rejected<<endl;
	
return 0;
}

/***************************/

// Usefull functions 

void Input(){
	cout << "Traveling Salesman Problem        " << endl;
	cout << "Genetic Algorithm             " << endl << endl;
	int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

void PrintConfig(vector<City> c){  
	ofstream WriteConfig;
  WriteConfig.open("config.final", ofstream::app);
	int wd = 12;
	for(int i=0; i<32; i++)  WriteConfig << c[i].x << setw(wd)<< c[i].y <<endl;  
	WriteConfig << c[0].x << setw(wd)<< c[0].y <<endl;
  WriteConfig.close();
}



