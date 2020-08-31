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

void Input(void);
void PrintConfig(vector<City> c);

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
	
	for (int i=0; i<10000; i++){	
		//cout<<i<<endl;		
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
		tt.PrintFitness(i);
	}
	PrintConfig(tt.pop[0].chromo);	
	rnd.SaveSeed();
	
return 0;
}

/***************************/

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

