#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){
	int M = 1000000; // total number of throws
	int N = 100; // total number of blocks
	int L = M/N; // throws in each block

	Random rnd;
	int seed[4];
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
	}	else cerr << "PROBLEM: Unable to open seed.in" << endl;

//Saving numbers generated with Uniform Distribution in [0,1)
	double G[M]={};
	ofstream WriteNumbers;
	WriteNumbers.open("uniform_numbers");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double a = 0;
			a = rnd.Rannyu();
			WriteNumbers << a << endl;
			G[i] = a;
		}
	} else cerr << "PROBLEM: Unable to open uniform_numbers" << endl;
	WriteNumbers.close();

	rnd.SaveSeed();


/////MONTE CARLO WITH UNIFORM DDP IN [0,1]/////

//Computing G(x_i) = pi/2 * cos(pi*x_i/2)
	for(int i=0; i<M; i++){
		G[i] = M_PI*0.5*cos(M_PI*0.5*G[i]);
	}

//Computing the mean over N=100
	double mu[N];
	for(int n=0; n<N; n++){ // number of experiments
		double sum = 0;
		//double sum2 = 0;
		for(int i=0; i<L; i++){ // number of measures in a experiment
			sum += G[i+n*L]; 
			//sum2 += pow( G[i+n*L], 2 );
		}
		mu[n] = 1./L*sum;
	}

//Saving values for mu
	WriteNumbers.open("mu");
	if (WriteNumbers.is_open()){
		for(int i=0; i<N; i++){
			WriteNumbers << mu[i] << endl;
		}
	} else cerr << "PROBLEM: Unable to open mu" << endl;
	WriteNumbers.close();

	return 0;
}
