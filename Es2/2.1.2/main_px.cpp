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


//Saving numbers generated with p(x) = 2 * (1-x)
	double G_p[M]={};
	ofstream WriteNumbers_p;
	WriteNumbers_p.open("p_numbers");
	if (WriteNumbers_p.is_open()){
		for(int i=0; i<M; i++){
			double a = 0;
			a = rnd.CosineTry();
			WriteNumbers_p << a << endl;
			G_p[i] = a;
		}
	} else cerr << "PROBLEM: Unable to open p_numbers" << endl;
	WriteNumbers_p.close();

	rnd.SaveSeed();



/////MONTE CARLO WITH DDP p(x) /////

//Computing I
	for(int i=0; i<M; i++){
		G_p[i] = M_PI/4*cos(M_PI*G_p[i]/2)/(1-G_p[i]);
	}

//Computing the mean over N=100
	double mu_p[N];
	for(int n=0; n<N; n++){ // number of experiments
		double sum = 0;
		for(int i=0; i<L; i++){ // number of measures in a experiment
			sum += G_p[i+n*L]; 
		}
		mu_p[n] = 1./L*sum;
	}

//Saving values for mu
	WriteNumbers_p.open("mu_p");
	if (WriteNumbers_p.is_open()){
		for(int i=0; i<N; i++){
			WriteNumbers_p << mu_p[i] << endl;
		}
	} else cerr << "PROBLEM: Unable to open mu_p" << endl;
	WriteNumbers_p.close();

	return 0;
}
