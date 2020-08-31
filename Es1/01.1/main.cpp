#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

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
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

//Saving numbers generated with Uniform Distribution in [0,1)
	ofstream WriteNumbers;
	WriteNumbers.open("numbers");
	if (WriteNumbers.is_open()){
		for(int i=0; i<1000000; i++){
		WriteNumbers << rnd.Rannyu() << endl;
		}
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	rnd.SaveSeed();

//Evaluating chi2 to check if uniform number generation [0,1) is uniform even in sub-intervals
	int M=100;
	int n=10000;	
	double chi2[M] = {0};
	WriteNumbers.open("chi2");
	if (WriteNumbers.is_open()){
		for(int j=0; j<M; j++){
			double count[M] = {0};
			for(int i=0; i<n; i++){
				double r= rnd.Rannyu();
				for(int k=0; k<M; k++){
					if(r>double(k)/M && r<double(k+1)/M)	count[k]++;
				}
			}
			for(int i=0; i<M; i++){
			chi2[j]+=pow(double(count[i]-n/M),2)/(n/M);				
			}
			WriteNumbers<<chi2[j]<<endl;
		}
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	rnd.SaveSeed();
	return 0;
}
