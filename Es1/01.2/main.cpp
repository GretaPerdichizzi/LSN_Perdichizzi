#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){
	int M=10000;
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

//Uniform dice
	ofstream WriteNumbers;
	WriteNumbers.open("S1_uniform");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<1; j++){
				double r=rnd.Rannyu();
				sum+=r;
			}
			WriteNumbers << sum << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S2_uniform");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<2; j++){
				double r=rnd.Rannyu();
				sum+=r;
			}
			WriteNumbers << sum/2 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S10_uniform");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<10; j++){
				double r=rnd.Rannyu();
				sum+=r;
			}
			WriteNumbers << sum/10 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S100_uniform");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<100; j++){
				double r=rnd.Rannyu();
				sum+=r;
			}
			WriteNumbers << sum/100 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	rnd.SaveSeed();

//Exponential dice
	WriteNumbers.open("S1_exp");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<1; j++){
				double r=rnd.Exponential(1);
				sum+=r;
			}
			WriteNumbers << sum/2 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S2_exp");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<2; j++){
				double r=rnd.Exponential(1);
				sum+=r;
			}
			WriteNumbers << sum/2 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S10_exp");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<10; j++){
				double r=rnd.Exponential(1);
				sum+=r;
			}
			WriteNumbers << sum/10 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S100_exp");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<100; j++){
				double r=rnd.Exponential(1);
				sum+=r;
			}
			WriteNumbers << sum/100 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	rnd.SaveSeed();

//Gaussian dice
	WriteNumbers.open("S1_gauss");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<1; j++){
				double r=rnd.Lorentz(1.,0);
				sum+=r;
			}
			WriteNumbers << sum << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S2_gauss");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<2; j++){
				double r=rnd.Lorentz(1., 0);
				sum+=r;
			}
			WriteNumbers << sum/2 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S10_gauss");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<10; j++){
				double r=rnd.Lorentz(1.,0);
				sum+=r;
			}
			WriteNumbers << sum/10 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	WriteNumbers.open("S100_gauss");
	if (WriteNumbers.is_open()){
		for(int i=0; i<M; i++){
			double sum=0;
			for(int j=0; j<100; j++){
				double r=rnd.Lorentz(1.,0);
				sum+=r;
			}
			WriteNumbers << sum/100 << endl;
		}	
	} else cerr << "PROBLEM: Unable to open numbers" << endl;
	WriteNumbers.close();

	rnd.SaveSeed();
	rnd.SaveSeed();
	return 0;
}
