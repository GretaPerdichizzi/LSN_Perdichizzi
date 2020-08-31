#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "pricing.h"

using namespace std;
 
int main (int argc, char *argv[]){
	int M = 1000000;
	Random rnd;
	int seed[4];
	int p1, p2;
	double S_dir;
	double S_discr;

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

	Pricing price;

	price.SetParams(1.);
	
	ofstream WriteCallOption_dir;
	WriteCallOption_dir.open("c_dir");

	ofstream WritePutOption_dir;
	WritePutOption_dir.open("p_dir");

	ofstream WriteCallOption_discr;
	WriteCallOption_discr.open("c_discr");

	ofstream WritePutOption_discr;
	WritePutOption_discr.open("p_discr");

	if (WriteCallOption_dir.is_open() && WritePutOption_dir.is_open() && WriteCallOption_discr.is_open() && WritePutOption_discr.is_open()){
		for (int i=0; i<M; i++){
			//direct
			S_dir = price.AssetPrice(0, rnd);			
			if((S_dir-price.K)>0){
				WriteCallOption_dir << price.CallOption(S_dir) << "\t";
				WritePutOption_dir << 0 << "\t";
			}
			else{
				WriteCallOption_dir << 0 << "\t";
				WritePutOption_dir << price.PutOption(S_dir) << "\t";	
			}
			//discretized
			S_discr = price.AssetPrice(1, rnd);
			if((S_discr-price.K)>0){
				WriteCallOption_discr << price.CallOption(S_discr) << "\t";
				WritePutOption_discr << 0 << "\t";
			}
			else{
				WriteCallOption_discr << 0 << "\t";
				WritePutOption_discr << price.PutOption(S_discr) << "\t";	
			}
		}
	}

	else cerr << "PROBLEM: Unable to open writing files" << endl;
	WriteCallOption_dir.close();
	WritePutOption_dir.close();
	WriteCallOption_discr.close();
	WritePutOption_discr.close();

return 0;
}
