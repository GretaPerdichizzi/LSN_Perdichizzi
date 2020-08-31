#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

int main (void){

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
	
	ofstream WriteSQ, WriteRO;
	WriteSQ.open("config.square");
	WriteRO.open("config.round");

	for (int i=0; i<32; i++){
		double x,y;	
		double theta;	
		x=rnd.Rannyu(-1, 1);
		y=rnd.Rannyu(-1, 1);
		theta=rnd.Rannyu(0,360);
		WriteSQ<<x<<setw(12)<<y<<endl;
		WriteRO<<cos(theta)<<setw(12)<<sin(theta)<<endl;
	}
	
	WriteSQ.close();
	WriteRO.close();

return 0;
}
