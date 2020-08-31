#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
	int step = 100; // total number of steps in each walk
	int L = 10000; //simulation repetitions
	
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

	//********DISCRETE SAMPLING**********//

	int a=1;				    // lenght of the step
	ofstream WriteRN;
	WriteRN.open("discrete_lenghts");
	double RW [step]={0};
	double RW2 [step]={0};
	if (WriteRN.is_open()){
		for(int n=0; n<L; n++){
			int count=0;
			int walk[] = {0, 0, 0};   // the origin in 0,0,0
			for (int i=n*step; i<step*(n+1); i++){
				//Selecting the direction and the orientation of the walk
				double r = rnd.Rannyu();
				for (int dir=0; dir<3; dir++){
					double temp = dir;
					if(r>=temp/6 && r<=(temp+1)/6){
						walk[dir]=walk[dir]+a;
					}
					if(r>=(temp+3)/6 && r<=(temp+4)/6){
						walk[dir]=walk[dir]-a;
					}
				}
				RW[count]+=sqrt(pow(walk[0],2)+pow(walk[1],2)+pow(walk[2],2));
				RW2[count]+=pow(walk[0],2)+pow(walk[1],2)+pow(walk[2],2);
				count++;
			}	
		}
		for(int i=0; i<step; i++){
			WriteRN << i <<"\t" << RW[i]/L << "\t" << sqrt(RW2[i]/L/L-pow(RW[i]/L,2)/L) << endl;
		}
	}

	else cerr << "PROBLEM: Unable to open discrete_lenghts" << endl;

	WriteRN.close();
	rnd.SaveSeed();	

	//********CONTINUUM SAMPLING**********//

	ofstream WriteRW;
	WriteRW.open("continuum_lenghts");
	double RW_continuum [step]={0};
	double RW2_continuum [step]={0};
	if (WriteRW.is_open()){
		for(int n=0; n<L; n++){
			int count=0;
			double walk[] = {0, 0, 0};   // the origin in 0,0,0
			for (int i=n*step; i<step*(n+1); i++){
				//Selecting the direction and the orientation of the walk
				double phi = rnd.Rannyu(0,2*M_PI);
				double theta = rnd.Theta3D();
				walk[0]+=sin(theta)*cos(phi);
				walk[1]+=sin(theta)*sin(phi);
				walk[2]+=cos(theta);
				RW_continuum[count]+=sqrt(pow(walk[0],2)+pow(walk[1],2)+pow(walk[2],2));
				RW2_continuum[count]+=pow(walk[0],2)+pow(walk[1],2)+pow(walk[2],2);
				count++;
			}	
		}
		for(int i=0; i<step; i++){
			WriteRW << i <<"\t" << RW_continuum[i]/L << "\t" << sqrt(RW2_continuum[i]/L/L-pow(RW_continuum[i]/L,2)/L) << endl;
		}
	}

	else cerr << "PROBLEM: Unable to open continuum_lenghts" << endl;

	WriteRW.close();
	rnd.SaveSeed();	


return 0;
}
