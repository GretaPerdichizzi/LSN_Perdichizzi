#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "metropolis.h"

using namespace std;
 
int main (int argc, char *argv[]){
	int M = 1000000;
	int L = 10000; //blocks lenght or metropolis steps
	int N = M/L; //blocks number
	double xn, yn, zn;
	double x1, y1, z1;
	double q;
	int accept=0, reject=0;

	//generate and initialize class Random rnd 
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
	
	//Uniform generation for T(x'|xn)
	cout<<"UNIFORM GENERATION"<<endl<<endl;
	//psi 100
	ofstream WritePsi100Unif;
	WritePsi100Unif.open("psi100_uniform.dat");
	if (WritePsi100Unif.is_open()){
		for(int j=0; j<N; j++){
			//new simulation starting from initial point
			xn=0;
			yn=0;
			zn=0;
			
			//Metropolis evaluation cycle
			for (int i=0; i<L; i++){
				x1=rnd.Rannyu(xn-1.29, xn+1.29);
				y1=rnd.Rannyu(yn-1.29, yn+1.29);
				z1=rnd.Rannyu(zn-1.29, zn+1.29);

				q=q_eval_100(xn, yn, zn, x1, y1, z1);

				metropolis(rnd, q, xn, yn, zn, x1, y1, z1, accept, reject);
				xn=x1;
				yn=y1;
				zn=z1;
				WritePsi100Unif << xn << "\t" << yn << "\t" << zn << "\t" << sqrt(xn*xn+yn*yn+zn*zn) <<endl;	
			}
		}
	}
	
	else cerr << "PROBLEM: Unable to open writing files" << endl;
	cout<<"100"<<endl;
	cout<<"accept="<<accept<<endl;
	cout<<"reject="<<reject<<endl<<endl;
	WritePsi100Unif.close();

	//psi 210
	accept=0;
	reject=0;
	ofstream WritePsi210Unif;
	WritePsi210Unif.open("psi210_uniform.dat");
	if (WritePsi210Unif.is_open()){
		for(int j=0; j<N; j++){
			//new simulation starting from initial point
			xn=1.73;
			yn=0.;
			zn=1.;
			
			//Metropolis evaluation cycle
			for (int i=0; i<L; i++){
				x1=rnd.Rannyu(xn-3, xn+3);
				y1=rnd.Rannyu(yn-3, yn+3);
				z1=rnd.Rannyu(zn-3, zn+3);

				q=q_eval_210(xn, yn, zn, x1, y1, z1);

				metropolis(rnd, q, xn, yn, zn, x1, y1, z1, accept, reject);
				xn=x1;
				yn=y1;
				zn=z1;
				WritePsi210Unif << xn << "\t" << yn << "\t" << zn << "\t" << sqrt(xn*xn+yn*yn+zn*zn) <<endl;	
			}
		}
	}
	
	else cerr << "PROBLEM: Unable to open writing files" << endl;
	cout<<"210"<<endl;
	cout<<"accept="<<accept<<endl;
	cout<<"reject="<<reject<<endl<<endl;
	WritePsi210Unif.close();

	//Gauss generation for T(x'|xn)
	cout<<"GAUSSIAN GENERATION"<<endl<<endl;
	//psi 100
	accept=0;
	reject=0;
	ofstream WritePsi100Gauss;
	WritePsi100Gauss.open("psi100_gauss.dat");
	if (WritePsi100Gauss.is_open()){
		for(int j=0; j<N; j++){
			//new simulation starting from initial point
			xn=0;
			yn=0;
			zn=0;
			
			//Metropolis evaluation cycle
			for (int i=0; i<L; i++){
				x1=rnd.Gauss(xn, 0.75);
				y1=rnd.Gauss(yn, 0.75);
				z1=rnd.Gauss(zn, 0.75);

				q=q_eval_100(xn, yn, zn, x1, y1, z1);

				metropolis(rnd, q, xn, yn, zn, x1, y1, z1, accept, reject);
				xn=x1;
				yn=y1;
				zn=z1;
				WritePsi100Gauss << xn << "\t" << yn << "\t" << zn << "\t" << sqrt(xn*xn+yn*yn+zn*zn) <<endl;	
			}
		}
	}
	
	else cerr << "PROBLEM: Unable to open writing files" << endl;
	cout<<"100"<<endl;
	cout<<"accept="<<accept<<endl;
	cout<<"reject="<<reject<<endl<<endl;
	WritePsi100Gauss.close();

	//psi 210
	accept=0;
	reject=0;
	ofstream WritePsi210Gauss;
	WritePsi210Gauss.open("psi210_gauss.dat");
	if (WritePsi210Gauss.is_open()){
		for(int j=0; j<N; j++){
			//new simulation starting from initial point
			xn=1.73;
			yn=0;
			zn=1.;
			
			//Metropolis evaluation cycle
			for (int i=0; i<L; i++){
				x1=rnd.Gauss(xn, 1.9);
				y1=rnd.Gauss(yn, 1.9);
				z1=rnd.Gauss(zn, 1.9);

				q=q_eval_210(xn, yn, zn, x1, y1, z1);

				metropolis(rnd, q, xn, yn, zn, x1, y1, z1, accept, reject);
				xn=x1;
				yn=y1;
				zn=z1;
				WritePsi210Gauss << xn << "\t" << yn << "\t" << zn << "\t" << sqrt(xn*xn+yn*yn+zn*zn) <<endl;	
			}
		}
	}
	
	else cerr << "PROBLEM: Unable to open writing files" << endl;
	cout<<"210"<<endl;
	cout<<"accept="<<accept<<endl;
	cout<<"reject="<<reject<<endl<<endl;
	WritePsi210Gauss.close();

return 0;
}
