#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "QMC_1d_free_particle.h"

using namespace std;

int main() {
    Input();
		H_mean(x0, mu, sigma);
    return 0;
}



void Input(void) {
    ifstream ReadInput;
    cout << "Free particle in 1D         " << endl;
    cout << "Potential V(x) = x^4 - 5/2 x^2    " << endl << endl;
    cout << "Hamiltonian H = T + V(x)" << endl << endl;
    cout << "Trial wave function Psi(x) = exp(- (x - mean)^2 / sigma^2 + exp(- (x + mean)^2 / sigma^2" << endl << endl;
		cout << "------------------------------------------------" <<endl<<endl;
		cout << "Simulated annealing to find optimal mu and sigma       " << endl << endl;

    //Initializing random numbers
		int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
}

//Evaluating trial wave function psi
double Psi(double x, double mu, double sigma) {
    double psi;
    psi=exp(-pow((x-mu),2)/(2.*pow(sigma,2)))+exp(-pow((x+mu),2)/(2.*pow(sigma,2)));
    return psi;
}

//Evaluating square trial wave function psi2
double Psi2(double x, double mu, double sigma) {
    double psi;
    psi=exp(-pow((x-mu),2)/(2.*pow(sigma,2)))+exp(-pow((x+mu),2)/(2.*pow(sigma,2)));
    return pow(psi,2);
}


//Evaluating second derivate of psi
double Der2Psi(double x, double mu, double sigma) {
    double der2psi;
		double exp1, exp2;
		exp1 = exp(-pow((x-mu),2)/2./pow(sigma,2));
		exp2 = exp(-pow((x+mu),2)/2./pow(sigma,2));
		der2psi = (-exp1-exp2+(exp1*(pow((x-mu),2))/pow(sigma,2))+(exp2*(pow((x+mu),2))/pow(sigma,2)))/pow(sigma,2);
		return der2psi;
}

double DerPsi2(double x, double mu, double sigma){
	double derpsi2;
	double exp1, exp2;
	exp1 = exp(-pow((x-mu),2)/2./pow(sigma,2));
	exp2 = exp(-pow((x+mu),2)/2./pow(sigma,2));
	derpsi2 = 2*(exp1+exp2)*(-(x-mu)*exp1/pow(sigma,2)-(x+mu)*exp2/pow(sigma,2));
	return derpsi2;
}


//Evaluating potential energy
double PotentialEnergy(double x) {
    return pow(x,2)*(pow(x,2)-5./2.);
}

//Evaluating local energy
double LocalEnergy(double x, double mu, double sigma){
	double der2psi, psi, K, V;
	der2psi = Der2Psi (x, mu, sigma);
	psi = Psi(x, mu, sigma);
	K = -0.5*der2psi/psi;
	V = PotentialEnergy(x);

	return (K+V);
}


//Evaluating mean of H
void H_mean(double x0, double mu, double sigma){
		double x=x0;
		double blk_av=0, blk_av2=0, err_H=0;
		
		ofstream WriteH, Writex;
		WriteH.open("output.H.0",ios::app);
		Writex.open("output.x.0",ios::app);		
		for (int j=1; j<=Nblocks; j++){
			accept = 0;
			reject = 0;
			double E = LocalEnergy(x0, mu, sigma);
			for (int i=0; i<Nsteps; i++){
				double x1=rnd.Rannyu(x+2.5, x-2.5);		
				double q=Psi2(x1, mu, sigma)/Psi2(x, mu, sigma);
				double r=rnd.Rannyu();
				
				if (r <= q){
					accept+=1;
					x=x1;
					Writex << x << endl;
				}
				else reject+=1;	

				E += LocalEnergy(x, mu, sigma);
			}
			blk_av += E/double(Nsteps);
    	blk_av2 += pow(E/double(Nsteps),2);
	    err_H=Error(blk_av,blk_av2,j);
			WriteH << setw(wd) << j <<  setw(wd) << blk_av/double(j) << setw(wd) << err_H << setw(wd) << mu << setw(wd) << sigma << endl;
		}	
		WriteH.close();
		Writex.close();
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

