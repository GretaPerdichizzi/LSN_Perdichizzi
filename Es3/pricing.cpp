#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "pricing.h"
#include "random.h"

using namespace std;

Pricing :: Pricing() : Random() {};

Pricing :: ~Pricing() {};

void Pricing :: SetParams(double Tfinal){
	S0 = 100;
	K = 100;
	r = 0.1;
	sigma = 0.25;
	T=Tfinal;
	tstep = T/100;
}

double Pricing :: AssetPrice(bool discr, Random & rnd){
	double S;
	double St1=S0;
	double W;
	double N=100;
	if (discr==0){
		W = rnd.Gauss(0,1.);	
		S=S0*exp((r-0.5*pow(sigma,2))*T + sigma*W*sqrt(T));		
	return S;
	}
	else if (discr==1){
		for (int i=0; i<N; i++){
			W = rnd.Gauss(0,1);
			S=St1*exp((r-0.5*pow(sigma,2))*tstep + sigma*W*sqrt(tstep));
			St1=S;
		}	
	return S;
	}
return 0;
}

double Pricing :: CallOption(double S){
	double C=0;	
	return C=exp(-r*T)*(S-K);
}

double Pricing :: PutOption(double S){
	double P=0;	
	return P=exp(-r*T)*(K-S);
}




