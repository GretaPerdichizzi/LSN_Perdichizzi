#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "metropolis.h"
#include "random.h"

using namespace std;


void metropolis (Random rnd, double q, double xn, double yn, double zn, double &x1, double &y1, double &z1, int &accept, int &reject){
	if (q > 1){
		accept+=1;
	}
	else if (q <= 1){
		double r=0;
		r=rnd.Rannyu();
		if (r <= q){
			accept+=1;
		}
		else {
			reject+=1;
			x1=xn;
			y1=yn;
			z1=zn;		
		}
	}
}

double q_eval_100 (double xn, double yn, double zn, double x1, double y1, double z1){
	double q=0;
	return q=exp(2*(sqrt(xn*xn+yn*yn+zn*zn)-sqrt(x1*x1+y1*y1+z1*z1)));
}

double q_eval_210 (double xn, double yn, double zn, double x1, double y1, double z1){
	double q=0;
	return q=z1*z1/(zn*zn)*exp(sqrt(xn*xn+yn*yn+zn*zn)-sqrt(x1*x1+y1*y1+z1*z1));
}


