#ifndef __1D_FREE_PARTICLE__
#define __1D_FREE_PARTICLE__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

double x0 = 0;
double mu = 0.81;
double sigma = 0.65;
int accept, reject;

const int Nsteps=1000; //step metropolis to evaluate H_mean into the simulated annealing
const int Nblocks=100;

int wd = 12;

//functions
void Input(void);
double Psi(double, double, double);
double Psi2(double, double, double);
double DerPsi2(double x, double mu, double sigma);
double Der2Psi(double, double, double);
double PotentialEnergy(double);
double LocalEnergy(double, double, double);
void H_mean(double x0, double mu, double sigma);
double Error(double sum, double sum2, int iblk);
#endif
