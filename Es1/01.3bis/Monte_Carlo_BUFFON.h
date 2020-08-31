/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __BUFFON__
#define __BUFFON__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const double L=0.25;
const double d=0.8;
int n_hits,n_throws;
double pi;

// averages
double blk_av,accepted,attempted;
double glob_av,glob_av2;
double err_p;

// simulation
int nstep=100000;
int nblk=100;

//functions
void SetSeed(void);
void Throw(void);
void Measure(void);
void Averages(int);
void Reset(int);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
