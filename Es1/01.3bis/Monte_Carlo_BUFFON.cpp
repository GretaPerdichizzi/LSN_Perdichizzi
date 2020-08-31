#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <stdlib.h>
#include <iomanip>
#include "Monte_Carlo_BUFFON.h"

using namespace std;

int main() {
    
    cout << endl << " BUFFON EXPERIMENT " << endl << endl;
    cout << " Estimation of Pi Greek " << endl;
    cout << " pi = 2*L* n_throws/(n_hits * d) " << endl;
    cout << " Using L=0.25  d=0.5  " << endl << endl;
    
    SetSeed();
    
    for(int iblk=1; iblk <= nblk; ++iblk) {
        Reset(iblk);
        for(int istep=1; istep <= nstep; ++istep) {
            Throw();
        }
        Measure();
        Averages(iblk);   //Print results for current block
    }

    return 0;
}



void SetSeed(){
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
}



void Throw() {
    double L_center, x, y;
    double xtip;
    
    n_throws+=1;
    
    //simulo dove cade il centro dell'ago
    L_center = rnd.Rannyu(0,d/2);
    
    //REJECTION TECHNIQUE per simulare le coordinate x,y dell'angolo sulla circ coniometrica
    do {
        x = rnd.Rannyu();
        y = rnd.Rannyu();
        attempted += 1.0;
    } while(x*x+y*y>=1);
    
    accepted += 1.0;
    
    //valuto se l'ago cade sulla riga o no
    xtip = L_center - L/2. * x / sqrt(x*x+y*y);
    if(xtip<0) {
        n_hits+=1;
    }
}



void Measure(){  // calcolare il valore di p_greco nel blocco
    pi = 2*L* n_throws/(n_hits * d);
}



void Averages(int iblk){ //Print results for current block
    ofstream Pi;
    const int wd=15;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Pi.open("output_pi.dat",ios::app);
    glob_av  += pi;
    glob_av2 += pi*pi;
    err_p=Error(glob_av,glob_av2,iblk);
    Pi << setw(wd) << iblk <<  setw(wd) << pi << setw(wd) << glob_av/(double)iblk << setw(wd) << err_p << endl;
    Pi.close();

    cout << "----------------------------" << endl << endl;
}



void Reset(int iblk) {     //Reset block averages
    if(iblk == 1){
        glob_av = 0;
        glob_av2 = 0;
    }

    pi = 0;
    n_throws = 0;
    n_hits = 0;
    
    attempted = 0.0;
    accepted = 0.0;
}



double Error(double sum, double sum2, int iblk) {
    if(iblk==0) return 0;
    else {
        if(sum2/(double)iblk - pow(sum/(double)iblk,2)>0) return sqrt( (sum2/(double)iblk - pow(sum/(double)iblk,2)) / (double)iblk );
        else return sqrt( ( pow(sum/(double)iblk,2) - sum2/(double)iblk ) / (double)iblk );
    }
}
