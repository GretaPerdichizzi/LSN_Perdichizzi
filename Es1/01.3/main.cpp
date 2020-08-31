#include "random.h"
#include "main.h"

using namespace std;
 
int main (int argc, char *argv[]){
	double approx_p=0;
	Input();
	ofstream WritePI;
	WritePI.open("pi_eval_rej.out");
	for (int i=0; i<n; i++){
		approx_p=Simulate();
		WritePI<<approx_p<<endl;
	}
	WritePI.close();
	rnd.SaveSeed();
	return 0;
}

void Input(){
	cout << "Buffon needle problem        " << endl;
	cout << "Sampling pi             " << endl << endl;
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
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}


double Simulate(void){
	int hits=0;
	for(int i=0; i<m; i++){
		//double r = rnd.Rannyu(0,d);
		//double theta = rnd.Rannyu(0,360);
		//double z = cos(theta*180/M_PI)*L/2;

		double r=rnd.Rannyu(0,d/2);
		double x=0, y=0;
		 do{
			x = rnd.Rannyu();
			y = rnd.Rannyu();
		 } while(x*x+y*y>=1);
		double cross = r - L/2. * x / sqrt(x*x+y*y);
	  if(cross<0) {
			hits+=1;
	  }
		//if(r+z>=d || r-z<0)  hits+=1;
		else{}
	}
	//double approx_p=L/d*(double(m)/double(hits));
	double approx_p=2*L/d*(double(m)/double(hits));
	return approx_p;
}



