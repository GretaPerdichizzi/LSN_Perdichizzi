/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
	for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
		Reset(iblk);   //Reset block averages
		for(int istep=1; istep <= nstep; ++istep){
			Move();           //Move particles with Verlet algorithm
			if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
			if(istep%10 == 0){
					//cout<<"istep"<<istep<<endl;
		      Measure();     //Properties measurement
					Accumulate(); //Update block averages
		      //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
		      nconf += 1;
		   }
			if (iblk==nblk && istep == nstep-1) ConfOld();
		}
		Averages(iblk);   //Print results for current block
	}
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
	ReadInput >> nblk;
	ReadInput >> restart;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
	if (restart==0 )	cout << "This code starts with random r(t-dt) " << endl << endl;
	else if (restart==1 )	cout << "This code starts with r(t-dt) of the previous simulation " << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
	iw = 4; //Pressure
  n_props = 5; //Number of observables

	//measurement of g(r)
  igofr = 5;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;


//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
	
	if (restart == 1) {
		cout << "Read initial old configuration from file old.0 " << endl << endl;
		ReadConf.open("old.0");
		for (int i=0; i<npart; ++i){
		  ReadConf >> xold[i] >> yold[i] >> zold[i];
		  xold[i] *= box;
		  yold[i] *= box;
		  zold[i] *= box;
  }
  ReadConf.close();
	}

	else{}

//Prepare initial velocities
	double sumv[3] = {0.0, 0.0, 0.0};
	double sumv2 = 0.0, fs;
	if( restart == 0 ){
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		for (int i=0; i<npart; ++i){
			vx[i] = rand()/double(RAND_MAX) - 0.5;
			vy[i] = rand()/double(RAND_MAX) - 0.5;
			vz[i] = rand()/double(RAND_MAX) - 0.5;

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		sumv2 = 0.0;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
	}
	if( restart == 1 ){
		cout << "Loading previous configurations and doing a 1st Verlet step, calculate velocities " << endl << endl;
		Move();
		sumv2 = 0.0;
		for (int i=0; i<npart; ++i){
			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
	}
	sumv2 /= (double)npart;

 	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;

		xold[i] = Pbc(x[i] - vx[i] * delta);
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
  }
	Measure();
	//Print initial values for the potential energy and virial
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Virial                   = " << walker[iw]/(double)npart << endl;
  cout << "Pressure                 = " << rho * temp + (walker[iw]) / vol << endl << endl;

	return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v = 0.0, w = 0.0, t=0.0;
  double vij, wij;
  double dx, dy, dz, dr;

	/* Uncomment to print instant values
  ofstream Epot, Ekin, Etot, Temp;
  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
	*/
	
	//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;
	
	//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
		//update of the histogram of g(r)
		bin=(int)floorf(dr/bin_size);
		//cout<<"bin = "<<bin<<endl;
		if(bin > 0 && bin<100){
		walker[igofr+bin]+=2;
		}
     if(dr < rcut){
		    vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
				//Contribution of energy and virial
		    v += vij;
				w += wij;
     }
		}
  }  
	walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;        

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	
 	/*Uncomment to print instant values
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
	*/
	
	walker[ik]=t;
	walker[it]=(2.0 / 3.0) * t;
	walker[ie]=(t+4*v);
	
	/*Uncomment to print instant values
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
	
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
	*/
  return;
}

void Reset(int iblk) //Reset block averages
{ 
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}




void Averages(int iblk) //Print results for current block
{
    
   double err, r, gdir;
   ofstream Gofr, Gave, Epot, Pres, Ekin, Temp, Etot;
	 ifstream GofrIn;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    
    Epot.open("output.epot.0",ios::app);
    Pres.open("output.pres.0",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);
		GofrIn.open("output.gofr.0", ios::app);

  	Ekin.open("output.ekin.0",ios::app);
  	Temp.open("output.temp.0",ios::app);
 		Etot.open("output.etot.0",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

		stima_kin = blk_av[ik]/blk_norm/(double)npart; //Ekin
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

		stima_temp = blk_av[it]/blk_norm/(double)npart; //Temp
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

		stima_etot = blk_av[ie]/blk_norm/(double)npart; //Etot
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

		
//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Pressure
    Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;
//Kinetic energy per particle
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Total energy per particle
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;

//g(r)
		for (int k=igofr; k<igofr+nbins; ++k) 
		{
			//printing average value of g(r) in each block with all the g(r) mean values in each blk in a line
			r = (k-5)*bin_size;
			//normalizing each bin
			double DeltaV=0.0;
			DeltaV = 4/3*M_PI*(pow((r+bin_size),3)-pow(r,3));
			gdir = walker[k]/(DeltaV*npart*rho);
			Gofr << setw(wd) << gdir << endl;
			glob_av[k]+=gdir;
			glob_av2[k]+=gdir*gdir;
		}
		
		if (iblk==nblk){
			for(int k=igofr; k<igofr+nbins; ++k){
				err=Error(glob_av[k],glob_av2[k],iblk);
				r=(k-5)*bin_size;
				Gave << setw(wd) << r << setw(wd) << glob_av[k]/(double)nblk << setw(wd) << err << endl;
			}
		}
		
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
		Gave.close();
		GofrIn.close();
		Ekin.close();
  	Temp.close();
 		Etot.close();
}



/*-------Print configurations---------*/


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfOld(void){ //Write last but one configuration
  ofstream WriteConf;

  cout << "Print last but one configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
