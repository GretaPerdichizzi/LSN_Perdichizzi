#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include "Population.h"

using namespace std;

Population :: Population(){}

Population :: ~Population(){}

//Set random
void Population :: SetRnd(Random &random){
	rnd = &random;
	rnd->SetRandom(seed,prime[0],prime[1]);
}

//****************************************************************//

//Initialize population
void Population :: GetPop(){
	Chromosome init;
	pop.push_back(init);
	pop[0].World_Init();
	for (int i=0; i<n_chromo-1; i++){
		init.chromo=pop[0].chromo;
		random_shuffle(init.chromo.begin()+1, init.chromo.end());
		pop.push_back(init);
		//cout<<pop[0].chromo[0].x<<endl;
	}
}

//****************************************************************//

//(Sorting population by L from low to high
void Population :: Sorting(){
	bool swapped;
	for(unsigned int i=0; i<pop[0].chromo.size();i++){
		swapped = false;
		for(unsigned int j=0; j<n_chromo-i-1; j++){
			if(pop[j].L1()>pop[j+1].L1()){
				swap(pop[j].chromo, pop[j+1].chromo);
				swapped=true;
			}
		}
		if(swapped == false)
			break;
	}
}

//****************************************************************//

//Rigged-roulette searching algorithm
int Population :: Search (Random &random){
	rnd = &random; 
	//Selecting a chromosome with rigged-roulette
	double r = rnd->Rannyu();	
	int j = (int)(n_chromo*pow(r,p));
	return j;
}

//****************************************************************//

//Check function
bool Population :: CheckFunct(){
	for(unsigned int i=0; i<n_chromo; ++i){
		//check if the 1st remains unchanged
		if(pop[i].chromo[0].x!=pop[0].City1.x && pop[i].chromo[0].y!=pop[0].City1.y){
			cout <<"The first is changed!"<<endl;
			cout <<"Original 1st city x = " << pop[0].City1.x << " New 1st city x =" <<	pop[i].chromo[0].x <<endl;
			cout <<"Original 1st city x = " << pop[0].City1.y << " New 1st city x =" <<	pop[i].chromo[0].y <<endl;		 
			return false;
		}
		for(unsigned int j=0; j<n_cities-1; ++j){
			for(unsigned int s=j+1; s<n_cities; ++s){ 
				//check if there are copies
				if(pop[i].chromo[j].x==pop[i].chromo[s].x && pop[i].chromo[j].y==pop[i].chromo[s].y){
					cout<<"There are copies here :" << "pop["<<i<<"].chromo["<<j<<"].x ="<<pop[i].chromo[j].x<<" and pop["<<i<<"].chromo["<<s<<"].x =" <<pop[i].chromo[s].x<<endl;
				return false;
				}
			}
		}
	}
	return true;
}

//****************************************************************//

//Evaluate the mean of the cost function over the best half of population and print the mean and the best

void Population :: PrintFitness(int gen, int size, int rank){
	for (int j=0; j<size; j++){
		if(rank==j){
		double sum=0;
		for(int i=0; int (i<n_chromo/2); i++ ){
			sum+=pop[i].L1();
		}
		sum/=(n_chromo/2);
		ofstream WriteL;
		WriteL.open("node"+j+"fitness_mean.out", ofstream::app);
		WriteL << gen << setw(12) << sum << setw(12) << pop[0].L1() << endl;
		WriteL.close();
		}
	}
}

//****************************************************************//

//Mutate the gene
void Population :: Mutate(Random &random, int j){
	rnd = &random; 
	double p[4];

	for(int i=0; i<4; i++) p[i] = rnd->Rannyu();
	
	if(p[0]<p_mutations) pop[j].PairPerm(random);
	if(p[1]<p_mutations) pop[j].Shift(random);
	if(p[2]<p_mutations) pop[j].GroupPerm(random);
	if(p[3]<p_mutations) pop[j].Inversion(random);
	
	else ;
}

//****************************************************************//

//Crossover function
void Population :: Crossover (Random &random, int parent1, int parent2){
	rnd = &random;
	Chromosome son1, son2;
	vector <City> tail1, tail2;
	
	if(rnd->Rannyu()<=p_crossover){
		int n = (int) rnd->Rannyu(1, n_cities); //select where to do the cut
		for (int i=0; i<n_cities; i++){
			for(int j=n ; j<n_cities; j++){
				if(pop[parent2].chromo[i].x==pop[parent1].chromo[j].x && pop[parent2].chromo[i].y==pop[parent1].chromo[j].y) {
					tail1.push_back(pop[parent2].chromo[i]);
				}
				if(pop[parent1].chromo[i].x==pop[parent2].chromo[j].x && pop[parent1].chromo[i].y==pop[parent2].chromo[j].y) {
					tail2.push_back(pop[parent1].chromo[i]);
				}
			}
		}
		for (int i=0; i<n; i++){
			son1.chromo.push_back(pop[parent1].chromo[i]);
			son2.chromo.push_back(pop[parent2].chromo[i]);
		}
		for (int i=0; i<n_cities-n; i++){
			son1.chromo.push_back(tail1[i]);
			son2.chromo.push_back(tail2[i]);
		}
	}
	else{	
		son1.chromo=pop[parent1].chromo;
		son2.chromo=pop[parent2].chromo;
	}
	
	//substitute sons to the end of the population
	pop[n_chromo-2].chromo=son1.chromo;
	pop[n_chromo-1].chromo=son2.chromo;
}

//****************************************************************//

