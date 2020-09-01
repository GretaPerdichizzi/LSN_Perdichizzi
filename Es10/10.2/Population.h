#ifndef __POPULATION__
#define __POPULATION__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "Chromosome.h"
#include "random.h"

using namespace std;

class Population{

	private:
		static const int n_chromo=100;	
		static const int n_cities=32;
		vector <double> fitness;	
		Random *rnd;
		int seed[4];
		
		double p=2;
		double p_crossover=0.8;
		double p_mutations=0.15;

	public:
		vector<Chromosome> pop;
		vector<double> best_path_x;
		vector<double> best_path_y;
		int prime[2];

		// Constructor
		Population();

		//Destructor
		~Population();

		//Copy Constructor
		Population(const Population&);

		//Methods
		void SetRnd(Random &random);
		void GetPop();
		void Sorting();
		bool CheckFunct();
		void PrintFitness(int gen, int size, int rank);
		int Search(Random &random);
		void Mutate(Random &random, int j);
		void GetBest();
		void CopyBest();
		
		//Crossover
		void Crossover(Random &random, int parent1, int parent2);
};


#endif
