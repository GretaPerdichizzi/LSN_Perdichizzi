#ifndef __CHROMOSOME__
#define __CHROMOSOME__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "City.h"
#include "random.h"

using namespace std;

class Chromosome{

	private:	
	static const int n_chromo=100;	
	static const int n_cities=32;
	vector<City> p_world;
	City coord;	
	Random *rnd;
	
	public:
		// Constructor
		Chromosome();

		//Destructor
		~Chromosome();

		//Variables
		vector<City> chromo;
		double lenghts;
		City City1;

		//Methods
		void World_Init();
		double L1();
		double Lval;

		//Mutation Methods
		void PairPerm(Random &random);
		void Shift(Random &random);
		void GroupPerm(Random &random);
		void Inversion(Random &random);

};

#endif

