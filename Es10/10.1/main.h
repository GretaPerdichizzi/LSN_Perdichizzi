#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"
using namespace std;

const int n_chromo=100;

//Random numbers
int seed[4];
Random rnd;

//Usefull functions
void Input(void);
void PrintConfig(vector<City> c);

