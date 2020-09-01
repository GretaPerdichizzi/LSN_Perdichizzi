#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include "random.h"
using namespace std;

//Random numbers
int seed[4];
int prime[8];
Random rnd;


void PrintConfig(vector<City> c, int size, int rank);
void Input(int size);

