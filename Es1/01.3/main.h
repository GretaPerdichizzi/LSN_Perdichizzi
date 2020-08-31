#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include "random.h"

Random rnd;
int seed[4];

const double L = 1.5; //needle lenght
const double d = 2.5; //distance between lines
int m = 10000; //number of throws
int n = 100; //number of blocks

void Input(void);
double Simulate(void);


