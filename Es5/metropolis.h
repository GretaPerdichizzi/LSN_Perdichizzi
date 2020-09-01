#include "random.h"

void metropolis (Random rnd, double q, double xn, double yn, double zn, double &x1, double &y1, double &z1, int &accept, int &reject);
double q_eval_100 (double xn, double yn, double zn, double x1, double y1, double z1);
double q_eval_210 (double xn, double yn, double zn, double x1, double y1, double z1);
