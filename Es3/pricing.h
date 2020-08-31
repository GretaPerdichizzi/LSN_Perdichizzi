#include "random.h"

#ifndef __Pricing__
#define __Pricing__

class Pricing : private Random {

private:
  int S0;
	double r, sigma, T, tstep;

protected:

public:
	int K;
  // constructors
  Pricing();
  // destructor
  ~Pricing();
  // methods
	void SetParams(double Tfinal);
 	double AssetPrice(bool discr, Random & rnd);
	double CallOption(double S);	
  double PutOption(double S);
};

#endif // __Pricing__

