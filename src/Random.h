#ifndef __RANDOM__
#define __RANDOM__

#include "FloatType.h"

//Seed with the time (a random seed) and return the time used to seed.
unsigned long SeedRandom();

//Seed with an explicit seed.
void SeedRandom(int val);

//Retrieve last seed used.
unsigned long GetLastRandomSeed();

//Return an int in the range (0 - range].  range must be (0 - 32768).
int RandZeroInt(int range);

//Return an int in the range (-range/2 - range/2).  range must be (0 - 32768).
//Return an int in the range [-range - range].
int RandNegInt(int range);

//Return a Float_t in the range (0 - range).
Float_t RandZeroFloat(Float_t range);

//Return a Float_t in the range (-range - range).
Float_t RandNegFloat(Float_t range);

//Return a Float_t from a standard normal (Gaussian) distribution
Float_t RandStdNorm();

//Return a Float_t from a general normal (Gaussian) distribution
Float_t RandGenNorm(Float_t mean, Float_t stddev);

#endif
