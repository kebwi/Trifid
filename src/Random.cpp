#include "Random.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

static const Float_t Pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862;
static const Float_t PiTimes2 = Pi * 2.0;

static unsigned long sseed;

//Seed with the time (a random seed) and return the time used to seed.
unsigned long SeedRandom()
{
	sseed = time(NULL);
	srand(sseed);
	return sseed;
}

//Seed with an explicit seed.
void SeedRandom(int val)
{
	srand(val);
}

//Retrieve last seed used.
unsigned long GetLastRandomSeed()
{
	return sseed;
}

//Return an int in the range (0 - range].  range must be (0 - 32768).
int RandZeroInt(int range)
{
	return rand() % range;
}

int RandNegInt(int range)
{
	//Return an int in the range (-range/2 - range/2)
	//range must be (0 - 32768)
	//if (range % 2 == 0)
	//	return 0;	//if the range isn't odd, this algorithm doesn't really work.
	//return rand() % range - range / 2;
	
	//Return an int in the range [-range - range]
	//range must be less than half the max size of an int, about 2 billion
	return rand() % ((range << 1) + 1) - range;
}

//Return a Float_t in the range [0 - range).
Float_t RandZeroFloat(Float_t range)
{
	//return ((Float_t)rand() / ((Float_t)(RAND_MAX + 1) / range));
	
	return (((Float_t)rand() / (Float_t)RAND_MAX) * range);
}

//Return a Float_t in the range [-range - range].
Float_t RandNegFloat(Float_t range)
{
	return (((Float_t)rand() / (Float_t)RAND_MAX * (range * 2.0)) - range);
}

//Return a Float_t from a normal (Gaussian) distribution
Float_t RandStdNorm()
{
	//Standard normal distribution
	
	//Box-Muller
	//http://en.wikipedia.org/wiki/Normal_distribution
	
	//Float_t a = RandZeroFloat(1.0);
	//Float_t b = RandZeroFloat(1.0);
	//Float_t c = sqrt(-2 * log(a)) * cos(2 * Pi * b);
	//return c;
	
	//return sqrt(-2.0 * log(RandZeroFloat(1.0))) * cos(PiTimes2 * RandZeroFloat(1.0));
	
	//================================================================================
	
	//Polar form of Box-Muller
	//http://www.taygeta.com/random/gaussian.html
	
	Float_t x1, x2, w, y1;
	static Float_t y2;
	static bool use_last = false;

	if (use_last)	//Use value from previous call
	{
		y1 = y2;
		use_last = false;
	}
	else
	{
		do
		{
			x1 = 2.0 * RandZeroFloat(1.0) - 1.0;
			x2 = 2.0 * RandZeroFloat(1.0) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
		
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = true;
	}

	return y1;//(y1 * stddev + mean);
}

//Return a Float_t from a general normal (Gaussian) distribution
Float_t RandGenNorm(Float_t mean, Float_t stddev)
{
	return RandStdNorm() * stddev + mean;
}
