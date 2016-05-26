#ifndef __FLOAT_TYPE__
#define __FLOAT_TYPE__

//When processing steerable pyramids...
//float and double run three to five times faster than long double.
//float doesn't run noticably faster (or slower) than double.
//Obviously, these figure may vary depending on the machine and the compiler optimizations.

//Use this single typedef to set the float type for the entire project.

//Leave only one of the following typedefs uncommented:
//typedef float Float_t;
typedef double Float_t;
//typedef long double Float_t;

#endif
