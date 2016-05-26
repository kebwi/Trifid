#ifndef __MAIN__
#define __MAIN__

#include "FloatType.h"
#include <vector>
#include <sys/time.h>

extern std::vector<Float_t> gManualShrinkageThresholds[3];

void die(const char *msg);
void warning(const char *msg);
void myAssert(bool expression);
void assertFloatEqual(const Float_t &f1, const Float_t &f2);
int roundToNearestInt(Float_t valueF, int minV, int maxV);

#define CURRENT_TIME_PER_SEC 1000
#define CURRENT_TIME_PER_SEC_F (Float_t)CURRENT_TIME_PER_SEC
unsigned long getCurrentTime_ms();

#endif
