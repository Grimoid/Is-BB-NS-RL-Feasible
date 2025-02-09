#ifndef GLOBALS_H
#define GLOBALS_H

#include <random>

extern std::default_random_engine generator;
extern std::uniform_real_distribution<double> distribution;

void setRandomSeed(unsigned int seed);


#endif
