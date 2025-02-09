#include "globals.h"

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);
void setRandomSeed(unsigned int seed) {
    generator.seed(seed);
}