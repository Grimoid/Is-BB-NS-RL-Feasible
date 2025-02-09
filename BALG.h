#ifndef BALG_H
#define BALG_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>
#include <tuple>
#include <cassert>
#include "globals.h"




class GenEnvironment {
public:
    virtual std::pair<double, std::nullptr_t> step(int action) = 0;
    int time;
};


class GeoBanditEnvironment : public GenEnvironment {
public:
    int time;
    int T;
    int K;
    double delta;
    double xi;
    std::vector<std::vector<double>> Prob;
    std::vector<std::vector<int>> Rewards;

    GeoBanditEnvironment(int T);
    std::pair<double, std::nullptr_t> step(int action) override;
};

std::tuple<std::vector<int>, std::vector<double>, int, std::vector<int>> UCB(const std::vector<std::vector<int>>& Table);
std::vector<double> ComputeCumRegret(const std::vector<std::vector<double>>& Instance, const std::vector<int>& Choices);

class BALG {
public:
    int T;
    int A;
    double rho;
    int no_of_restarts;
    std::vector<int> restats;
    std::vector<int> TotalNumber;
    std::vector<double> TotalSum;
    std::vector<int> ChosenArms;
    std::vector<double> ReceivedRewards;
    GeoBanditEnvironment* environment;

    BALG(int A, int T, double restart_coeff, GeoBanditEnvironment* environment);
    std::vector<int> run();
};

class RANDALG {
public:
    int T;
    int A;
    double restart_proba;
    std::vector<int> restats;
    std::vector<int> TotalNumber;
    std::vector<double> TotalSum;
    std::vector<int> ChosenArms;
    std::vector<double> ReceivedRewards;
    GeoBanditEnvironment* environment;

    RANDALG(int A, int T, double restart_proba, GeoBanditEnvironment* environment);
    std::vector<int> run();
};

#endif
