#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <vector>
#include <random>
#include <tuple>
#include <string>
#include <functional>

int randmax(const std::vector<double>& vec, int rank = 1);
double kl(double p, double q);
double klSG(double p, double q);
double klUp(double p, double level);

int klUCBBandit(const std::vector<int>& ArmList, const std::vector<double>& RewardList, int tloc);
int UCBBandit(const std::vector<int>& ArmList, const std::vector<double>& RewardList, int tloc);

std::vector<std::vector<int>> CreateRewards(const std::vector<std::vector<double>>& Problem);
std::vector<double> ComputeCumRegret(const std::vector<std::vector<double>>& Instance, const std::vector<int>& Choices);




std::tuple<std::vector<int>, std::vector<double>, int, std::vector<int>> QCD_PLUS(
    const std::vector<std::vector<int>>& Table, 
    int(*banditAlgorithm)(const std::vector<int>&, const std::vector<double>&, int), 
    double(*beta)(int, double), 
    double(*klChange)(double, double), 
    double delta=0.05, 
    int subSample=10, 
    int subSample2=5, 
    std::string verbose="off"
);



std::tuple<std::vector<int>, std::vector<double>, int, std::vector<int>> GLRklUCBGlobal(
    const std::vector<std::vector<int>>& Table,
    double(*beta)(int, double),
    double(*klChange)(double, double),
    double alpha0=0.1,
    double delta=0.05,
    std::string alphatype="auto",
    int subSample=10,
    int subSample2=5,
    std::string verbose="off"
);



#endif
