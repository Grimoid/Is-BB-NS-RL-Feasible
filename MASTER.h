#ifndef MASTER_H
#define MASTER_H

#include "UCB.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <random>
#include <functional>
#include "globals.h"

class MASTER {
public:
    MASTER(int K, double delta, int T, int n_init = 0);
    int select_arm();
    void update_state(int x, double y);
    void re_init();
    std::vector<int> run(const std::vector<std::vector<int>>& Environment);
    std::vector<UCB> policy;

protected:
    int K;
    double delta;
    int T;
    std::string model;
    int n_init;
    double c;
    int t;
    int n;
    int tn;
    double n_hat;
    int alg_index;
    std::vector<double> g_tilde;
    std::vector<int> alge;
    std::vector<double> reward;
    std::vector<std::vector<int>> active_list;
    std::function<double(double)> rho;
    std::function<double(double)> rho_hat;

    virtual void procedure1();
    int test1();
    int test2();
    void restart();
};



#endif
