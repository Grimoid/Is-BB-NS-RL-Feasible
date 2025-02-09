#ifndef UCB_H
#define UCB_H

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <random>
#include <functional>

class UCB {
public:
    UCB(int T, double delta, int s, int e);
    std::tuple<int, double, std::nullptr_t> select_arm();
    void update_state(int x, double y);
    void re_init();
    std::vector<int> return_params();
    friend std::ostream& operator<<(std::ostream& os, const UCB& ucb);

    int s; 
    int e; 
    int len; 

private:
    int T;
    double delta;
    int num_arms;
    int t;
    int ctr;
    int lazy_update_fr;
    double c;
    std::vector<int> counts;
    std::vector<double> values;
};



#endif
