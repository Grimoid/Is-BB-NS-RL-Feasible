#include "UCB.h"
#include "algorithms.h"
using namespace std;

UCB::UCB(int T, double delta, int s, int e) 
    : T(T), delta(delta), s(s), e(e), num_arms(5), t(0), ctr(0), lazy_update_fr(5), c(2.0) {
    len = e - s + 1;
    counts.resize(num_arms, 0);
    values.resize(num_arms, 0.0);
}

tuple<int, double, nullptr_t> UCB::select_arm() {
    vector<double> avg_rewards(num_arms, 0.0);
    for (int i = 0; i < num_arms; ++i) {
        if (counts[i] != 0) {
            avg_rewards[i] = values[i] / counts[i];
        }
    }

    vector<double> exploration_terms(num_arms, 0.0);
    for (int i = 0; i < num_arms; ++i) {
        if (counts[i] != 0) {
            exploration_terms[i] = c * sqrt(log(T / delta) / counts[i]);
        }
    }

    vector<double> ucb_values(num_arms, 0.0);
    for (int i = 0; i < num_arms; ++i) {
        ucb_values[i] = avg_rewards[i] + exploration_terms[i];
    }

    for (int i = 0; i < num_arms; ++i) {
        if (counts[i] == 0) {
            ucb_values[i] = c * sqrt(log(T / delta));
        }
    }

    auto max_it = max_element(ucb_values.begin(), ucb_values.end());
    int chosen_arm = distance(ucb_values.begin(), max_it);
    return make_tuple(chosen_arm, *max_it, nullptr);
}

void UCB::update_state(int x, double y) {
    counts[x] += 1;
    values[x] += y;
    t += 1;
}

void UCB::re_init() {
    t = 0;
    ctr = 0;
    fill(counts.begin(), counts.end(), 0);
    fill(values.begin(), values.end(), 0.0);
}

vector<int> UCB::return_params() {
    return {s + 1, e + 1, static_cast<int>(log2(len))};
}

ostream& operator<<(ostream& os, const UCB& ucb) {
    os << "[start: " << ucb.s << ", end: " << ucb.e << ", m: " << static_cast<int>(log2(ucb.len)) << "]";
    return os;
}






