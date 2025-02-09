#include "BALG.h"
#include "problems.h"
#include "globals.h"
#include "algorithms.h"




GeoBanditEnvironment::GeoBanditEnvironment(int T) : time(0), T(T), K(5), delta(1.0 / T), xi(0.9) {
    double Deltamin = 0.1;
    double Deltamax = 0.4;
    int num_breaks = 0;

    while (num_breaks == 0) {
        std::tie(Prob, std::ignore) = ProblemGeo(T, K, xi, Deltamin, Deltamax);
        num_breaks = Prob.size();
    }
    Rewards = CreateRewards(Prob);
}

std::pair<double, std::nullptr_t> GeoBanditEnvironment::step(int action) {
    time++;
    return {Rewards[action][time - 1], nullptr};
}

std::tuple<std::vector<int>, std::vector<double>, int, std::vector<int>> UCB(const std::vector<std::vector<int>>& Table) {
    int K = Table.size();
    int T = Table[0].size();
    int episode = 1;
    std::vector<int> ChangePoints;
    std::vector<int> counts(K, 1);
    std::vector<double> values(K, 0);
    std::vector<int> ChosenArms(T);
    std::vector<double> ReceivedRewards(T);

    for (int t = 0; t < T; ++t) {
        std::vector<double> avg_rewards(K);
        std::vector<double> exploration_terms(K);
        for (int i = 0; i < K; ++i) {
            avg_rewards[i] = values[i] / counts[i];
            exploration_terms[i] = std::sqrt(2 * std::log(std::accumulate(counts.begin(), counts.end(), 0)) / counts[i]);
        }
        std::vector<double> ucb_values(K);
        for (int i = 0; i < K; ++i) {
            ucb_values[i] = avg_rewards[i] + exploration_terms[i];
        }
        for (int i = 0; i < K; ++i) {
            if (counts[i] == 0) {
                ucb_values[i] = std::numeric_limits<double>::infinity();
            }
        }
        int arm = std::distance(ucb_values.begin(), std::max_element(ucb_values.begin(), ucb_values.end()));
        int reward = Table[arm][t];
        ReceivedRewards[t] = reward;
        ChosenArms[t] = arm;
        counts[arm]++;
        values[arm] += reward;
    }
    return std::make_tuple(ChosenArms, ReceivedRewards, episode, ChangePoints);
}



BALG::BALG(int A, int T, double restart_coeff, GeoBanditEnvironment* environment)
    : A(A), T(T), environment(environment) {
    no_of_restarts = std::pow(T, restart_coeff);
    rho = (float)no_of_restarts/(T*std::sqrt(std::log(T)));

    int nbreaks = 0;
    while (nbreaks < 1) {
        restats.clear();
        for (int i = 0; i <= T; ++i) {
            if (distribution(generator) <= rho) restats.push_back(i);
        }
        nbreaks = restats.size();
    }
    TotalNumber.resize(A, 0);
    TotalSum.resize(A, 0.0);
    ChosenArms.resize(T);
    ReceivedRewards.resize(T);


}

std::vector<int> BALG::run() {
    int restart_ind = 0;
    for (int t = 1; t <= T; ++t) {
        if (restart_ind < restats.size() && t - 1 == restats[restart_ind]) {
            std::fill(TotalNumber.begin(), TotalNumber.end(), 0);
            std::fill(TotalSum.begin(), TotalSum.end(), 0.0);
            restart_ind++;
        }

        std::vector<double> indices(A, 0);
        for (int i = 0; i < A; ++i) {
            if (TotalNumber[i] > 0) {
                indices[i] = klUp(TotalSum[i] / TotalNumber[i], std::log(t + 1) / TotalNumber[i]);
            }
            else {
            indices[i] = std::numeric_limits<double>::infinity();
        }
        }

        int I = randmax(indices);
        auto [rew, _] = environment->step(I);

        ChosenArms[t - 1] = I;
        ReceivedRewards[t - 1] = rew;
        TotalNumber[I]++;
        TotalSum[I] += rew;
    }
    return ChosenArms;
}

RANDALG::RANDALG(int A, int T, double restart_proba, GeoBanditEnvironment* environment)
    : A(A), T(T), restart_proba(restart_proba), environment(environment) {
    std::binomial_distribution<> d(T, restart_proba);
   
    int nbreaks = 0;
    while (nbreaks < 1) {
        restats.clear();
        for (int i = 0; i <= T; ++i) {
            if (distribution(generator) <= restart_proba) restats.push_back(i);
        }
        nbreaks = restats.size();
    }


    TotalNumber.resize(A, 0);
    TotalSum.resize(A, 0.0);
    ChosenArms.resize(T);
    ReceivedRewards.resize(T);
}

std::vector<int> RANDALG::run() {
    int restart_ind = 0;

    for (int t = 1; t <= T; ++t) {
        if (restart_ind < restats.size() && t - 1 == restats[restart_ind]) {
            std::fill(TotalNumber.begin(), TotalNumber.end(), 0);
            std::fill(TotalSum.begin(), TotalSum.end(), 0.0);
            restart_ind++;
        }
        std::vector<double> indices(A, 0);
        for (int i = 0; i < A; ++i) {
            if (TotalNumber[i] > 0) {
                indices[i] = klUp(TotalSum[i] / TotalNumber[i], std::log(t + 1) / TotalNumber[i]);
            }
            else {
            indices[i] = std::numeric_limits<double>::infinity();
        }
        }

        int I = randmax(indices);
        auto [rew, _] = environment->step(I);

        ChosenArms[t - 1] = I;
        ReceivedRewards[t - 1] = rew;
        TotalNumber[I]++;
        TotalSum[I] += rew;
    }
    return ChosenArms;
}

