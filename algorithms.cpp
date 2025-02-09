#include "algorithms.h"
#include "globals.h"
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <string>
#include <functional>
#include <iostream>

using namespace std;


template <typename T>
void print1DMatrix(const T& matrix) {
    for (const auto& element : matrix) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}


int randmax(const vector<double>& vec, int rank) {
    vector<double> sorted_vec = vec;
    sort(sorted_vec.rbegin(), sorted_vec.rend());
    double m = sorted_vec[rank - 1];
    vector<int> Ind;
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] == m) Ind.push_back(i);
    }
    uniform_int_distribution<int> dist(0, Ind.size() - 1);
    return Ind[dist(generator)];
}

double kl(double p, double q) {
    if (p == q) return 0;
    p = max(p, numeric_limits<double>::epsilon());
    p = min(p, 1 - numeric_limits<double>::epsilon());
    return p * log(p / q) + (1 - p) * log((1 - p) / (1 - q));
}

double klSG(double p, double q) {
    return 2 * pow(p - q, 2);
}

double klUp(double p, double level) {
    double lM = p;
    double uM = min(1.0, p + sqrt(level / 2));
    for (int i = 0; i < 50; ++i) {
        double qM = (uM + lM) / 2;
        if (kl(p, qM) > level) uM = qM;
        else lM = qM;
    }
    return uM;
}


vector<vector<int>> CreateRewards(const vector<vector<double>>& Problem) {
    int K = Problem.size();
    int T = Problem[0].size();
    vector<vector<int>> Table(K, vector<int>(T, 0));
    for (int i = 0; i < K; ++i) {
        for (int t = 0; t < T; ++t) {
            Table[i][t] = distribution(generator) < Problem[i][t];
        }
    }
    return Table;
}

std::vector<double> findColumnMaxima(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty() || matrix[0].empty()) {
        return {}; 
    }

    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size();
    std::vector<double> columnMaxima(numCols, std::numeric_limits<double>::lowest()); 

    for (size_t col = 0; col < numCols; ++col) {
        for (size_t row = 0; row < numRows; ++row) {
            if (matrix[row][col] > columnMaxima[col]) {
                columnMaxima[col] = matrix[row][col];
            }
        }
    }

    return columnMaxima;
}

vector<double> ComputeCumRegret(const vector<vector<double>>& Instance, const vector<int>& Choices) {
    int K = Instance.size();
    int T = Instance[0].size();
    vector<double> Regret(T, 0.0);
    double regret = 0;

    vector <double> Maxima(T,0.0);
    
    Maxima=findColumnMaxima(Instance);

    for (int t = 0; t < T; ++t) {
        double check;
        check=Maxima[t] - Instance[Choices[t]][t];
        regret += check;
        Regret[t] = regret;
    }
    return Regret;
}



int klUCBBandit(const vector<int>& ArmList, const vector<double>& RewardList, int tloc) {
    int nArms = ArmList.size();
    vector<double> indices(nArms, 0.0);
    for (int i = 0; i < nArms; ++i) {
        if (ArmList[i] > 0) {
            indices[i] = klUp(RewardList[i] / ArmList[i], log(tloc + 1) / ArmList[i]);
        }
        else return 1234567;
    }
    return randmax(indices);
}

int UCBBandit(const vector<int>& ArmList, const vector<double>& RewardList, int tloc) {
    int num_arms = ArmList.size();
    vector<double> ucb_values(num_arms, 0.0);
    int total_pulls = accumulate(ArmList.begin(), ArmList.end(), 0);

    for (int i = 0; i < num_arms; ++i) {
        if (ArmList[i] > 0) {
            double avg_reward = RewardList[i] / ArmList[i];
            double exploration_term = sqrt(2 * log(total_pulls) / ArmList[i]);
            ucb_values[i] = avg_reward + exploration_term;
        }
    }

    for (int i = 0; i < num_arms; ++i) {
        if (ArmList[i] == 0) {
            ucb_values[i] = numeric_limits<double>::infinity();
        }
    }

    return distance(ucb_values.begin(), max_element(ucb_values.begin(), ucb_values.end()));
}


tuple<vector<int>, vector<double>, int, vector<int>> QCD_PLUS(const vector<vector<int>>& Table, int(*banditAlgorithm)(const vector<int>&, const vector<double>&, int), double(*beta)(int, double), double(*klChange)(double, double), double delta, int subSample, int subSample2, string verbose) {
    int K = Table.size();
    int T = Table[0].size();
    int episode = 1;
    vector<int> ChangePoints;
    vector<int> ChosenArms(T, 0);
    vector<double> ReceivedRewards(T, 0.0);
    int t = 0;
    int tau = 0;

    while (t < T) {
        bool newEpisode = false;
        
        vector<vector<double>> SUMS_CD(K);
        int tloc = 0;

        vector<int> TotalNumber_CD(K, 0);
        vector<double> TotalSum_CD(K, 0.0);

        vector<int> TotalNumber_B(K, 0);
        vector<double> TotalSum_B(K, 0.0);
        while (!newEpisode && t < T) {
            int I = banditAlgorithm(TotalNumber_B, TotalSum_B, tloc);
            if(I==1234567){
            if (accumulate(TotalNumber_B.begin(), TotalNumber_B.end(), 0) == 0) {
                for (int arm = 0; arm < K; ++arm) {
                    int rew = Table[arm][t];
                    TotalNumber_B[arm] += 1;
                    TotalSum_B[arm] += rew;
                }
            }
            I = banditAlgorithm(TotalNumber_B, TotalSum_B, tloc);

            }
            

            double rew = Table[I][t];
            ChosenArms[t] = I;
            ReceivedRewards[t] = rew;
            t++;
            tloc++;

            TotalNumber_B[I] += 1;
            TotalSum_B[I] += rew;
            
            TotalNumber_CD[I] += 1;
            TotalSum_CD[I] += rew;
            SUMS_CD[I].push_back(TotalSum_CD[I]);

            int check = 0;
            if (tloc % subSample == 1) {
                int s = 1;
                int nb = TotalNumber_CD[I];
                vector<double> sums = SUMS_CD[I];
                while (s < nb && check == 0) {
                    if (s % subSample2 == 1) {
                        int draw1 = s;
                        int draw2 = nb - s;
                        double mu1 = sums[s - 1] / draw1;
                        double mu2 = (sums[nb - 1] - sums[s - 1]) / draw2;
                        double mu = sums[nb - 1] / nb;
                        if (draw1 * klChange(mu1, mu) + draw2 * klChange(mu2, mu) > beta(nb, delta)) {
                            newEpisode = true;
                            check++;
                        }
                    }
                    s++;
                }
                if (check > 0) {
                    episode++;
                    ChangePoints.push_back(t);
                    tau = t;

                    if (verbose == "on") {
                        cout << "detected a change on " << I + 1 << " at t=" << t << endl;
                    }
                }
            }
        }
    }

    return make_tuple(ChosenArms, ReceivedRewards, episode, ChangePoints);
}




tuple<vector<int>, vector<double>, int, vector<int>> GLRklUCBGlobal(
    const vector<vector<int>>& Table,
    double(*beta)(int, double),
    double(*klChange)(double, double),
    double alpha0,
    double delta,
    string alphatype,
    int subSample,
    int subSample2,
    string verbose
) {
    int K = Table.size();
    int T = Table[0].size();
    int episode = 1;
    vector<int> ChangePoints;
    vector<int> ChosenArms(T, 0);
    vector<double> ReceivedRewards(T, 0.0);


    int(*banditAlgorithm)(const vector<int>&, const vector<double>&, int)=klUCBBandit;
    int t = 0;
    int tau = 0;

    while (t < T) {
        bool newEpisode = false;
        double alpha = alpha0 * sqrt(episode);
        if (alphatype == "constant") alpha = alpha0;
        int ExploRange;
        if(alpha!=0)ExploRange = ceil(K / alpha);
        else ExploRange =1;
        vector<vector<double>> SUMS_CD(K);
        int tloc = 0;

        vector<int> TotalNumber(K, 0);
        vector<double> TotalSum(K, 0.0);
        while (!newEpisode && t < T) {
            int I = (t - tau) % ExploRange;
            if (I >= K || alpha==0) {
                I = banditAlgorithm(TotalNumber, TotalSum, tloc);
                if(I==1234567){
                if (accumulate(TotalNumber.begin(), TotalNumber.end(), 0) == 0) {
                    for (int arm = 0; arm < K; ++arm) {
                        int rew = Table[arm][t];
                        TotalNumber[arm] += 1;
                        TotalSum[arm] += rew;
                    }
                }
                I = banditAlgorithm(TotalNumber, TotalSum, tloc);

                }
            }

            double rew = Table[I][t];
            ChosenArms[t] = I;
            ReceivedRewards[t] = rew;
            t++;
            tloc++;

            
            TotalNumber[I] += 1;
            TotalSum[I] += rew;
                
            SUMS_CD[I].push_back(TotalSum[I]);

            int check = 0;
            if (tloc % subSample == 1) {
                int s = 1;
                int nb = TotalNumber[I];
                vector<double> sums = SUMS_CD[I];
                while (s < nb && check == 0) {
                    if (s % subSample2 == 1) {
                        int draw1 = s;
                        int draw2 = nb - s;
                        double mu1 = sums[s - 1] / draw1;
                        double mu2 = (sums[nb - 1] - sums[s - 1]) / draw2;
                        double mu = sums[nb - 1] / nb;
                        if (draw1 * klChange(mu1, mu) + draw2 * klChange(mu2, mu) > beta(nb, delta)) {
                            newEpisode = true;
                            check++;
                        }
                    }
                    s++;
                }
                if (check > 0) {
                    episode++;
                    ChangePoints.push_back(t);
                    tau = t;
                    if (verbose == "on") {
                        cout << "detected a change on " << I + 1 << " at t=" << t << endl;
                    }
                }
            }
        }
    }

    return make_tuple(ChosenArms, ReceivedRewards, episode, ChangePoints);
}


