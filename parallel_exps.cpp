#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <numeric>
#include <tuple>
#include <chrono>
#include <string>
#include <fstream>
#include <functional>
#include "globals.h"
#include "problems.h"
#include "algorithms.h"
#include "UCB.h"
#include "MASTER.h"
#include "BALG.h"
#include <omp.h>
#include <sys/stat.h> 
#include <cstdlib>    
#include <unistd.h>   



using namespace std;


auto beta_func = [](int n, double delta) { return  log(((4 * pow(n, 1.5)) / delta)) ; };


struct Policy {
    function<tuple<vector<int>, vector<double>, int, vector<int>>(const vector<vector<int>> &)> algorithm;
    string name;
};

tuple<vector<Policy>, vector<string>> init_algos(int T, int K, const vector<int> &BreakPoints, double Deltamin, double tweak, double xi, double explor_freq, bool verbose = false) {
    double alpha0 = tweak * sqrt(K * log(T) / T);
    double delta_GLR = 1.0 / std::pow(T,0.5);
    double delta_GLR_BK=1.0/sqrt(T);



    
    auto GLRKLUCB = [=](const vector<vector<int>> &Table) {
        return GLRklUCBGlobal(Table, beta_func, kl, alpha0, delta_GLR_BK);
    };





    auto QCD_PLUS_KLUCB = [=](const vector<vector<int>> &Table) {
        return QCD_PLUS(Table, klUCBBandit, beta_func, kl, delta_GLR);
    };


    auto QCD_PLUS_UCB = [=](const vector<vector<int>> &Table) {
        return QCD_PLUS(Table, UCBBandit, beta_func, kl, delta_GLR);
    };



    auto RR_OPT = [=](const vector<vector<int>> &Table) {
        GeoBanditEnvironment temp_environment(T);
        temp_environment.Rewards = Table;
        double restart_coeff = 1 - xi / 2;
        BALG balg(K, T, restart_coeff, &temp_environment);
        auto ChosenArms = balg.run();
        return make_tuple(ChosenArms, vector<double>(ChosenArms.size(), 0.0), 0, vector<int>());
    };
    auto RANDALG_p_0_05 = [=](const vector<vector<int>> &Table) {
        GeoBanditEnvironment temp_environment(T);
        temp_environment.Rewards = Table;
        double restart_proba = 0.05;
        RANDALG randalg(K, T, restart_proba, &temp_environment);
        auto ChosenArms = randalg.run();
        return make_tuple(ChosenArms, vector<double>(ChosenArms.size(), 0.0), 0, vector<int>());
    };

    auto MASTER_ALGO = [=](const vector<vector<int>> &Table) {
        int n_init = static_cast<int>(floor(log2(T))) - 4;
        MASTER master(K, 1.0 / T, T, n_init);
        auto ChosenArms = master.run(Table);
        return make_tuple(ChosenArms, vector<double>(ChosenArms.size(), 0.0), 0, vector<int>());
    };

    
    

     vector<Policy> policies = {
         {MASTER_ALGO, "MASTER_ALGO"},
        {GLRKLUCB,"GLRKLUCB"},
         {RANDALG_p_0_05, "RANDALG_p_0_05"},
         {RR_OPT, "RR_OPT"},
         {QCD_PLUS_KLUCB, "QCD_PLUS_KLUCB"},
         {QCD_PLUS_UCB, "QCD_PLUS_UCB"}
    };

    vector<string> names = {
         "MASTER_ALGO",
         "GLRKLUCB",
         "RANDALG_p_0_05",
         "RR_OPT",
         "QCD_PLUS_KLUCB",
         "QCD_PLUS_UCB"
    };

    return make_tuple(policies, names);
}


template <typename T>
void print1DMatrix(const T& matrix) {
    for (const auto& element : matrix) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}

void save_results(const string &name, const vector<vector<double>> &Regret, const vector<int> &Restarts, double avg_time, double stddev_time) {
    ofstream file(name);
    if (file.is_open()) {
        file << "Regret:\n";
        for (const auto& row : Regret) {
            for (const auto& val : row) {
                file << val << " ";
            }
            file << "\n";
        }
        file << "Restarts:\n";
        for (const auto& val : Restarts) {
            file << val << " ";
        }
        file << "\n";
        
        file << "MTPE:\n" << avg_time << "\n";
        file << "STDEVPE:\n" << stddev_time << "\n";
        
        file.close();
    } else {
        cerr << "Unable to open file";
    }
}

void BayesianExpes(int N, int K, int T, double xi, double Deltamin, double Deltamax, double tweak, double explor_freq, double Deltachange = 0.1, bool saveres = true, const string &PB = "Normal",const string &directory = "./cpp_results/") {
    string file_name = directory + PB + "_Geo_beta_tweak_" + to_string(tweak);
    int num_breaks = 0;
    vector<tuple<vector<vector<double>>, vector<int>>> Problems;
    vector<vector<int>> BreakPoint_list;

    

    while (num_breaks == 0) {
        if (PB == "Normal") {
            Problems.push_back(ProblemGeo(T, K, xi, Deltamin, Deltamax));
        } else if (PB == "Worst") {
            Problems.push_back(WorstProblemGeo(T, K, xi, Deltamin, Deltamax));
        }else if (PB == "Normal_Unif") {
            Problems.push_back(ProblemGeoUnif(T, K, xi, Deltamin, Deltamax));
        }else if (PB == "Worst_Unif") {
            Problems.push_back(WorstProblemGeoUnif(T, K, xi));
        }
        num_breaks = get<1>(Problems.back()).size();
    }
    BreakPoint_list.push_back(get<1>(Problems.back()));

    auto [policies, names] = init_algos(T, K, get<1>(Problems[0]), Deltamin, tweak, xi, explor_freq);
    int lP = policies.size();
    vector<int> tsave;
    for (int i = 1; i < T; i += T / 500) {
        tsave.push_back(i);
    }
    int ts = tsave.size();

    // Generate the Problems in parallel
    #pragma omp parallel for
    for (int n = 1; n < N; ++n) {
        int num_breaks = 0;
        while (num_breaks == 0) {
            if (PB == "Normal") {
                #pragma omp critical
                Problems.push_back(ProblemGeo(T, K, xi, Deltamin, Deltamax));
            } else if (PB == "Worst") {
                #pragma omp critical
                Problems.push_back(WorstProblemGeo(T, K, xi, Deltamin, Deltamax));
            }
            #pragma omp critical
            num_breaks = get<1>(Problems.back()).size();
        }
        #pragma omp critical
        BreakPoint_list.push_back(get<1>(Problems.back()));
    }

    vector<vector<vector<int>>> Tables;
    for (const auto &Pb : Problems) {
        Tables.push_back(CreateRewards(get<0>(Pb)));
    }

    // Parallelize the policy execution
    #pragma omp parallel for
    for (int imeth = 0; imeth < lP; ++imeth) {
        auto start = chrono::high_resolution_clock::now();
        auto policy = policies[imeth];
        string policy_name = names[imeth];
        string name = file_name + "_" + policy_name + "_T_" + to_string(T) + "_N_" + to_string(N) + "_K_" + to_string(K) + "_xi_" + to_string(xi) + ".txt";
        vector<vector<double>> Regret(N, vector<double>(ts, 0.0));
        vector<int> Restarts(N, 0);

        std::cout << "Running Experiments for: " << policy_name << endl;
        vector<double> times(N, 0.0);  // to store time taken for each experiment
        // Parallelize the experiments across different problems
        #pragma omp parallel for
        for (int n = 0; n < N; ++n) {
            vector<int> ChosenArms;
            vector<int> cpd;
            vector<vector<int>> curr_table;

            tuple<vector<vector<double>>, vector<int>> curr_problem;
            if(PB=="Normal_Unif"||PB=="Worst_Unif"){
                curr_table=Tables[0];
                curr_problem=Problems[0];
                }

            else {
                curr_table=Tables[n];
                curr_problem=Problems[n];
                }

            auto start_exp = chrono::high_resolution_clock::now();
            if (policy_name == "RR_OPT") {
                GeoBanditEnvironment temp_environment(T);
                temp_environment.Rewards = curr_table;
                double restart_coeff = 1 - xi / 2;
                BALG balg(K, T, restart_coeff, &temp_environment);
                ChosenArms = balg.run();
            } else if (policy_name.find("RANDALG") != string::npos) {
                GeoBanditEnvironment temp_environment(T);
                temp_environment.Rewards = curr_table;
                double restart_proba = stod(policy_name.substr(policy_name.find_last_of('_') + 1));
                if(restart_proba==5.0)restart_proba=0.05;
                RANDALG randalg(K, T, restart_proba, &temp_environment);
                ChosenArms = randalg.run();
            } else if (policy_name == "MASTER_ALGO") {
                int n_init = static_cast<int>(floor(log2(T)));
                MASTER master(K, 1.0 / T, T, n_init);
                ChosenArms = master.run(curr_table);
            } else {
                if (policy.algorithm) {
                    tie(ChosenArms, ignore, ignore, cpd) = policy.algorithm(curr_table);
                }
            }

            auto reg = ComputeCumRegret(get<0>(curr_problem), ChosenArms);
            for (int t = 0; t < ts; ++t) {
                Regret[n][t] = reg[tsave[t]];
            }
            auto end_exp = chrono::high_resolution_clock::now();
            chrono::duration<double> elapsed = end_exp - start_exp;
            times[n] = elapsed.count(); 
            Restarts[n] = cpd.size();
            if (n % 500 == 0) {
                std::cout << "Results for " << policy_name << " iteration:" << n << endl;
                std::cout << "Regret is " << accumulate(Regret[n].begin(), Regret[n].end(), 0.0) / Regret[n].size() << endl;
                std::cout << endl;
            }
        }


        double avg_time = accumulate(times.begin(), times.end(), 0.0) / N;
        double sum_sq_diff = 0.0;
        for (const auto& t : times) {
        sum_sq_diff += (t - avg_time) * (t - avg_time);
         }
        double stddev_time = sqrt(sum_sq_diff / N);

        std::cout << "Results for " << policy_name << endl;
        std::cout << "Mean final regret is " << accumulate(Regret.back().begin(), Regret.back().end(), 0.0) / Regret.back().size() << endl;
        std::cout << "Mean number of restarts is " << accumulate(Restarts.begin(), Restarts.end(), 0.0) / Restarts.size() << endl;

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        std::cout << "Elapsed time is " << elapsed.count() << " seconds" << endl;


        std::cout << "Mean time per experiment: " << avg_time << " seconds" << endl;
        std::cout << "Standard deviation of time: " << stddev_time << " seconds" << endl;

        if (saveres) {
            #pragma omp critical
            save_results(name, Regret, Restarts, avg_time, stddev_time);
        }
    }
}




int main(int argc, char* argv[]) {
    int T = std::stoi(argv[1]);
    double xi = std::stod(argv[2]);
    std::string PB = argv[3];

    bool Save = true;
    vector<int> N_list = {4000}; // Assuming N_list is fixed or internal
    double tweak = 1.0;
    int K = 5;
    double Deltamin = 0.1;
    double Deltamax = 0.4;
    double Deltachange = 0.1;
    double explor_freq = 0.3;

    for (const auto &N : N_list) {

        string folder_name = "./cpp_results_" + PB + "_" + to_string(T) + "_" + to_string(N) + "/";
        
        // Check if the folder exists, if not create it
        if (access(folder_name.c_str(), F_OK) == -1) {
            mkdir(folder_name.c_str(), 0777);
        }

        std::cout << "Folder created: " << folder_name << endl;

        std::cout << "Arms: " << K << ", Horizon: " << T << ", tweak: " << tweak << ", xi: " << xi << ", PB: " << PB << endl;
        auto start = chrono::high_resolution_clock::now();
        BayesianExpes(N, K, T, xi, Deltamin, Deltamax, tweak, explor_freq, Deltachange, Save, PB, folder_name);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        std::cout << "Total elapsed time is " << elapsed.count() << " seconds" << endl;
    }
    return 0;
}




