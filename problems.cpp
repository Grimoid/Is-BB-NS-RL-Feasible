#include "problems.h"
#include "globals.h"
#include <vector>
#include <tuple>
#include <random>
#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace std;





int max_index(const vector<vector<double>>& matrix, int col) {
    int max_idx = 0;
    double max_val = matrix[0][col];
    for (int i = 1; i < matrix.size(); ++i) {
        if (matrix[i][col] > max_val) {
            max_val = matrix[i][col];
            max_idx = i;
        }
    }
    return max_idx;
}

vector<int> create_change_points(int T, double rho) {
    vector<int> change_points;
    int nbreaks = 0;
    while (nbreaks < 1) {
        change_points.clear();
        for (int i = 1; i < T; ++i) {
            if (distribution(generator) <= rho) change_points.push_back(i);
        }
        nbreaks = change_points.size();
    }
    change_points.insert(change_points.begin(), 0);
    change_points.push_back(T);
    return change_points;
}


tuple<vector<vector<double>>, vector<int>> WorstProblemGeo(int T, int K, double ksi, double Deltamin, double Deltamax, bool worst_case) {
    double rho = pow(T, -ksi);
    vector<int> change_points = create_change_points(T, rho);
    int nbreaks = change_points.size() - 1;

    vector<vector<double>> MeansMatrix(K, vector<double>(nbreaks, 0.0));
    for (int a = 0; a < K; ++a) {
        MeansMatrix[a][0] = 0.3+(distribution(generator)) * 0.005;
    }


    for (int k = 1; k < nbreaks; ++k) {

        int worst_arm=0;
        double worst_value=100000000;

        int best_arm=0;
        double best_value=-10;
        for(int a=0;a<K;++a){
            if(MeansMatrix[a][k-1]<=worst_value){
                worst_arm=a;
                worst_value=MeansMatrix[a][k-1];
            }
            if(MeansMatrix[a][k-1]>=best_value){
            best_arm=a;
            best_value=MeansMatrix[a][k-1];
            }
        }   
       

        double max_mean_possible = 1.0;  
        double min_mean_possible = 0.0;  
        if(best_value<1){


        double temp_val=best_value+(distribution(generator)) * 0.05;
        while(temp_val>1){

            temp_val=best_value+(distribution(generator)) * 0.05;
            if(temp_val>=0.99){best_value=1;break;}
        }
        MeansMatrix[worst_arm][k] = temp_val;

        for (int a = 0; a < K; ++a) {
            if (a == worst_arm) continue;
            MeansMatrix[a][k] = MeansMatrix[a][k-1];
            
        }
        }
        if(best_value==1){
            int worst_arm = std::min_element(MeansMatrix.begin(), MeansMatrix.end(), 
                                     [k](const std::vector<double>& a, const std::vector<double>& b) {
                                         return a[k-1] < b[k-1];
                                     }) - MeansMatrix.begin();

            int best_arm = std::max_element(MeansMatrix.begin(), MeansMatrix.end(), 
                                            [k](const std::vector<double>& a, const std::vector<double>& b) {
                                                return a[k-1] < b[k-1];
                                            }) - MeansMatrix.begin();

            std::vector<double> vals(K);
            for (int i = 0; i < K; ++i) {
                vals[i] = distribution(generator)*0.005+  0.3;
            }

            std::sort(vals.begin(), vals.end());

            for (int a = 0; a < K; ++a) {
                if (vals[a] == *std::min_element(vals.begin(), vals.end())) {
                    MeansMatrix[best_arm][k] = *std::min_element(vals.begin(), vals.end());
                } else if (vals[a] == *std::max_element(vals.begin(), vals.end())) {
                    MeansMatrix[worst_arm][k] = *std::max_element(vals.begin(), vals.end());
                } else {
                    MeansMatrix[a][k] = vals[a];
                }
            }
            
        }
    }
    return ProblemGeneral(MeansMatrix, change_points);
}

tuple<vector<vector<double>>, vector<int>> WorstProblemGeoUnif(int T, int K, double ksi) {
    double rho = pow(T, -ksi);
    int nbreaks=0;
    vector<int> change_points;

    while(nbreaks<rho*T){
        change_points = create_change_points(T, rho);
        nbreaks = change_points.size()-1;
    }
    nbreaks = change_points.size() - 1;

    vector<vector<double>> MeansMatrix(K, vector<double>(nbreaks, 0.0));
    for (int a = 0; a < K; ++a) {
        MeansMatrix[a][0] = 0.3+(distribution(generator)) * 0.0045+0.0005;
    }




    for (int k = 1; k < nbreaks; ++k) {

        int worst_arm=0;
        double worst_value=100000000;

        int best_arm=0;
        double best_value=-10;
        for(int a=0;a<K;++a){
            if(MeansMatrix[a][k-1]<=worst_value){
                worst_arm=a;
                worst_value=MeansMatrix[a][k-1];
            }
            if(MeansMatrix[a][k-1]>=best_value){
            best_arm=a;
            best_value=MeansMatrix[a][k-1];
            }
        }   
       

        double max_mean_possible = 1.0;  
        double min_mean_possible = 0.0;  
        if(best_value<1){

        double temp_val=best_value+(distribution(generator)) * 0.045+0.005;
        while(temp_val>1){

            temp_val=best_value+(distribution(generator)) * 0.045+0.005;
            if(temp_val>=0.99){best_value=1;break;}
        }
        MeansMatrix[worst_arm][k] = temp_val;

        for (int a = 0; a < K; ++a) {
            if (a == worst_arm) continue;
            MeansMatrix[a][k] = MeansMatrix[a][k-1];
            
        }
        }
        if(best_value==1){
            int worst_arm = std::min_element(MeansMatrix.begin(), MeansMatrix.end(), 
                                     [k](const std::vector<double>& a, const std::vector<double>& b) {
                                         return a[k-1] < b[k-1];
                                     }) - MeansMatrix.begin();

            int best_arm = std::max_element(MeansMatrix.begin(), MeansMatrix.end(), 
                                            [k](const std::vector<double>& a, const std::vector<double>& b) {
                                                return a[k-1] < b[k-1];
                                            }) - MeansMatrix.begin();

            std::vector<double> vals(K);
            for (int i = 0; i < K; ++i) {
                vals[i] = distribution(generator)*0.0045+0.0005+  0.3;
            }

            std::sort(vals.begin(), vals.end());

            for (int a = 0; a < K; ++a) {
                if (vals[a] == *std::min_element(vals.begin(), vals.end())) {
                    MeansMatrix[best_arm][k] = *std::min_element(vals.begin(), vals.end());
                } else if (vals[a] == *std::max_element(vals.begin(), vals.end())) {
                    MeansMatrix[worst_arm][k] = *std::max_element(vals.begin(), vals.end());
                } else {
                    MeansMatrix[a][k] = vals[a];
                }
            }

        }
    }
    return ProblemUnif(T, MeansMatrix);
}




tuple<vector<vector<double>>, vector<int>> ProblemGeo(int T, int K, double ksi, double Deltamin, double Deltamax, bool worst_case) {
    double rho = pow(T, -ksi);
    vector<int> change_points = create_change_points(T, rho);
    int nbreaks = change_points.size() - 1;

    vector<vector<double>> MeansMatrix(K, vector<double>(nbreaks, 0.0));
    for (int a = 0; a < K; ++a) {
        MeansMatrix[a][0] = (distribution(generator) + distribution(generator)) * 0.5;
    }

    for (int k = 1; k < nbreaks; ++k) {
        vector<int> changing_arms(K);
        iota(changing_arms.begin(), changing_arms.end(), 0);
        shuffle(changing_arms.begin(), changing_arms.end(), generator);
        int result_length = uniform_int_distribution<int>(1, K)(generator);
        changing_arms.resize(result_length);

        for (int a = 0; a < K; ++a) {
            double prevmean = MeansMatrix[a][k - 1];
            if (find(changing_arms.begin(), changing_arms.end(), a) != changing_arms.end()) {
                double gap = Deltamin + (Deltamax - Deltamin) * distribution(generator);
                if (distribution(generator) <= 0.5) gap = -gap;
                MeansMatrix[a][k] = (prevmean + gap > 1 || prevmean + gap < 0) ? prevmean - gap : prevmean + gap;
            } else {
                MeansMatrix[a][k] = prevmean;
            }
        }
    }
    return ProblemGeneral(MeansMatrix, change_points);
}

void PrintMeans(const vector<vector<double>>& MeansMatrix) {
    for (size_t i = 0; i < MeansMatrix.size(); ++i) {
        cout << "Arm " << i + 1 << ": ";
        for (size_t j = 0; j < MeansMatrix[i].size(); ++j) {
            cout << MeansMatrix[i][j] << " ";
        }
        cout << endl;
    }
}


tuple<vector<vector<double>>, vector<int>> ProblemGeoUnif(int T, int K, double ksi, double Deltamin, double Deltamax) {
    double rho = pow(T, -ksi);
    int nbreaks=0;
    vector<int> change_points;

    while(nbreaks<rho*T){
        change_points = create_change_points(T, rho);
        nbreaks = change_points.size()-1;
    }
    nbreaks = change_points.size() - 1;
    vector<vector<double>> MeansMatrix(K, vector<double>(nbreaks, 0.0));
    for (int a = 0; a < K; ++a) {
        MeansMatrix[a][0] = (distribution(generator) + distribution(generator)) * 0.5;
    }

    for (int k = 1; k < nbreaks; ++k) {
        vector<int> changing_arms(K);
        iota(changing_arms.begin(), changing_arms.end(), 0);
        shuffle(changing_arms.begin(), changing_arms.end(), generator);
        int result_length = uniform_int_distribution<int>(1, K)(generator);
        changing_arms.resize(result_length);

        for (int a = 0; a < K; ++a) {
            double prevmean = MeansMatrix[a][k - 1];
            if (find(changing_arms.begin(), changing_arms.end(), a) != changing_arms.end()) {
                double gap = Deltamin + (Deltamax - Deltamin) * distribution(generator);
                if (distribution(generator) <= 0.5) gap = -gap;
                MeansMatrix[a][k] = (prevmean + gap > 1 || prevmean + gap < 0) ? prevmean - gap : prevmean + gap;
            } else {
                MeansMatrix[a][k] = prevmean;
            }
        }
    }

    return ProblemUnif(T, MeansMatrix);
}







tuple<vector<vector<double>>, vector<int>> ProblemUnif(int T, const vector<vector<double>>& MeansMatrix) {
    int K = MeansMatrix.size();
    int Episodes = MeansMatrix[0].size();
    int nbBreaks = Episodes - 1;
    vector<vector<double>> Problem(K, vector<double>(T, 0.0));
    int part = round(T / (nbBreaks + 1));
    for (int c = 0; c <= nbBreaks; ++c) {
        for (int arm = 0; arm < K; ++arm) {
            fill(Problem[arm].begin() + c * part, Problem[arm].begin() + (c + 1) * part, MeansMatrix[arm][c]);
        }
    }
    if (Episodes * part < T) {
        for (int arm = 0; arm < K; ++arm) {
            fill(Problem[arm].begin() + Episodes * part, Problem[arm].end(), MeansMatrix[arm][Episodes - 1]);
        }
    }
    vector<int> BreakPoints;
    for (int i = 1; i <= nbBreaks; ++i) {
        BreakPoints.push_back(i * part);
    }
    return make_tuple(Problem, BreakPoints);
}



void PlotProblem(const vector<vector<double>>& Pb) {
    int K = Pb.size();
    int T = Pb[0].size();
    for (int k = 0; k < K; ++k) {
        cout << "Arm " << k << ": ";
        for (int t = 0; t < T; ++t) {
            cout << Pb[k][t] << " ";
        }
        cout << endl;
    }
}


tuple<vector<vector<double>>, vector<int>> ProblemGeneral(const vector<vector<double>>& MeansMatrix, const vector<int>& BreakPoints) {
    int L = BreakPoints.size();
    int K = MeansMatrix.size();
    int l = MeansMatrix[0].size();

    if (L > (l + 1)) {
        throw std::invalid_argument("wrong dimensions");
    } else {
        int T = BreakPoints[L - 1];
        vector<vector<double>> Problem(K, vector<double>(T, 0.0));
        
        for (int k = 0; k < L - 1; ++k) {
            for (int arm = 0; arm < K; ++arm) {
                int start = BreakPoints[k];
                int end = BreakPoints[k + 1];
                for (int t = start; t < end; ++t) {
                    Problem[arm][t] = MeansMatrix[arm][k];
                }
            }
        }

        vector<int> newBreakPoints(BreakPoints.begin() + 1, BreakPoints.end() - 1);
        return make_tuple(Problem, newBreakPoints);
    }
}