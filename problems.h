#ifndef PROBLEMS_H
#define PROBLEMS_H

#include <vector>
#include <tuple>
#include <string>

std::vector<int> create_change_points(int T, double rho);



std::tuple<std::vector<std::vector<double>>, std::vector<int>> WorstProblemGeoUnif(int T, int K, double ksi=0.85);
std::tuple<std::vector<std::vector<double>>, std::vector<int>> ProblemGeoUnif(int T, int K, double ksi, double Deltamin=0.1, double Deltamax=0.4);

std::tuple<std::vector<std::vector<double>>, std::vector<int>> ProblemGeo(int T, int K, double ksi=0.85, double Deltamin=0.1, double Deltamax=0.4, bool worst_case=false);

std::tuple<std::vector<std::vector<double>>, std::vector<int>> WorstProblemGeo(int T, int K, double ksi=0.85, double Deltamin=0.1, double Deltamax=0.4, bool worst_case=false);



std::tuple<std::vector<std::vector<double>>, std::vector<int>> ProblemUnif(int T, const std::vector<std::vector<double>>& MeansMatrix);
std::tuple<std::vector<std::vector<double>>, std::vector<int>> ProblemGeneral(const std::vector<std::vector<double>>& MeansMatrix, const std::vector<int>& BreakPoints);
void PlotProblem(const std::vector<std::vector<double>>& Pb);

#endif
