#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <vector>

inline double AVERAGEX(const std::vector<std::vector<double>> &Var, int i, int j){
    return 0.5 * (Var[j][i] + Var[j][i + 1]);
}

inline double AVERAGEY(const std::vector<std::vector<double>> &Var, int i, int j){
    return 0.5 * (Var[j][i] + Var[j + 1][i]);
}

inline double AVERAGEXY(const std::vector<std::vector<double>> &Var, int i, int j){
    return 0.25 * (Var[j][i] + Var[j][i + 1] + Var[j + 1][i] + Var[j + 1][i + 1]);
}

inline double LAPLACIAN(const std::vector<std::vector<double>> &Var, int i, int j, double Dx, double Dy){
    return (Var[j][i + 1] - 2 * Var[j][i] + Var[j][i - 1]) / (Dx * Dx) + (Var[j + 1][i] - 2 * Var[j][i] + Var[j - 1][i]) / (Dy * Dy);
}



#endif //UTILS_HPP_