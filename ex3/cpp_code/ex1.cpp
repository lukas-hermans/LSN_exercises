#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "random/random.h"
#include "tools/tools.h"

int main()
{
    double S_0 = 100;    // asset price at t=0
    double T = 1;        // delivery time
    double K = 100;      // strike price
    double r = 0.1;      // risk-free interest rate
    double sigma = 0.25; // volatility

    int N_steps = 100; // number of steps for discretized calculation (time of each step: T / N_steps)

    int M = 5E6;   // number of throws
    int N = 10000; // number of blocks
    int L = M / N; // number of throws per block

    int seed[4] = {4, 13, 9, 17}; // seed for Rannyu generator
    Random rnd = Random(seed);    // instance of random number generator

    std::vector<double> C_direct = rnd.blocking_method_direct("call", N, L, S_0, T, K, r, sigma);
    std::vector<double> P_direct = rnd.blocking_method_direct("put", N, L, S_0, T, K, r, sigma);
    std::vector<double> C_discretized = rnd.blocking_method_discretized("call", N, L, N_steps, S_0, T, K, r, sigma);
    std::vector<double> P_discretized = rnd.blocking_method_discretized("put", N, L, N_steps, S_0, T, K, r, sigma);

    std::vector<std::vector<double>> C_direct_prog = tools::prog_mean(C_direct);
    std::vector<std::vector<double>> P_direct_prog = tools::prog_mean(P_direct);
    std::vector<std::vector<double>> C_discretized_prog = tools::prog_mean(C_discretized);
    std::vector<std::vector<double>> P_discretized_prog = tools::prog_mean(P_discretized);

    tools::write_data(tools::vec_arange(L, M, L), C_direct_prog[0], C_direct_prog[1], "../data/C_direct_prog.txt", "M, C(M), C_error(M)");
    tools::write_data(tools::vec_arange(L, M, L), P_direct_prog[0], P_direct_prog[1], "../data/P_direct_prog.txt", "M, P(M), P_error(M)");
    tools::write_data(tools::vec_arange(L, M, L), C_discretized_prog[0], C_discretized_prog[1], "../data/C_discretized_prog.txt", "M, C(M), C_error(M)");
    tools::write_data(tools::vec_arange(L, M, L), P_discretized_prog[0], P_discretized_prog[1], "../data/P_discretized_prog.txt", "M, P(M), P_error(M)");

    return 0;
}