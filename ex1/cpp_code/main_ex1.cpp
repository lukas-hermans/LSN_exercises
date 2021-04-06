#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random/random.h"
#include "tools/tools.h"

int main()
{
    // input ex. 01.1 & 01.2
    int M = 100000; // number of throws
    int N = 100;    // number of blocks
    int L = M / N;  // number of throws per block

    // input ex. 01.3
    int N_intervals = 100;         // number of subintervals of [0, 1)
    double int_length = 1 / 100.0; // length of one interval
    int repeat = 100;              // number of repetition for calculation of chi2
    int n = 10000;                 // number of draws per repition

    int seed[4] = {4, 13, 9, 17}; // seed for Rannyu generator
    Random rnd = Random(seed);    // prepared instance for Rannyu generator

    // calculate random numbers and their variance for each block
    std::vector<double> r = rnd.blocking_method(N, L);
    std::vector<double> r_var = rnd.blocking_method_var(N, L);

    // calculate progressive means and their errors (as a function of the block number N)
    std::vector<std::vector<double>> r_prog = tools::prog_mean(r);
    std::vector<std::vector<double>> r_var_prog = tools::prog_mean(r_var);

    // calculate chi2 for number of repitions specified by "repeat"
    std::vector<double> chi2;
    std::vector<int> counter(N_intervals, 0.0);
    double draw;
    for (int i = 0; i < repeat; i++)
    {
        std::fill(counter.begin(), counter.end(), 0);
        for (int j = 0; j < n; j++)
        {
            draw = rnd.Rannyu();
            tools::count(counter, draw, N_intervals, int_length);
        }
        chi2.push_back(tools::chi2(counter, N_intervals, n));
    }

    rnd.SaveSeed();

    tools::write_data("../data/1_r_vs_N.txt", "N, r_mean(N), r_mean_std(N)", r_prog[0], r_prog[1]);
    tools::write_data("../data/1_var_vs_N.txt", "N, r_var(N), r_var_std(N)", r_var_prog[0], r_var_prog[1]);
    tools::write_data("../data/1_chi2.txt", "repeat, chi2, ignore", chi2);

    return 0;
}