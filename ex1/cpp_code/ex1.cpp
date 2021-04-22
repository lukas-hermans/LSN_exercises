#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <functional>
#include "tools/random.h"
#include "tools/blocking.h"

double rnd_variable(Random &rnd);
double rnd_variance(Random &rnd);
void count(std::vector<int> &counter, double draw, int N_intervals, double int_length);
double chi2(std::vector<int> counter, int N_intervals, int N_draws);

int main()
{
    // input ex. 01.1.1 & 01.1.2
    int M = 150000; // number of throws
    int N = 100;    // number of blocks
    int L = M / N;  // number of throws per block

    // input ex. 01.1.3
    int N_repetitions = 100;               // number of repetitions for calculation of chi2
    int N_draws = 10000;                   // number of draws per repetition
    int N_intervals = 100;                 // number of subintervals of [0, 1)
    double int_length = 1.0 / N_intervals; // length of each subinterval

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed);              // prepared instance for Rannyu generator

    // calculate progressive mean & error for Rannyu generator as well as its variance & error
    std::vector<std::vector<double>> r_prog = blocking::prog_mean(rnd_variable, rnd, N, L);
    std::vector<std::vector<double>> r_var_prog = blocking::prog_mean(rnd_variance, rnd, N, L);

    // calculate chi2 for number of repitions specified by "N_repetitions"
    std::vector<double> chi2_vec;
    std::vector<int> counter(N_intervals);
    double draw;
    for (int i = 0; i < N_repetitions; i++)
    {
        std::fill(counter.begin(), counter.end(), 0); // set all entries of counter vector to zero

        for (int j = 0; j < N_draws; j++)
        {
            draw = rnd.Rannyu();
            count(counter, draw, N_intervals, int_length);
        }

        chi2_vec.push_back(chi2(counter, N_intervals, N_draws));
    }

    rnd.SaveSeed();

    blocking::write_data(r_prog, "../data/r_vs_M.txt", "M, r_mean(M), r_mean_error(M)");
    blocking::write_data(r_var_prog, "../data/r_var_vs_M.txt", "M, r_var(M), r_var_error(M)");
    blocking::write_data(chi2_vec, "../data/chi2.txt", "chi2");

    return 0;
}

// Samples a random variable r using Rannyu generator.
double rnd_variable(Random &rnd)
{
    double r = rnd.Rannyu();
    return r;
}

// Samples (r - 0.5)^2 using Rannyu generator.
double rnd_variance(Random &rnd)
{
    return pow(rnd.Rannyu() - 0.5, 2);
}

// Searches for the subinterval in which the draw lies and increases the corresponding entry of the counter vector by 1.
// Uses a sequential search algorithm.
// Note that a binary search would be more efficient, but for the used number of subintervals the computational difference is negligible.
void count(std::vector<int> &counter, double draw, int N_intervals, double int_length)
{
    for (int i = 0; i < N_intervals; i++)
    {
        if (draw >= int_length * i && draw < int_length * (i + 1))
        {
            counter[i] += 1;
            return;
        }
    }
}

// Computes the value of chi2 given an observed count of N_draws random numbers over the N_intervals subintervals.
double chi2(std::vector<int> counter, int N_intervals, int N_draws)
{
    int N_exp = N_draws / N_intervals; // expected number of draws in each subinterval (if uniform)
    double chi2 = 0;
    for (int i = 0; i < N_intervals; i++) // sum over all subintervals to compute chi2
    {
        chi2 += pow(counter[i] - N_exp, 2) / (double)N_exp;
    }

    return chi2;
}