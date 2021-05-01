#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <functional>
#include <cmath>
#include "../../tools/random.h"
#include "../../tools/blocking.h"

double C_direct(Random &rnd, double S_0, double T, double K, double r, double sigma);
double P_direct(Random &rnd, double S_0, double T, double K, double r, double sigma);
double C_discretized(Random &rnd, double N_steps, double S_0, double T, double K, double r, double sigma);
double P_discretized(Random &rnd, double N_steps, double S_0, double T, double K, double r, double sigma);

int main()
{
    double S_0 = 100;    // asset price at t=0
    double T = 1;        // delivery time
    double K = 100;      // strike price
    double r = 0.1;      // risk-free interest rate
    double sigma = 0.25; // volatility

    int N_steps = 100; // number of steps for discretized calculation (time of each step: T / N_steps)

    int M = 100000; // number of throws
    int N = 100;    // number of blocks
    int L = M / N;  // number of throws per block

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed);              // instance of random number generator

    std::vector<std::vector<double>> C_direct_prog = blocking::prog_mean(std::bind(C_direct, std::placeholders::_1, S_0, T, K, r, sigma), rnd, N, L);
    std::vector<std::vector<double>> P_direct_prog = blocking::prog_mean(std::bind(P_direct, std::placeholders::_1, S_0, T, K, r, sigma), rnd, N, L);
    std::vector<std::vector<double>> C_discretized_prog = blocking::prog_mean(std::bind(C_discretized, std::placeholders::_1, N_steps, S_0, T, K, r, sigma), rnd, N, L);
    std::vector<std::vector<double>> P_discretized_prog = blocking::prog_mean(std::bind(P_discretized, std::placeholders::_1, N_steps, S_0, T, K, r, sigma), rnd, N, L);

    blocking::write_data(C_direct_prog, "../data/C_direct_prog.txt", "M, C(M), C_error(M)");
    blocking::write_data(P_direct_prog, "../data/P_direct_prog.txt", "M, P(M), P_error(M)");
    blocking::write_data(C_discretized_prog, "../data/C_discretized_prog.txt", "M, C(M), C_error(M)");
    blocking::write_data(P_discretized_prog, "../data/P_discretized_prog.txt", "M, P(M), P_error(M)");

    return 0;
}

double C_direct(Random &rnd, double S_0, double T, double K, double r, double sigma)
{
    double Z = rnd.Gauss(0, 1);
    double S_T = S_0 * exp((r - pow(sigma, 2) / 2) * T + sigma * Z * sqrt(T));
    return exp(-r * T) * max(0.0, S_T - K);
}

double P_direct(Random &rnd, double S_0, double T, double K, double r, double sigma)
{
    double Z = rnd.Gauss(0, 1);
    double S_T = S_0 * exp((r - pow(sigma, 2) / 2) * T + sigma * Z * sqrt(T));
    return exp(-r * T) * max(0.0, K - S_T);
}

double C_discretized(Random &rnd, double N_steps, double S_0, double T, double K, double r, double sigma)
{
    double Z;
    double S_disc = S_0;
    for (int i = 0; i < N_steps; i++)
    {
        Z = rnd.Gauss(0, 1);
        S_disc = S_disc * exp((r - 0.5 * pow(sigma, 2)) * T / N_steps + sigma * Z * sqrt(T / N_steps));
    }
    double S_T = S_disc;

    return exp(-r * T) * max(0.0, S_T - K);
}

double P_discretized(Random &rnd, double N_steps, double S_0, double T, double K, double r, double sigma)
{
    double Z;
    double S_disc = S_0;
    for (int i = 0; i < N_steps; i++)
    {
        Z = rnd.Gauss(0, 1);
        S_disc = S_disc * exp((r - 0.5 * pow(sigma, 2)) * T / N_steps + sigma * Z * sqrt(T / N_steps));
    }
    double S_T = S_disc;

    return exp(-r * T) * max(0.0, K - S_T);
}
