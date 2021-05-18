#include <vector>
#include <cmath>
#include <fstream>
#include "../../tools/random.h"
#include "../../tools/blocking.h"
#include "../../tools/metropolis.h"

std::vector<double> T_uni(Random &rnd, std::vector<double> x_in, std::vector<double> x_min, std::vector<double> x_max);
double pdf(std::vector<double> x, double mu, double sigma);
double sample_func(std::vector<double> x, double mu, double sigma);
std::vector<double> vec_linspace(double from, double to, int steps);

int main()
{
    int M = 250000;    // number of throws
    int N = 100;       // number of blocks
    int N_equi = 2000; // number of steps of Metropolis algorithm to equilibrate

    int nbins = 100; // number of bins for histogram
    int xlow = -5;   // starting point of histogram (w.r.t. x)
    int xup = 5;     // ending point of histogram (w.r.t. x)

    int seed[4] = {2780, 1540, 2047, 641}; // seed for Rannyu generator
    Random rnd = Random(seed);             // prepared instance for Rannyu generator

    std::vector<double> x0 = {0.0}; //starting point
    std::vector<double> x_min = {-2.5};
    std::vector<double> x_max = {2.5};

    Metropolis metropolis(rnd); // instance for Metropolis algorithm
    metropolis.T = std::bind(T_uni, std::placeholders::_1, std::placeholders::_2, x_min, x_max);
    metropolis.set_x0(x0);

    // // minimization of energy (in dependence on mu and sigma)
    // std::vector<double> mu = vec_linspace(0.6, 1.0, 200);
    // std::vector<double> sigma = vec_linspace(0.4, 0.8, 200);

    // std::ofstream file("../data/opt_params.txt");

    // file << "mu, sigma, E" << std::endl;

    // // compute energy for every combination of mu and sigma and print result
    // for (int i = 0; i < (int)mu.size(); i++)
    // {
    //     for (int j = 0; j < (int)sigma.size(); j++)
    //     {
    //         metropolis.clean_cache();

    //         metropolis.sample_func = std::bind(sample_func, std::placeholders::_1, mu[i], sigma[j]);
    //         metropolis.p = std::bind(pdf, std::placeholders::_1, mu[i], sigma[j]);

    //         metropolis.do_many_steps(M, N, N_equi);

    //         // print results for current parameters
    //         std::cout << "mu = " << mu[i] << ", sigma = " << sigma[j] << std::endl;
    //         std::cout << "acceptance rate: " << metropolis.acc / (double)M << std::endl;
    //         std::cout << "<H> = " << metropolis.blocking_data[1][N - 1] << std::endl;
    //         std::cout << std::endl;

    //         file << mu[i] << ", " << sigma[j] << ", " << metropolis.blocking_data[1][N - 1] << std::endl;
    //     }
    // }

    // ideal parameters
    double mu_ideal = 0.786;
    double sigma_ideal = 0.622;

    metropolis.save_step = 1; // used for histogram in results presentation

    metropolis.clean_cache();
    metropolis.sample_func = std::bind(sample_func, std::placeholders::_1, mu_ideal, sigma_ideal);
    metropolis.p = std::bind(pdf, std::placeholders::_1, mu_ideal, sigma_ideal);

    metropolis.do_many_steps(M, N, N_equi);
    std::cout << "acceptance rate: " << metropolis.acc / (double)M << std::endl;

    blocking::write_data(metropolis.blocking_data, "../data/E_vs_M.txt", "M, E_mean(M), E_error(M)");
    metropolis.write_hist("../data/hist.txt", nbins, xlow, xup, M, N);

    rnd.SaveSeed();

    return 0;
}

// Uniform trial transition probability given a vector x_in.
// Specfied by x_min and x_max (vectors; for every direction).
std::vector<double> T_uni(Random &rnd, std::vector<double> x_in, std::vector<double> x_min, std::vector<double> x_max)
{
    std::vector<double> x;
    x.push_back(x_in[0] + rnd.Rannyu(x_min[0], x_max[0]));

    return x;
}

// Probability density function of trial wavefunction.
double pdf(std::vector<double> x, double mu, double sigma)
{
    return pow(exp(-pow(x[0] - mu, 2) / (2 * pow(sigma, 2))) + exp(-pow(x[0] + mu, 2) / (2 * pow(sigma, 2))), 2);
}

// Function in the integral to sample total energy.
// Depends on trial wavefunction (and potential).
double sample_func(std::vector<double> x, double mu, double sigma)
{
    double psi = exp(-pow(x[0] - mu, 2) / (2 * pow(sigma, 2))) + exp(-pow(x[0] + mu, 2) / (2 * pow(sigma, 2)));
    double dpsi_dx2 = 1 / pow(sigma, 2) * (pow((x[0] - mu) / sigma, 2) - 1) * exp(-pow(x[0] - mu, 2) / (2 * pow(sigma, 2))) + 1 / pow(sigma, 2) * (pow((x[0] + mu) / sigma, 2) - 1) * exp(-pow(x[0] + mu, 2) / (2 * pow(sigma, 2)));
    double V = pow(x[0], 4) - 5 / 2.0 * pow(x[0], 2);

    return (-0.5 * dpsi_dx2 + V * psi) / psi;
}

std::vector<double> vec_linspace(double from, double to, int steps)
{
    std::vector<double> vec;

    for (int i = 0; i < steps; i++)
    {
        vec.push_back(from + i * (to - from) / steps);
    }

    return vec;
}