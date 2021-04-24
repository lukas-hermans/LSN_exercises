#include <vector>
#include <functional>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../../tools/random.h"
#include "../../tools/blocking.h"
#include "../../tools/metropolis.h"

// TO DO: progressive mean + error, output exp. value and positions, clear method

std::vector<double> T_uni(Random &rnd, std::vector<double> x_in, std::vector<double> x_min, std::vector<double> x_max);
std::vector<double> T_gaussian(Random &rnd, std::vector<double> x_in, std::vector<double> mu, std::vector<double> sigma);
double pdf_100(std::vector<double> x);
double pdf_210(std::vector<double> x);
double r(std::vector<double> x);

int main()
{
    int M = 1000000; // number of throws
    int N = 100;     // number of blocks

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed);              // prepared instance for Rannyu generator

    std::vector<double> x0 = {0, 0, 0}; // starting position for Metropolis algorithm

    // uniform trial transition probability (set to achieve approx. 50% acceptance rate)
    std::vector<double> x_min_100 = {-1.2, -1.2, -1.2};
    std::vector<double> x_max_100 = {1.2, 1.2, 1.2};
    std::vector<double> x_min_210 = {-3, -3, -3};
    std::vector<double> x_max_210 = {3, 3, 3};

    // gaussian trial transition probability (set to achieve approx. 50% acceptance rate)
    std::vector<double> mu_100 = {0, 0, 0};
    std::vector<double> sigma_100 = {0.75, 0.75, 0.75};
    std::vector<double> mu_210 = {0, 0, 0};
    std::vector<double> sigma_210 = {1.9, 1.9, 1.9};

    Metropolis metropolis(rnd); // instance for Metropolis algorithm
    metropolis.save_step = 50;

    // ************************** //
    // *** uniform transition *** //
    // ************************** //
    metropolis.set_x0(x0);
    metropolis.sample_func = r;

    // ground state
    metropolis.T = std::bind(T_uni, std::placeholders::_1, std::placeholders::_2, x_min_100, x_max_100);
    metropolis.p = pdf_100;
    metropolis.do_many_steps(M, N);
    std::cout << "acceptance rate for uni_100: " << metropolis.acc / (double)M << std::endl;
    metropolis.write_x_history("../data/uni_100.xyz");
    blocking::write_data(metropolis.blocking_data, "../data/uni_100.txt", "M, r_mean(M), r_error(M)");
    metropolis.clean_cache();

    // excited state
    metropolis.T = std::bind(T_uni, std::placeholders::_1, std::placeholders::_2, x_min_210, x_max_210);
    metropolis.p = pdf_210;
    metropolis.do_many_steps(M, N);
    std::cout << "acceptance rate for uni_210: " << metropolis.acc / (double)M << std::endl;
    metropolis.write_x_history("../data/uni_210.xyz");
    blocking::write_data(metropolis.blocking_data, "../data/uni_210.txt", "M, r_mean(M), r_error(M)");
    metropolis.clean_cache();

    // *************************** //
    // *** gaussian transition *** //
    // *************************** //

    // ground state
    metropolis.T = std::bind(T_gaussian, std::placeholders::_1, std::placeholders::_2, mu_100, sigma_100);
    metropolis.p = pdf_100;
    metropolis.do_many_steps(M, N);
    std::cout << "acceptance rate for gaussian_100: " << metropolis.acc / (double)M << std::endl;
    metropolis.write_x_history("../data/gaussian_100.xyz");
    blocking::write_data(metropolis.blocking_data, "../data/gaussian_100.txt", "M, r_mean(M), r_error(M)");
    metropolis.clean_cache();

    // excited state
    metropolis.T = std::bind(T_gaussian, std::placeholders::_1, std::placeholders::_2, mu_210, sigma_210);
    metropolis.p = pdf_210;
    metropolis.do_many_steps(M, N);
    std::cout << "acceptance rate for gaussian_210: " << metropolis.acc / (double)M << std::endl;
    metropolis.write_x_history("../data/gaussian_210.xyz");
    blocking::write_data(metropolis.blocking_data, "../data/gaussian_210.txt", "M, r_mean(M), r_error(M)");
    metropolis.clean_cache();

    return 0;
}

// Uniform trial transition probability given a vector x_in.
// Specfied by x_min and x_max (vectors; for every direction).
std::vector<double> T_uni(Random &rnd, std::vector<double> x_in, std::vector<double> x_min, std::vector<double> x_max)
{
    std::vector<double> x;
    for (int i = 0; i < (int)x_in.size(); i++)
    {
        x.push_back(x_in[i] + rnd.Rannyu(x_min[i], x_max[i]));
    }
    return x;
}

// Gaussian trial transition probability given a vector x_in.
// Specfied by mu and sigma (vectors; for every direction).
std::vector<double> T_gaussian(Random &rnd, std::vector<double> x_in, std::vector<double> mu, std::vector<double> sigma)
{
    std::vector<double> x;
    for (int i = 0; i < (int)x_in.size(); i++)
    {
        x.push_back(x_in[i] + rnd.Gauss(mu[i], sigma[i]));
    }
    return x;
}

// Calculates value of pdf of hydrogen orbital 100 at position x.
double pdf_100(std::vector<double> x)
{
    return 1 / M_PI * exp(-2 * r(x));
}

// Calculates value of pdf of hydrogen orbital 210 at position x.
double pdf_210(std::vector<double> x)
{
    return 1 / (32 * M_PI) * x[2] * x[2] * exp(-r(x));
}

// Calculates radius from given position vector.
double r(std::vector<double> x)
{
    return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}