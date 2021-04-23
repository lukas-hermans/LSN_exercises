#include <vector>
#include <functional>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../../tools/random.h"
#include "../../tools/blocking.h"
#include "../../tools/metropolis.h"

// TO DO: progressive mean + error, output exp. value and positions

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

    // uniform trial transition probability
    std::vector<double> x_min = {-0.5, -0.5, -0.5};
    std::vector<double> x_max = {0.5, 0.5, 0.5};

    // gaussian trial transition probability
    std::vector<double> mu = {0, 0, 0};
    std::vector<double> sigma = {1, 1, 1};

    Metropolis metropolis(rnd);
    metropolis.T = std::bind(T_uni, std::placeholders::_1, std::placeholders::_2, x_min, x_max);
    metropolis.p = pdf_210;
    metropolis.set_x0(x0);
    metropolis.sample_func = r;

    metropolis.do_many_steps(M, N);

    blocking::write_data({metropolis.samples_block}, "../data/test.txt", "dummy");

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
    double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    return 1 / M_PI * exp(-2 * r);
}

// Calculates value of pdf of hydrogen orbital 210 at position x.
double pdf_210(std::vector<double> x)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    return 1 / (32 * M_PI) * x[2] * x[2] * exp(-r);
}

// calculates radius from given position vector
double r(std::vector<double> x)
{
    return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}