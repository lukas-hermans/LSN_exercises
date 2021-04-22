#include <vector>
#include <fstream>
#include <string>
#include <functional>
#include "tools/random.h"

void write_dist(std::function<double(Random &)> dist_func, Random &rnd, int N_repetitions, std::vector<int> N_terms, std::string path);
double uni_variable(Random &rnd);
double exp_variable(Random &rnd, double lambda);
double lorentz_variable(Random &rnd, double gamma, double mu);

int main()
{
    std::vector<int> N_terms = {1, 2, 10, 100}; // number of terms in sum of random variables
    int N_repetitions = 10000;                  // number of repetitions of sum calculation (=> in total: N_repetitions * N draws)

    double lambda = 1;        // lambda parameter for exponential draws
    double mu = 0, gamma = 1; // mu and gamma parameter for Lorentzian draw

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd(seed);                       // prepared instance for Rannyu generator

    write_dist(uni_variable, rnd, N_repetitions, N_terms, "../data/uniform.txt");
    write_dist(std::bind(exp_variable, std::placeholders::_1, lambda), rnd, N_repetitions, N_terms, "../data/exponential.txt");
    write_dist(std::bind(lorentz_variable, std::placeholders::_1, gamma, mu), rnd, N_repetitions, N_terms, "../data/lorentzian.txt");

    return 0;
}

// Writes draws from a sum of N_terms from random distribution specified by dist_func for N_repetitions into a file specified by path.
void write_dist(std::function<double(Random &)> dist_func, Random &rnd, int N_repetitions, std::vector<int> N_terms, std::string path)
{
    std::ofstream file(path);

    // write header of file
    for (int n : N_terms)
    {
        file << "N_terms = " << n;
        if (n != N_terms.back()) // no comma after last column
        {
            file << ", ";
        }
    }

    // draw from distribution for all entries in N_terms with N_repetitions
    // each repetition is a line of the file
    double sum;
    for (int i = 0; i < N_repetitions; i++)
    {
        file << std::endl;

        for (int n : N_terms)
        {
            sum = 0;
            for (int j = 0; j < n; j++) // calculate sum of n terms
            {
                sum += dist_func(rnd);
            }
            file << sum / n;

            if (n != N_terms.back()) // no comma after last column
            {
                file << ", ";
            }
        }
    }
}

// Samples a uniformally distributed random variable.
double uni_variable(Random &rnd)
{
    return rnd.Rannyu();
}

// Samples an exponentially distributed random variable.
double exp_variable(Random &rnd, double lambda)
{
    return rnd.exponential_draw(lambda);
}

// Samples an Lorentzian distributed random variable.
double lorentz_variable(Random &rnd, double gamma, double mu)
{
    return rnd.lorentzian_draw(gamma, mu);
}