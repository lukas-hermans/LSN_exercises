#include <vector>
#include <fstream>
#include <iostream>
#include "random/random.h"

int main()
{
    std::vector<int> N = {1, 2, 10, 100}; // number of terms in sum of random variables
    int N_repetitions = 10000;            // number of repetitions of sum calculation (=> in total: N_repetitions * N draws)

    // parameters for sampling from exp. and Lorentzian distribution
    double lambda = 1;
    double mu = 0;
    double gamma = 1;

    int seed[4] = {4, 13, 9, 17}; // seed for Rannyu generator
    Random rnd = Random(seed);    // prepared instance for Rannyu generator

    // create output files and add header
    std::ofstream file_uniform("../data/uniform.txt"), file_exponential("../data/exponential.txt"), file_lorentzian("../data/lorentzian.txt");
    for (int n : N)
    {
        file_uniform << "N=" << n << ", ";
        file_exponential << "N=" << n << ", ";
        file_lorentzian << "N=" << n << ", ";
    }

    // write sum of random variables to corresponding files
    double sum_uniform, sum_exponential, sum_lorentzian;
    int counter; // count size of sum
    for (int i = 0; i < N_repetitions; i++)
    {
        file_uniform << std::endl;
        file_exponential << std::endl;
        file_lorentzian << std::endl;

        sum_uniform = sum_exponential = sum_lorentzian = 0;
        counter = 0;
        for (int n : N) // loop over elements in vector N
        {
            while (counter < n) // sum of random variables with n terms (n are elements of vector N)
            {
                sum_uniform += rnd.Rannyu();
                sum_exponential += rnd.exponential_draw(lambda);
                sum_lorentzian += rnd.lorentzian_draw(gamma, mu);
                counter++;
            }
            file_uniform << sum_uniform / n << ", ";
            file_exponential << sum_exponential / n << ", ";
            file_lorentzian << sum_lorentzian / n << ", ";
        }
    }

    return 0;
}