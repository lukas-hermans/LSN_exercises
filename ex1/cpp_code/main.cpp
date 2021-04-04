#include <iostream>
#include <fstream>
#include "random/random.h"

int main()
{

    int M = 100000; // number of throws
    int N = 100;    // number of blocks
    int L = M / N;  // number of throws per block

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for pseudo-random number generator
    Random rnd = generate_rnd(seed);        // prepared instance for pseudo-random number generation

    double sum;
    double A[N];
    double A_sq[N];
    for (int i = 0; i < N; i++)
    {
        sum = 0;
        for (int j = 0; j < L; j++)
        {
            sum += rnd.Rannyu();
        }
        A[i] = sum / L;
        A_sq[i] = A[i] * A[i];
    }

    // double A_mean[N];
    // double rolling_mean;
    // double A_sigma[N];
    // double rolling_sigma;
    // for (int i = 0; i < N; i++)
    // {
    //     rolling_mean
    // }

    std::ofstream file("../data/01.1_r_vs_M.txt");
    file << "M, r" << std::endl;
    for (int i = 0; i < N; i++)
    {
        file << A[i] << ", " << A_sq[i] << std::endl;
    }

    file.close();

    return 0;
}