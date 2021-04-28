#include <vector>
#include <cmath>
#include "../../tools/random.h"
#include "../../tools/blocking.h"

double g_uniform(Random &rnd);
double f_importance(Random &rnd);

int main()
{
    int M = 100000; // number of throws
    int N = 100;    // number of blocks
    int L = M / N;  // number of throws per block

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed);              // instance of random number generator

    std::vector<std::vector<double>> I_uniform_prog = blocking::prog_mean(g_uniform, rnd, N, L);
    std::vector<std::vector<double>> I_importance_prog = blocking::prog_mean(f_importance, rnd, N, L);

    blocking::write_data(I_uniform_prog, "../data/I_uniform_prog.txt", "M, I_prog(M), error(M)");
    blocking::write_data(I_importance_prog, "../data/I_importance_prog.txt", "M, I_prog(M), error(M)");

    return 0;
}

// Samples g(x) = pi/2*cos(pi*x/2) where x is distributed uniformly on [0, 1).
double g_uniform(Random &rnd)
{
    double x = rnd.Rannyu(); // sample from p(x) = 1
    return M_PI / 2 * cos(M_PI * x / 2);
}

// Samples f(x) = pi/2*cos(pi*x/2) / d(x) where x is distributed according to d(x) = 2*(1-x) on [0, 1).
double f_importance(Random &rnd)
{
    double x = rnd.importance_draw(); // sample from d(x) = 2*(1-x)
    return M_PI / 2 * cos(M_PI * x / 2) / (2 * (1 - x));
}