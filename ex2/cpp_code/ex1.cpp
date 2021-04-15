#include <vector>
#include "random/random.h"
#include "tools/tools.h"

int main()
{
    int M = 100000; // number of throws
    int N = 100;    // number of blocks
    int L = M / N;  // number of throws per block

    int seed[4] = {4, 13, 92, 17}; // seed for Rannyu generator
    Random rnd = Random(seed);     // instance of random number generator

    std::vector<double> I_uniform = rnd.blocking_method_uniform(N, L);
    std::vector<double> I_p = rnd.blocking_method_p(N, L);

    std::vector<std::vector<double>> I_uniform_prog = tools::prog_mean(I_uniform);
    std::vector<std::vector<double>> I_p_prog = tools::prog_mean(I_p);

    tools::write_data(tools::vec_arange(L, M, L), I_uniform_prog[0], I_uniform_prog[1], "../data/I_uniform_prog.txt", "I, I(M), I_error(M)");
    tools::write_data(tools::vec_arange(L, M, L), I_p_prog[0], I_p_prog[1], "../data/I_p_prog.txt", "I, I(M), I_error(M)");

    return 0;
}