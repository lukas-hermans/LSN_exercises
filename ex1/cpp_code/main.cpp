#include <iostream>
#include "random/random.h"

int main()
{
    int seed[4] = {0000, 0000, 0000, 0001}; // seed for pseudo-random number generator
    int M = 100000;                         // number of throws
    int N = 100;                            // number of blocks
    int L = M / N;                          // number of throws per block

    Random rnd = generate_rnd(seed);
    double a;
    for (int i = 0; i < 100; i++)
    {
        a = rnd.Rannyu();
    }
    return 0;
}