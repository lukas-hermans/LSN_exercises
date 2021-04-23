#include <vector>
#include <functional>
#include <cmath>
#include <iostream>
#include "random.h"
#include "metropolis.h"

// rnd: Random instance for random number generation
Metropolis::Metropolis(Random &rnd) : rnd(rnd)
{
}

Metropolis::~Metropolis()
{
}

// Set starting position.
void Metropolis::set_x0(std::vector<double> x0)
{
    x.push_back(x0);
}

// Calculate one step of Metroplis algorithm.
void Metropolis::do_step()
{
    std::vector<double> x_trial = T(rnd, x.back()); // trial for new step

    double alpha = min(1.0, p(x_trial) / p(x.back())); // acceptance probability (simplified version!)
    double r = rnd.Rannyu();                           // random number between 0 and 1

    (r <= alpha) ? (x.push_back(x_trial)) : (x.push_back(x.back())); // accept or reject trial position
}

void Metropolis::do_many_steps(int M_throws, int N_blocks)
{
    double sum;
    int L = M_throws / N_blocks; // number of throws per block
    for (int i = 0; i < N_blocks; i++)
    {
        sum = 0;
        for (int j = 0; j < L; j++)
        {
            this->do_step();
            sum += sample_func(x.back());
        }
        sum /= L;
        samples_block.push_back(sum);
    }
}

// std::vector<double> x = T(rnd, x_in); // calculate trial for next position of Metropolis algorithm

// double alpha = min(1.0, p(x) / p(x_in)); // acceptance probability (simplified version!)
// double r = rnd.Rannyu();                 // random number between 0 and 1

// return (r <= alpha) ? x : x_in; // accept or reject trial position
