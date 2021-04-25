#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <fstream>
#include <iostream>
#include "random.h"
#include "blocking.h"
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
    this->x0 = x0;
    x_current = x0;
}

// Calculate one step of Metroplis algorithm.
void Metropolis::do_step()
{
    std::vector<double> x_trial = T(rnd, x_current); // trial for new step

    double alpha = min(1.0, p(x_trial) / p(x_current)); // acceptance probability (simplified version!)
    double r = rnd.Rannyu();                            // random number between 0 and 1

    if (r <= alpha) // accept trial
    {
        x_current = x_trial;
        if (is_equi == false)
        {
            if (step_count % save_step == 0)
            {
                x.push_back(x_current);
            }
            acc += 1;
        }
    }
    if (r > alpha) // reject trial
    {
        if (is_equi == false)
        {
            if (step_count % save_step == 0)
            {
                x.push_back(x_current);
            }
        }
    }

    if (is_equi == false)
    {
        step_count += 1;
    }
}

// Computes progressive mean and its error of sample_func using blocking method with N blocks of L steps of Metropolis algorithm.
// N_equi is the number of steps chosen to equilibrate the Metropolis algorithm.
// Same as prog_mean function in blocking.cpp but sample_func depends on memory of x (position history).
// Stores results in blocking_data vector.
void Metropolis::do_many_steps(int M_throws, int N_blocks, int N_equi)
{
    // equilibrate the Metropolis algorithm
    is_equi = true;
    for (int i = 0; i < N_equi; i++)
    {
        this->do_step();
    }
    is_equi = false;

    int L = M_throws / N_blocks; // number of throws per block
    std::vector<double> av, av2;

    // calculate av and av2 of given sample_func for each block
    double sum = 0;
    for (int i = 0; i < N_blocks; i++)
    {
        for (int j = 0; j < L; j++)
        {
            this->do_step();
            sum += sample_func(x_current);
        }
        sum /= L;

        av.push_back(sum);            // add sample value of current block
        av2.push_back(pow(av[i], 2)); // add squared sample value of current block

        sum = 0;
    }

    // calculate progressive mean and its error from av and av2
    double sum_prog = 0, sum2_prog = 0;
    for (int i = 0; i < N_blocks; i++)
    {
        for (int j = 0; j < i + 1; j++)
        {
            sum_prog += av[j];
            sum2_prog += av2[j];
        }
        blocking_data[0].push_back((i + 1) * L);
        blocking_data[1].push_back(sum_prog / (i + 1));
        blocking_data[2].push_back(blocking::error(blocking_data[1][i], sum2_prog / (i + 1), i + 1));

        sum_prog = sum2_prog = 0;
    }
}

// Write position history x to a file specified by path.
void Metropolis::write_x_history(std::string path)
{
    std::ofstream file(path);

    for (int i = 0; i < (int)x.size(); i++) // loop over all stored positions
    {
        for (int j = 0; j < (int)x[0].size(); j++) // loop over number of dimensions
        {
            file << x[i][j];
            if (j < (int)x[0].size() - 1)
            {
                file << ", ";
            }
        }
        file << std::endl;
    }
}

// Resets cache to start with a new Metropolis simulation.
void Metropolis::clean_cache()
{
    x.resize(1);
    x_current = x0;
    for (int i = 0; i < (int)blocking_data.size(); i++)
    {
        blocking_data[i].resize(0);
    }
    step_count = 0;
    acc = 0;
}
