#include <vector>
#include <cmath>
#include <iostream>
#include <functional>
#include "../../tools/random.h"
#include "../../tools/blocking.h"

double r_sq(std::vector<double> xyz);
void discrete_step(Random &rnd, std::vector<double> &xyz);
void continuous_step(Random &rnd, std::vector<double> &xyz);
std::vector<std::vector<double>> do_walk(Random &rnd, std::function<void(Random &, std::vector<double> &)> step_func, int N_steps, int N_blocks, int L);

int main()
{
    int N_repetitions = 10000;        // number of repetitions of random walk
    int N_steps = 100;                // number of steps of random walker
    int N_blocks = 100;               // number of blocks used for calculation of uncertainties (blocking method!)
    int L = N_repetitions / N_blocks; // number of repetitions per block

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed);              // instance of random number generator

    std::vector<std::vector<double>> discrete_walk = do_walk(rnd, discrete_step, N_steps, N_blocks, L);
    std::vector<std::vector<double>> continuous_walk = do_walk(rnd, continuous_step, N_steps, N_blocks, L);

    blocking::write_data(discrete_walk, "../data/discrete_walk.txt", "step, r_mean, r_error");
    blocking::write_data(continuous_walk, "../data/continuous_walk.txt", "step, r_mean, r_error");

    return 0;
}

// Computes squared radius of given position vector.
double r_sq(std::vector<double> xyz)
{
    return xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
}

// Computes a discrete step on a 3D lattice (with lattice constant a=1) and changes xyz correspondingly.
void discrete_step(Random &rnd, std::vector<double> &xyz)
{
    int move_dir = floor(rnd.Rannyu(1, 7)); // random number in {1,2,3,4,5,6}

    if (move_dir == 1) // move in +x direction
        xyz[0] += 1;
    if (move_dir == 2) // -x
        xyz[0] -= 1;
    if (move_dir == 3) // +y
        xyz[1] += 1;
    if (move_dir == 4) // -y
        xyz[1] -= 1;
    if (move_dir == 5) // +z
        xyz[2] += 1;
    if (move_dir == 6) // -z
        xyz[2] -= 1;
}

// Computes a continuous step of length a=1 and changes xyz correspondingly.
void continuous_step(Random &rnd, std::vector<double> &xyz)
{
    double phi = rnd.Rannyu(0, 2 * M_PI);
    double theta = rnd.Rannyu(0, M_PI);

    xyz[0] += cos(phi) * sin(theta);
    xyz[1] += sin(phi) * sin(theta);
    xyz[2] += cos(theta);
}

// Computes mean radius after 1, 2, 3, ..., N_steps steps and its error using N_blocks with L repetitions for each block.
// Returns a vector that contains three nested vectors:
//      1. Number of steps: 1, 2, ..., N_steps
//      2. Mean radius
//      3. Error
std::vector<std::vector<double>> do_walk(Random &rnd, std::function<void(Random &, std::vector<double> &)> step_func, int N_steps, int N_blocks, int L)
{
    std::vector<double> xyz;                     // current position of walker
    std::vector<double> step_list = {0};         // vector with 1, 2, ..., N_steps
    std::vector<double> r2_block(N_steps);       // vector to store mean squared radius for current block
    std::vector<double> r_av(N_steps + 1, 0.0);  // vector to store root of mean squared radius
    std::vector<double> r_av2(N_steps + 1, 0.0); // vector to store mean squared radius

    // compute mean squared radius for 1, 2, ..., N_steps using N_blocks with L repetitions
    for (int i = 0; i < N_blocks; i++) // loop over all blocks
    {
        std::fill(r2_block.begin(), r2_block.end(), 0.0);

        //  compute mean squared radius for 1, 2, ..., N_steps for current block
        for (int j = 0; j < L; j++) // loop over repetitions for current block
        {
            xyz = {0, 0, 0};

            for (int k = 0; k < N_steps; k++) // do a random walk of in total N_steps for current repetition
            {
                step_func(rnd, xyz); // do a step
                r2_block[k] += r_sq(xyz) / L;
            }
        }

        for (int k = 0; k < N_steps; k++)
        {
            if (i == N_blocks - 1)
                step_list.push_back(k + 1);

            r_av[k + 1] += sqrt(r2_block[k]) / N_blocks; // add root of mean squared radius for k + 1 step to r_av
            r_av2[k + 1] += r2_block[k] / N_blocks;      // add mean squared radius for k + 1 step to r_av2
        }
    }

    std::vector<double> r_error = blocking::error(r_av, r_av2, N_blocks); // error on mean radius

    return {step_list, r_av, r_error};
}