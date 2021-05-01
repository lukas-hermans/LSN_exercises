#include <vector>
#include <cmath>
#include <iostream>
#include <functional>
#include "../../tools/random.h"
#include "../../tools/blocking.h"

double r_sq(std::vector<double> xyz);
void discrete_step(Random &rnd, std::vector<double> &xyz);
void continuous_step(Random &rnd, std::vector<double> &xyz);
std::vector<std::vector<double>> do_walk(Random &rnd, std::function<void(Random &, std::vector<double> &)> step_func, int N_repetitions, int N_steps);

int main()
{
    int N_repetitions = 10000; // number of repetitions of random walk
    int N_steps = 100;         // number of steps of random walker

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed);              // instance of random number generator

    std::vector<std::vector<double>> discrete_walk = do_walk(rnd, discrete_step, N_repetitions, N_steps);
    std::vector<std::vector<double>> continuous_walk = do_walk(rnd, continuous_step, N_repetitions, N_steps);

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
    double theta = acos(1 - 2 * rnd.Rannyu(0, 1));

    xyz[0] += cos(phi) * sin(theta);
    xyz[1] += sin(phi) * sin(theta);
    xyz[2] += cos(theta);
}

// Computes mean radius after 1, 2, 3, ..., N_steps steps and its error using N_blocks with L repetitions for each block.
// Returns a vector that contains three nested vectors:
//      1. Number of steps: 1, 2, ..., N_steps
//      2. Mean radius
//      3. Error
std::vector<std::vector<double>> do_walk(Random &rnd, std::function<void(Random &, std::vector<double> &)> step_func, int N_repetitions, int N_steps)
{
    std::vector<double> xyz;
    std::vector<double> r2_av(N_steps, 0.0);
    std::vector<double> r2_av2(N_steps, 0.0);

    for (int j = 0; j < N_repetitions; j++)
    {
        xyz = {0.0, 0.0, 0.0};            // starting point: origin
        for (int k = 0; k < N_steps; k++) // do random walk with N_steps
        {
            step_func(rnd, xyz);
            r2_av[k] += r_sq(xyz);
            r2_av2[k] += r_sq(xyz) * r_sq(xyz);
        }
    }

    std::vector<double> step_list;
    std::vector<double> r2_error;
    for (int i = 0; i < N_steps; i++)
    {
        r2_av[i] /= N_repetitions;
        r2_av2[i] /= N_repetitions;
        step_list.push_back(i + 1);
        r2_error.push_back(blocking::error(r2_av[i], r2_av2[i], N_repetitions));
    }

    // add starting point
    step_list.insert(step_list.begin(), 0.0);
    r2_av.insert(r2_av.begin(), 0.0);
    r2_error.insert(r2_error.begin(), 0.0);

    return {step_list, r2_av, r2_error};
}