#include <functional>
#include <vector>
#include <cmath>
#include "../../tools/random.h"
#include "../../tools/blocking.h"

double do_needle_throw(Random &rnd, double L_needle, double d);

int main()
{
    double L_needle = 0.7; //length of needle
    double d = 1;          // distance between straight lines on horizontal plane

    int N = 100;    // number of blocks
    int M = 100000; // total number of throws
    int L = M / N;  // throws per block

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed);              // prepared instance for Rannyu generator

    std::vector<std::vector<double>> pi_prog = blocking::prog_mean(std::bind(do_needle_throw, std::placeholders::_1, L_needle, d), rnd, N, L);
    blocking::write_data(pi_prog, "../data/pi_buffon.txt", "M, N_hit / N_thr(M), pi_error(M)");

    return 0;
}

// Simulation of Buffon's experiment to estimate the number pi.
// L_needle is the length of the needle, while d describes the distance between the two straigth lines.
double do_needle_throw(Random &rnd, double L_needle, double d)
{
    double cos_phi = rnd.sample_cos();
    double delta_x = rnd.Rannyu(0, d / 2);

    return (delta_x <= L_needle / 2 * cos_phi) ? 1 : 0; // return 1 if hit, else 0
}