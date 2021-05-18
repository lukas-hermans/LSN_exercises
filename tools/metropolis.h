#ifndef __Metropolis__
#define __Metropolis__

#include <vector>
#include <string>
#include <functional>
#include "random.h"

class Metropolis
{
public:
    Random &rnd; // prepared instance of Random

    std::function<std::vector<double>(Random &, std::vector<double>)> T; // transition probability
    std::function<double(std::vector<double>)> p;                        // sample pdf

    bool is_equi = true;                // boolean to indicate if algorithm is in equilibrating phase
    int step_count = 0;                 // counts how many steps have already been completed
    std::vector<double> x0, x_current;  // start position and current position
    int save_step = 99999999;           // difference between steps that should be save (not too large, to avoid memory overflow)
    std::vector<std::vector<double>> x; // position history
    int acc = 0;                        // number of accepted trials

    std::function<double(std::vector<double>)> sample_func;                               // function that should be sampled
    std::vector<std::vector<double>> blocking_data = std::vector<std::vector<double>>(3); // stores data from use of Metropolis algorithm with blocking method (column 1: number of throws, 2: progressive mean, 3. error)

    Metropolis(Random &rnd);
    ~Metropolis();

    void set_x0(std::vector<double> x0);
    void do_step();
    void do_many_steps(int M_throws, int N_blocks, int N_equi);
    void write_x_history(std::string path);
    void write_hist(std::string path, int nbins, double xlow, double xup, int M_throws, int N_blocks);
    void clean_cache();
};

#endif // __Metropolis__