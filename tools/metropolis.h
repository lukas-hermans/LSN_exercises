#ifndef __Metropolis__
#define __Metropolis__

#include <vector>
#include <functional>
#include "random.h"

class Metropolis
{
public:
    Random &rnd; // prepared instance of Random

    std::function<std::vector<double>(Random &, std::vector<double>)> T; // transition probability
    std::function<double(std::vector<double>)> p;                        // sample pdf
    std::vector<std::vector<double>> x;                                  // position history

    std::function<double(std::vector<double>)> sample_func; // function that should be sampled

    std::vector<double> samples_block;      // vector of samples for all blocks
    std::vector<double> prog_samples_error; // error on progressive mean of samples (using blockking method)

    Metropolis(Random &rnd);
    ~Metropolis();

    void set_x0(std::vector<double> x0);
    void do_step();
    void do_many_steps(int M_throws, int N_blocks);
};

#endif // __Metropolis__