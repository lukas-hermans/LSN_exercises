#ifndef __blocking__
#define __blocking__

#include <string>
#include <vector>
#include <iostream>
#include <functional>
#include "random.h"

namespace blocking
{
    double error(double av, double av_2, int n);
    std::vector<double> error(std::vector<double> av, std::vector<double> av_2, int n);
    std::vector<std::vector<double>> prog_mean(std::function<double(Random &)> func, Random &rnd, int N, int L);
    void write_data(std::vector<std::vector<double>> data, std::string path, std::string header);
    void write_data(std::vector<double> data, std::string path, std::string header); // overloading of function above
} // namespace blocking

#endif // __blocking__