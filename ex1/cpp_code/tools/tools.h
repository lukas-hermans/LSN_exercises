#ifndef __tools__
#define __tools__

#include <string>
#include <vector>

namespace tools
{
    double error(double av, double av_2, double n);
    std::vector<double> square_vec(std::vector<double> vec);
    std::vector<std::vector<double>> prog_mean(std::vector<double> vec);
    void count(std::vector<int> &counter, double draw, int N_intervals, double int_length);
    double chi2(std::vector<int> counter, int N_intervals, int n);
    std::vector<int> vec_arange(int start, int stop, int step = 1);
    void write_data(std::vector<int> row_num, std::vector<double> column1, std::vector<double> column2, std::string path = "data.txt", std::string header = "no header specified");
} // namespace tools

#endif // __tools__