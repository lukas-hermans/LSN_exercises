#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include "tools.h"

namespace tools
{
    double error(double av, double av_2, double n)
    {
        if (n == 1)
        {
            return 0;
        }
        else
        {
            return sqrt((av_2 - pow(av, 2)) / (n - 1));
        }
    }

    // Computes elementwise square of a vector.
    std::vector<double> square_vec(std::vector<double> vec)
    {
        for (int i = 0; i < (int)vec.size(); i++)
        {
            vec[i] = pow(vec[i], 2);
        }
        return vec;
    }

    // Computes progressive mean and its a error for a vector vec and returns it as a vector of two vectors.
    std::vector<std::vector<double>> prog_mean(std::vector<double> vec)
    {
        int N = vec.size();
        std::vector<double> vec_2 = square_vec(vec);

        std::vector<double> vec_prog(N, 0.0);
        std::vector<double> vec_2_prog(N, 0.0);
        std::vector<double> vec_error_prog;

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i + 1; j++)
            {
                vec_prog[i] += vec[j];
                vec_2_prog[i] += vec_2[j];
            }
            vec_prog[i] /= i + 1;
            vec_2_prog[i] /= i + 1;
            vec_error_prog.push_back(error(vec_prog[i], vec_2_prog[i], i + 1));
        }

        return std::vector<std::vector<double>>{vec_prog, vec_error_prog};
    }

    // Write data consisting of two columns into a file in a given path.
    void write_data(std::string path, std::string header, std::vector<double> column1, std::vector<double> column2)
    {
        int N = column1.size();
        std::ofstream file(path);
        file << header << std::endl;
        for (int i = 0; i < N; i++)
        {
            file << i + 1 << ", " << column1[i] << ", " << column2[i] << std::endl;
        }
        file.close();
    }

    // Searches for the subinterval in which the draw lies and increases the corresponding entry of the counter vector by 1.
    void count(std::vector<int> &counter, double draw, int N_intervals, double int_length)
    {
        for (int i = 0; i < N_intervals; i++)
        {
            if (draw >= int_length * i && draw < int_length * (i + 1))
            {
                counter[i] += 1;
                return;
            }
        }
    }

    // Computes the value of chi2 given an observed count of n random numbers over the subintervals.
    double chi2(std::vector<int> counter, int N_intervals, int n)
    {
        int n_exp = n / N_intervals; // expected number of draws in each subinterval
        double chi2 = 0;
        for (int i = 0; i < N_intervals; i++)
        {
            chi2 += pow(counter[i] - n_exp, 2) / (double)n_exp;
        }

        return chi2;
    }

} // namespace tools
