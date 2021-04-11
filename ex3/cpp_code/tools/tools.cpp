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

    // Returns a vector of an integer sequence specified by start, stop (is included) and step.
    std::vector<int> vec_arange(int start, int stop, int step)
    {
        std::vector<int> vec;
        for (int i = start; i < stop + 1; i += step)
            vec.push_back(i);
        return vec;
    }

    // Computes progressive mean and its error for a vector vec and returns it as a vector of two vectors.
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
    void write_data(std::vector<int> row_num, std::vector<double> column1, std::vector<double> column2, std::string path, std::string header)
    {
        int N = column1.size();
        std::ofstream file(path);
        file << header << std::endl;
        for (int i = 0; i < N; i++)
        {
            file << row_num[i] << ", " << column1[i] << ", " << column2[i] << std::endl;
        }
        file.close();
    }

} // namespace tools
