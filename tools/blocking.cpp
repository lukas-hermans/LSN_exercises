#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>
#include "random.h"
#include "blocking.h"

namespace blocking
{
    double error(double av, double av_2, int n)
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

    // Computes progressive mean and its error using blocking method with N blocks of L throws.
    // Returns a vector that contains three nested vectors:
    //      1. Number of throws M
    //      2. Progressive mean
    //      3. Error
    std::vector<std::vector<double>> prog_mean(std::function<double(Random &)> func, Random &rnd, int N, int L)
    {
        std::vector<double> av, av2;

        // calculate av and av2 of given function for each block
        double sum = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < L; j++) // loop over throws in i-th block
            {
                sum += func(rnd);
            }
            av.push_back(sum / L);
            av2.push_back(pow(av[i], 2));

            sum = 0;
        }

        std::vector<std::vector<double>> data(3);

        // calculate progressive mean and its error from av and av2
        double sum_prog = 0, sum2_prog = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i + 1; j++)
            {
                sum_prog += av[j];
                sum2_prog += av2[j];
            }
            data[0].push_back((i + 1) * L);
            data[1].push_back(sum_prog / (i + 1));
            data[2].push_back(error(data[1][i], sum2_prog / (i + 1), i + 1));

            sum_prog = sum2_prog = 0;
        }

        return data;
    }

    // Write data from several vectors in a vector into a file in a given path.
    void write_data(std::vector<std::vector<double>> data, std::string path, std::string header)
    {
        std::ofstream file(path);
        file << header << std::endl;

        for (int i = 0; i < (int)data[0].size(); i++) // loop over all elements in each nested vector
        {
            for (int j = 0; j < (int)data.size(); j++) // loop over data columns
            {
                file << data[j][i];
                if (j < (int)data.size() - 1) // donnot put comma for final data columns
                {
                    file << ", ";
                }
            }
            file << std::endl;
        }

        file.close();
    }

    // Overloading of write_data function above if data is just one vector and not a nested one.
    void write_data(std::vector<double> data, std::string path, std::string header)
    {
        std::ofstream file(path);
        file << header << std::endl;

        for (int i = 0; i < (int)data.size(); i++) // loop over all elements in vector
        {
            file << data[i] << std::endl;
        }

        file.close();
    }

} // namespace blocking
