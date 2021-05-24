#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <exception>
#include <numeric>
#include "../../tools/random.h"
#include "salesman.h"

genetic_salesman::genetic_salesman(string save_path, Random &rnd, int pop_size, int n_gens, std::vector<point> map, double p_c, double p_m) : save_path(save_path), rnd(rnd), pop_size(pop_size), n_gens(n_gens), map(map), p_c(p_c), p_m(p_m)
{
    n_cities = map.size();
    city_list = vec_arange(1, n_cities);
    pop_fitness_l1 = std::vector<double>(pop_size, 0.0);
    pop_fitness_l2 = std::vector<double>(pop_size, 0.0);
}

genetic_salesman::~genetic_salesman()
{
}

// Generates the start popoulation where every combination of cities (that occur exactly ones) has equal probability.
void genetic_salesman::gen_start_pop()
{
    std::vector<int> remaining_cities, new_individual;
    int add_city_index;

    for (int i = 0; i < pop_size; i++) // loop over population (i.e. over all individuals in population)
    {
        remaining_cities = city_list; // initially, the salesman has not visited any city
        new_individual = {};

        for (int j = 1; j < n_cities; j++) // loop over all cities (except first that we fix) to create new individual
        {
            add_city_index = (int)rnd.Rannyu(0, remaining_cities.size()); // random index of remaining_city vector (note that (int) alwys rounds down)

            new_individual.push_back(remaining_cities[add_city_index]);
            remaining_cities.erase(remaining_cities.begin() + add_city_index); // drop drawn city (salesman visits every city exactly once)
        }

        this->check_individual(new_individual);
        pop_l1.push_back(new_individual); // add new individual to population
        pop_l2.push_back(new_individual); // add new individual to population
    }

    this->compute_pop_fitness(); // compute fitness of current population
    this->write_data(0);         // write data

    std::cout << "Generated start population successfully" << std::endl;
}

void genetic_salesman::evolute_pop()
{
    std::cout << std::endl
              << "Evolution starts" << std::endl;

    for (int igen = 1; igen < n_gens + 1; igen++) // loop over all generations
    {
        if (igen % 10 == 0)
        {
            std::cout << "Generation: " << igen << "/" << n_gens << std::endl;
        }

        pop_new_l1 = {};
        pop_new_l2 = {};

        // create new population using selection, crossover and mutation until pop_size is reached
        while ((int)pop_new_l2.size() != pop_size)
        {
            this->selection(); // select new pair of parents (mother and father)
            this->crossover(); // perform crossover between parents
            this->mutation();  // perform random mutations on the parents

            // add parents (that are now children) to new population
            this->check_individual(mother_l1);
            pop_new_l1.push_back(mother_l1);
            this->check_individual(father_l1);
            pop_new_l1.push_back(father_l1);
            this->check_individual(mother_l2);
            pop_new_l2.push_back(mother_l2);
            this->check_individual(father_l2);
            pop_new_l2.push_back(father_l2);

            if ((int)pop_new_l1.size() > pop_size || (int)pop_new_l2.size() > pop_size)
            {
                throw std::runtime_error("ERROR: NEW POPULATION TOO LARGE! PLEASE SELECT EVEN POPULATION SIZE!");
            }
        }

        pop_l1 = pop_new_l1;
        pop_l2 = pop_new_l2;

        this->compute_pop_fitness();
        this->write_data(igen);
    }

    std::cout << "Evolution finished successfully" << std::endl;
}

// Selects mother and father individual in current population based on their fitness and stores them in mother/father vector.
void genetic_salesman::selection()
{
    int mother_index_l1 = (int)(pop_size * pow(1 - rnd.Rannyu(), 2));
    int father_index_l1;
    do
    {
        father_index_l1 = (int)(pop_size * pow(1 - rnd.Rannyu(), 2));
    } while (father_index_l1 == mother_index_l1);

    int mother_index_l2 = (int)(pop_size * pow(1 - rnd.Rannyu(), 2));
    int father_index_l2;
    do
    {
        father_index_l2 = (int)(pop_size * pow(1 - rnd.Rannyu(), 2));
    } while (father_index_l2 == mother_index_l2);

    mother_l1 = pop_l1[pop_order_l1[mother_index_l1]];
    father_l1 = pop_l1[pop_order_l1[father_index_l1]];
    mother_l2 = pop_l2[pop_order_l2[mother_index_l2]];
    father_l2 = pop_l2[pop_order_l2[father_index_l2]];
}

// Performs crossover between mother and father.
void genetic_salesman::crossover()
{
    if (rnd.Rannyu() <= p_c)
    {
        int cut_index = rnd.Rannyu(1, n_cities - 2);

        // copy  of mother and father
        std::vector<int> mother_l1_copy = mother_l1;
        std::vector<int> father_l1_copy = father_l1;

        // last part of chromosomes to be crossed over
        std::vector<int> mother_l1_cutted(mother_l1.begin() + cut_index, mother_l1.end());
        std::vector<int> father_l1_cutted(father_l1.begin() + cut_index, father_l1.end());

        // index of cutted part (e.g. for the mother this is the index in the father chromosome)
        std::vector<int> mother_l1_cutted_index;
        std::vector<int> father_l1_cutted_index;
        for (int i = 0; i < (int)mother_l1_cutted.size(); i++)
        {
            mother_l1_cutted_index.push_back(std::distance(father_l1.begin(), std::find(father_l1.begin(), father_l1.end(), mother_l1_cutted[i])));
            father_l1_cutted_index.push_back(std::distance(mother_l1.begin(), std::find(mother_l1.begin(), mother_l1.end(), father_l1_cutted[i])));
        }

        // sort indexes
        std::sort(mother_l1_cutted_index.begin(), mother_l1_cutted_index.end());
        std::sort(father_l1_cutted_index.begin(), father_l1_cutted_index.end());

        // change mother and father corresponding to performed crossover
        for (int i = 0; i < (int)mother_l1_cutted.size(); i++)
        {
            mother_l1[cut_index + i] = father_l1_copy[mother_l1_cutted_index[i]];
            father_l1[cut_index + i] = mother_l1_copy[father_l1_cutted_index[i]];
        }
    }
    if (rnd.Rannyu() <= p_c)
    {
        int cut_index = rnd.Rannyu(1, n_cities - 2);

        // copy  of mother and father
        std::vector<int> mother_l2_copy = mother_l2;
        std::vector<int> father_l2_copy = father_l2;

        // last part of chromosomes to be crossed over
        std::vector<int> mother_l2_cutted(mother_l2.begin() + cut_index, mother_l2.end());
        std::vector<int> father_l2_cutted(father_l2.begin() + cut_index, father_l2.end());

        // index of cutted part (e.g. for the mother this is the index in the father chromosome)
        std::vector<int> mother_l2_cutted_index;
        std::vector<int> father_l2_cutted_index;
        for (int i = 0; i < (int)mother_l2_cutted.size(); i++)
        {
            mother_l2_cutted_index.push_back(std::distance(father_l2.begin(), std::find(father_l2.begin(), father_l2.end(), mother_l2_cutted[i])));
            father_l2_cutted_index.push_back(std::distance(mother_l2.begin(), std::find(mother_l2.begin(), mother_l2.end(), father_l2_cutted[i])));
        }

        // sort indexes
        std::sort(mother_l2_cutted_index.begin(), mother_l2_cutted_index.end());
        std::sort(father_l2_cutted_index.begin(), father_l2_cutted_index.end());

        // change mother and father corresponding to performed crossover
        for (int i = 0; i < (int)mother_l2_cutted.size(); i++)
        {
            mother_l2[cut_index + i] = father_l2_copy[mother_l2_cutted_index[i]];
            father_l2[cut_index + i] = mother_l2_copy[father_l2_cutted_index[i]];
        }
    }
}

// Performs several mutations each with probability p_m and changes mother and father correspondingly.
void genetic_salesman::mutation()
{
    // pair permutation
    if (rnd.Rannyu() <= p_m)
        this->pair_permutation(mother_l1);
    if (rnd.Rannyu() <= p_m)
        this->pair_permutation(mother_l2);
    if (rnd.Rannyu() <= p_m)
        this->pair_permutation(father_l1);
    if (rnd.Rannyu() <= p_m)
        this->pair_permutation(father_l2);

    // shift mutation
    if (rnd.Rannyu() <= p_m)
        this->shift_mutation(mother_l1);
    if (rnd.Rannyu() <= p_m)
        this->shift_mutation(mother_l2);
    if (rnd.Rannyu() <= p_m)
        this->shift_mutation(father_l1);
    if (rnd.Rannyu() <= p_m)
        this->shift_mutation(father_l2);

    // swap mutation
    if (rnd.Rannyu() <= p_m)
        this->swap_mutation(mother_l1);
    if (rnd.Rannyu() <= p_m)
        this->swap_mutation(mother_l2);
    if (rnd.Rannyu() <= p_m)
        this->swap_mutation(father_l1);
    if (rnd.Rannyu() <= p_m)
        this->swap_mutation(father_l2);

    // inversion mutation
    if (rnd.Rannyu() <= p_m)
        this->inversion_mutation(mother_l1);
    if (rnd.Rannyu() <= p_m)
        this->inversion_mutation(mother_l2);
    if (rnd.Rannyu() <= p_m)
        this->inversion_mutation(father_l1);
    if (rnd.Rannyu() <= p_m)
        this->inversion_mutation(father_l2);
}

// Performs a random pair permutation on given individual.
void genetic_salesman::pair_permutation(std::vector<int> &individual)
{
    int index_1, index_2;

    index_1 = (int)rnd.Rannyu(0, n_cities - 1);
    do
    {
        index_2 = (int)rnd.Rannyu(0, n_cities - 1);
    } while (index_1 == index_2);

    std::swap(individual[index_1], individual[index_2]);
}

// Performs a random shift mutation on given individual.
void genetic_salesman::shift_mutation(std::vector<int> &individual)
{
    int index_1 = (int)rnd.Rannyu(0, n_cities - 2);
    int index_2 = (int)rnd.Rannyu(index_1 + 1, n_cities - 1);
    int m = (int)rnd.Rannyu(1, index_2 - index_1 + 1);

    std::vector<int> sub_individual(individual.begin() + index_1, individual.begin() + index_2);
    std::rotate(sub_individual.begin(), sub_individual.begin() + m, sub_individual.end());

    std::copy(sub_individual.begin(), sub_individual.end(), individual.begin() + index_1);
}

// Performs a random swap mutation on given individual.
void genetic_salesman::swap_mutation(std::vector<int> &individual)
{
    int index_1 = (int)rnd.Rannyu(0, n_cities - 2);
    int index_2 = (int)rnd.Rannyu(index_1 + 1, n_cities - 1);

    int m = (int)rnd.Rannyu(0, (double)min({n_cities - 2 - index_2, index_2 - index_1}));

    for (int i = 0; i < m; i++)
    {
        std::swap(individual[index_1 + i], individual[index_2 + i]);
    }
}

// Performs a random inversion mutation on given individual.
void genetic_salesman::inversion_mutation(std::vector<int> &individual)
{
    int index_1 = (int)rnd.Rannyu(0, n_cities - 2);
    int index_2 = (int)rnd.Rannyu(index_1 + 1, n_cities - 1);

    std::reverse(individual.begin() + index_1, individual.begin() + index_2);
}

// Raises a warning if the individual does not fullfill requirement that every city is visited exavtly once.
void genetic_salesman::check_individual(std::vector<int> individual)
{
    std::sort(individual.begin(), individual.end()); // sort elements of vector in ascending order (only locally in this function)

    if (individual != city_list)
    {
        std::cout << "WARNING: POPULATION CONTAINS INCORRECT INDIVIDUAL:" << std::endl;
        std::cout << "Sorted individual: ";
        for (int i = 0; i < (int)individual.size(); i++)
        {
            std::cout << individual[i] << " ";
        }
        std::cout << std::endl;
    }
}

// Computes the fitness of the current population in stores it in vector pop_fitness_l1 and pop_fitness_l2.
// The corresponding order of individuals are stored in pop_order_l1 and pop_order_l2.
void genetic_salesman::compute_pop_fitness()
{
    std::vector<int> individual;
    double fitness;
    int city_index_1, city_index_2;

    for (int i = 0; i < pop_size; i++) // loop over all individuals in current population
    {
        // compute fitness for current individual (which includes start city that is not part of individual vector) using L1 norm
        individual = pop_l1[i];
        fitness = 0;
        city_index_1 = 0;
        city_index_2 = individual[0];
        fitness += sqrt(pow(map[city_index_1].x - map[city_index_2].x, 2) + pow(map[city_index_1].y - map[city_index_2].y, 2));
        for (int j = 0; j < (int)individual.size() - 1; j++) // loop over city sequence in individual (except last one)
        {
            city_index_1 = individual[j];
            city_index_2 = individual[j + 1];
            fitness += sqrt(pow(map[city_index_1].x - map[city_index_2].x, 2) + pow(map[city_index_1].y - map[city_index_2].y, 2));
        }
        city_index_1 = individual.back();
        city_index_2 = 0;
        fitness += sqrt(pow(map[city_index_1].x - map[city_index_2].x, 2) + pow(map[city_index_1].y - map[city_index_2].y, 2));

        pop_fitness_l1[i] = fitness;

        // compute fitness for current individual (which includes start city that is not part of individual vector) using L2 norm
        individual = pop_l2[i];
        fitness = 0;
        city_index_1 = 0;
        city_index_2 = individual[0];
        fitness += pow(map[city_index_1].x - map[city_index_2].x, 2) + pow(map[city_index_1].y - map[city_index_2].y, 2);
        for (int j = 0; j < (int)individual.size() - 1; j++) // loop over city sequence in individual (except last one)
        {
            city_index_1 = individual[j];
            city_index_2 = individual[j + 1];
            fitness += pow(map[city_index_1].x - map[city_index_2].x, 2) + pow(map[city_index_1].y - map[city_index_2].y, 2);
        }
        city_index_1 = individual.back();
        city_index_2 = 0;
        fitness += pow(map[city_index_1].x - map[city_index_2].x, 2) + pow(map[city_index_1].y - map[city_index_2].y, 2);

        pop_fitness_l2[i] = fitness;
    }

    // compute corresponding order vectors
    pop_order_l1 = vec_arange(0, pop_size);
    std::sort(pop_order_l1.begin(), pop_order_l1.end(), [&](int i, int j)
              { return pop_fitness_l1[i] < pop_fitness_l1[j]; });

    pop_order_l2 = vec_arange(0, pop_size);
    std::sort(pop_order_l2.begin(), pop_order_l2.end(), [&](int i, int j)
              { return pop_fitness_l2[i] < pop_fitness_l2[j]; });
}

// Writes best fitness and fitness for best half into file.
// Also writes best path (i.e. a certain individual) to file.
// Output for l1 and l2 norm.
void genetic_salesman::write_data(int igen)
{
    // reset files at beginning of new run and write header
    if (igen == 0)
    {
        std::ofstream l1(save_path + "l1.txt"), l1_path(save_path + "l1_path.txt"), l2(save_path + "l2.txt"), l2_path(save_path + "l2_path.txt");
        l1 << "igen, l1_best, l1_mean" << std::endl;
        l1.close();
        l1_path.close();
        l2 << "igen, l2_best, l2_mean" << std::endl;
        l2.close();
        l2_path.close();
    }

    std::ofstream l1, l1_path, l2, l2_path;
    l1.open(save_path + "l1.txt", std::ios_base::app);
    l1_path.open(save_path + "l1_path.txt", std::ios_base::app);
    l2.open(save_path + "l2.txt", std::ios_base::app);
    l2_path.open(save_path + "l2_path.txt", std::ios_base::app);

    // compute mean of l1 and l2 norm for best half of population
    double l1_mean = 0.0;
    double l2_mean = 0.0;
    for (int i = 0; i < pop_size / 2; i++)
    {
        l1_mean += pop_fitness_l1[pop_order_l1[i]];
        l2_mean += pop_fitness_l2[pop_order_l2[i]];
    }
    l1_mean /= pop_size / 2;
    l2_mean /= pop_size / 2;

    l1 << igen << ", " << pop_fitness_l1[pop_order_l1[0]] << ", " << l1_mean << std::endl;
    l2 << igen << ", " << pop_fitness_l2[pop_order_l2[0]] << ", " << l2_mean << std::endl;

    l1.close();
    l2.close();

    // write best path to filef for l1 and l2
    std::vector<int> l1_best = pop_l1[pop_order_l1[0]];
    std::vector<int> l2_best = pop_l2[pop_order_l2[0]];

    for (int i = 0; i < (int)l1_best.size(); i++)
    {
        l1_path << l1_best[i];
        l2_path << l2_best[i];

        if (i < (int)l1_best.size() - 1)
        {
            l1_path << ", ";
            l2_path << ", ";
        }
        else
        {
            l1_path << std::endl;
            l2_path << std::endl;
        }
    }

    l1_path.close();
    l2_path.close();
}

// Creates a vector of points (with x and y coordinate) of n_cities.
// Type can be "circle" or "square".
std::vector<point> make_map(Random &rnd, int n_cities, std::string type)
{
    std::vector<point> city_vector;
    point new_city;

    if (type == "circle")
    {
        double phi;
        for (int i = 0; i < n_cities; i++)
        {
            phi = rnd.Rannyu(0, 2 * M_PI);
            new_city.x = cos(phi);
            new_city.y = sin(phi);

            city_vector.push_back(new_city);
        }
    }

    if (type == "square")
    {
        for (int i = 0; i < n_cities; i++)
        {
            new_city.x = rnd.Rannyu(-1, 1);
            new_city.y = rnd.Rannyu(-1, 1);

            city_vector.push_back(new_city);
        }
    }

    return city_vector;
}

// Writes city map to a specified path.
void write_map(std::string path, std::vector<point> city_map)
{
    std::ofstream file(path);
    file << "x, y" << std::endl;

    for (int i = 0; i < (int)city_map.size(); i++)
    {
        file << city_map[i].x << ", " << city_map[i].y << std::endl;
    }
}

std::vector<int> vec_arange(int from, int to)
{
    std::vector<int> vec;
    for (int i = from; i < to; i++)
    {
        vec.push_back(i);
    }

    return vec;
}