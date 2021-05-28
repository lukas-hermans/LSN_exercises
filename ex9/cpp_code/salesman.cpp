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
    pop_fitness = std::vector<double>(pop_size, 0.0);
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
        pop.push_back(new_individual); // add new individual to population
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

        pop_new = {};

        // create new population using selection, crossover and mutation until pop_size is reached
        while ((int)pop_new.size() != pop_size)
        {
            this->selection(); // select new pair of parents (mother and father)
            this->crossover(); // perform crossover between parents
            this->mutation();  // perform random mutations on the parents

            // add parents (that are now children) to new population
            this->check_individual(mother);
            pop_new.push_back(mother);
            this->check_individual(father);
            pop_new.push_back(father);

            if ((int)pop_new.size() > pop_size)
            {
                throw std::runtime_error("ERROR: NEW POPULATION TOO LARGE! PLEASE SELECT EVEN POPULATION SIZE!");
            }
        }

        pop = pop_new;

        this->compute_pop_fitness();
        this->write_data(igen);
    }

    std::cout << "Evolution finished successfully" << std::endl;
}

// Selects mother and father individual in current population based on their fitness and stores them in mother/father vector.
void genetic_salesman::selection()
{
    int mother_index = (int)(pop_size * pow(rnd.Rannyu(), 2));
    int father_index;
    do
    {
        father_index = (int)(pop_size * pow(rnd.Rannyu(), 2));
    } while (father_index == mother_index);

    mother = pop[pop_order[mother_index]];
    father = pop[pop_order[father_index]];
}

// Performs crossover between mother and father.
void genetic_salesman::crossover()
{
    if (rnd.Rannyu() <= p_c)
    {
        int cut_index = rnd.Rannyu(1, n_cities - 2);

        // copy  of mother and father
        std::vector<int> mother_copy = mother;
        std::vector<int> father_copy = father;

        // last part of chromosomes to be crossed over
        std::vector<int> mother_cutted(mother.begin() + cut_index, mother.end());
        std::vector<int> father_cutted(father.begin() + cut_index, father.end());

        // index of cutted part (e.g. for the mother this is the index in the father chromosome)
        std::vector<int> mother_cutted_index;
        std::vector<int> father_cutted_index;
        for (int i = 0; i < (int)mother_cutted.size(); i++)
        {
            mother_cutted_index.push_back(std::distance(father.begin(), std::find(father.begin(), father.end(), mother_cutted[i])));
            father_cutted_index.push_back(std::distance(mother.begin(), std::find(mother.begin(), mother.end(), father_cutted[i])));
        }

        // sort indexes
        std::sort(mother_cutted_index.begin(), mother_cutted_index.end());
        std::sort(father_cutted_index.begin(), father_cutted_index.end());

        // change mother and father corresponding to performed crossover
        for (int i = 0; i < (int)mother_cutted.size(); i++)
        {
            mother[cut_index + i] = father_copy[mother_cutted_index[i]];
            father[cut_index + i] = mother_copy[father_cutted_index[i]];
        }
    }
}

// Performs several mutations each with probability p_m and changes mother and father correspondingly.
void genetic_salesman::mutation()
{
    // pair permutation
    if (rnd.Rannyu() <= p_m)
        this->pair_permutation(mother);
    if (rnd.Rannyu() <= p_m)
        this->pair_permutation(father);

    // shift mutation
    if (rnd.Rannyu() <= p_m)
        this->shift_mutation(mother);
    if (rnd.Rannyu() <= p_m)
        this->shift_mutation(father);

    // swap mutation
    if (rnd.Rannyu() <= p_m)
        this->swap_mutation(mother);
    if (rnd.Rannyu() <= p_m)
        this->swap_mutation(father);

    // inversion mutation
    if (rnd.Rannyu() <= p_m)
        this->inversion_mutation(mother);
    if (rnd.Rannyu() <= p_m)
        this->inversion_mutation(father);
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

// Computes the fitness of the current population in stores it in vector pop_fitness.
// The corresponding order of individuals are stored in pop_order.
void genetic_salesman::compute_pop_fitness()
{
    std::vector<int> individual;
    double fitness;
    int city_index_1, city_index_2;

    for (int i = 0; i < pop_size; i++) // loop over all individuals in current population
    {
        // compute fitness for current individual (which includes start city that is not part of individual vector) using L1 norm
        individual = pop[i];
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

        pop_fitness[i] = fitness;
    }

    // compute corresponding order vector
    pop_order = vec_arange(0, pop_size);
    std::sort(pop_order.begin(), pop_order.end(), [&](int i, int j)
              { return pop_fitness[i] < pop_fitness[j]; });
}

// Writes best fitness and fitness for best half into file.
// Also writes best path (i.e. a certain individual) to file.
// Output for l1 norm.
void genetic_salesman::write_data(int igen)
{
    // reset files at beginning of new run and write header
    if (igen == 0)
    {
        std::ofstream l1(save_path + "l1.txt"), l1_path(save_path + "l1_path.txt");
        l1 << "igen, l1_best, l1_mean" << std::endl;
        l1.close();
        l1_path.close();
    }

    std::ofstream l1, l1_path;
    l1.open(save_path + "l1.txt", std::ios_base::app);
    l1_path.open(save_path + "l1_path.txt", std::ios_base::app);

    // compute mean of l1 norm for best half of population
    double l1_mean = 0.0;
    for (int i = 0; i < pop_size / 2; i++)
    {
        l1_mean += pop_fitness[pop_order[i]];
    }
    l1_mean /= pop_size / 2;

    l1 << igen << ", " << pop_fitness[pop_order[0]] << ", " << l1_mean << std::endl;

    l1.close();

    // write best path to filef for l1
    std::vector<int> l1_best = pop[pop_order[0]];

    for (int i = 0; i < (int)l1_best.size(); i++)
    {
        l1_path << l1_best[i];

        if (i < (int)l1_best.size() - 1)
        {
            l1_path << ", ";
        }
        else
        {
            l1_path << std::endl;
        }
    }

    l1_path.close();
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