#ifndef __salesman__
#define __salesman__

#include <string>

// Point structure that contains a x and y coordinate in the 2D space.
struct point
{
    double x;
    double y;
};

class genetic_salesman
{
public:
    genetic_salesman(string save_path, Random &rnd, int pop_size, int n_gens, std::vector<point> map, double p_c, double p_m, int rank = -1);
    ~genetic_salesman();

    void gen_start_pop();

    void evolute_pop(int steps);
    std::vector<int> gen_individual();
    void sim_anneal(std::vector<double> beta_list);
    void selection();
    void crossover();
    void mutation();
    void pair_permutation(std::vector<int> &individual);
    void check_individual(std::vector<int> individual);
    void shift_mutation(std::vector<int> &individual);
    void swap_mutation(std::vector<int> &individual);
    void inversion_mutation(std::vector<int> &individual);

    void compute_pop_fitness();
    double compute_individual_fitness(std::vector<int> individual);

    void write_data(int igen);

    string save_path; // path to store output files

    Random rnd; // random number generator

    int n_cities; // number of cities on map
    int pop_size; // size of population
    int n_gens;   // number of generations

    std::vector<point> map; // map of cities

    double p_c; // probability for crossover between two parents
    double p_m; // probability for each mutation for each new individual

    std::vector<int> city_list; // list of all cities in ascending order

    int only_mother = 0; // 1 if only mother should be computed for selection & mutation

    std::vector<std::vector<int>> pop;     // vector of current population
    std::vector<std::vector<int>> pop_new; // vector of new population
    std::vector<double> pop_fitness;       // vector of fitness for all individuals in current population
    std::vector<int> pop_order;            // contains order of index of all individuals of population (first index has highest fitness)

    std::vector<int> mother, father; // parent individuals

    // for Metropolis algorithm in simulated annealing
    double fitness_old;
    double prob_acc; // acceptance probability
    double r;        // random number between 0 and 1 to draw from acceptance probability

    int rank;        // rank of node (-1 indicates single-thread computation)
    int current_gen; // number of current generation
};

std::vector<point> make_map(Random &rnd, int n_cities, std::string type);
void write_map(std::string path, std::vector<point> city_map);
std::vector<int> vec_arange(int from, int to);
std::vector<double> vec_linspace(double from, double to, int steps);

#endif // __salesman__