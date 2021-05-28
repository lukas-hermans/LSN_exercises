#ifndef __salesman__
#define __salesman__

// Point structure that contains a x and y coordinate in the 2D space.
struct point
{
    double x;
    double y;
};

class genetic_salesman
{
public:
    genetic_salesman(string save_path, Random &rnd, int pop_size, int n_gens, std::vector<point> map, double p_c, double p_m);
    ~genetic_salesman();

    void gen_start_pop();

    void evolute_pop();
    void selection();
    void crossover();
    void mutation();
    void pair_permutation(std::vector<int> &individual);
    void check_individual(std::vector<int> individual);
    void shift_mutation(std::vector<int> &individual);
    void swap_mutation(std::vector<int> &individual);
    void inversion_mutation(std::vector<int> &individual);

    void compute_pop_fitness();

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

    std::vector<std::vector<int>> pop;     // vector of current population
    std::vector<std::vector<int>> pop_new; // vector of new population
    std::vector<double> pop_fitness;       // vector of fitness for all individuals in current population
    std::vector<int> pop_order;            // contains order of index of all individuals of population (first index has highest fitness)

    std::vector<int> mother, father; // parent individuals
};

std::vector<point> make_map(Random &rnd, int n_cities, std::string type);
void write_map(std::string path, std::vector<point> city_map);
std::vector<int> vec_arange(int from, int to);

#endif // __salesman__