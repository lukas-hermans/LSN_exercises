#include <vector>
#include <string>
#include "../../tools/random.h"
#include "../../tools/salesman.h"

int main()
{
    int seed[4] = {1033, 1282, 549, 641}; // seed for Rannyu generator (should be same as for ex10 to generate identical city maps)
    Random rnd = Random(seed);            // prepared instance for Rannyu generator

    int n_cities = 32; // number of cities

    std::vector<point> circle_map = make_map(rnd, n_cities, "circle");
    write_map("../data/circle/map.txt", circle_map);

    std::vector<point> square_map = make_map(rnd, n_cities, "square");
    write_map("../data/square/map.txt", square_map);

    // parameters for simulated annealing
    std::vector<double> beta_list = vec_linspace(0, 1.0 / 0.01, 500); // list of beta values
    int pop_size = 1;
    int n_gens = 10000;
    double p_c = 0;
    double p_m = 0;

    // simulated annealing
    genetic_salesman circle_sim = genetic_salesman("../data/circle/", rnd, pop_size, n_gens, circle_map, p_c, p_m);
    circle_sim.gen_start_pop();
    circle_sim.sim_anneal(beta_list);
    genetic_salesman square_sim = genetic_salesman("../data/square/", rnd, pop_size, n_gens, square_map, p_c, p_m);
    square_sim.gen_start_pop();
    square_sim.sim_anneal(beta_list);

    return 0;
}