#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../../tools/random.h"
#include "../../tools/salesman.h"

int main()
{
    int seed[4] = {1033, 1282, 549, 641}; // seed for Rannyu generator
    Random rnd = Random(seed);            // prepared instance for Rannyu generator

    int n_cities = 32;  // number of cities
    int pop_size = 500; // size of population
    int n_gens = 250;   // number of generations

    double p_c = 0.5; // crossover probability
    double p_m = 0.1; // mutation probability

    std::vector<point> circle_map = make_map(rnd, n_cities, "circle");
    write_map("../data/circle/map.txt", circle_map);

    std::vector<point> square_map = make_map(rnd, n_cities, "square");
    write_map("../data/square/map.txt", square_map);

    genetic_salesman circle_gen = genetic_salesman("../data/circle/", rnd, pop_size, n_gens, circle_map, p_c, p_m);
    circle_gen.gen_start_pop();
    circle_gen.evolute_pop(n_gens);

    genetic_salesman square_gen = genetic_salesman("../data/square/", rnd, pop_size, n_gens, square_map, p_c, p_m);
    square_gen.gen_start_pop();
    square_gen.evolute_pop(n_gens);

    rnd.SaveSeed();

    return 0;
}
