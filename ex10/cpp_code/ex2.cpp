#include "../../tools/random.h"
#include "../../tools/salesman.h"
#include <iostream>
#include "mpi.h"

int main(int argc, char *argv[])
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n_cities = 32;  // number of cities
    int pop_size = 500; // size of population
    int n_gens = 250;   // number of generations
    int N_migr = 50;    // number of generations between exchange between continents
    int n_steps;
    if (N_migr > 0)
        n_steps = n_gens / N_migr; // number of steps (i.e. generations) between migration processes
    else
        n_steps = n_gens;
    int max_exchange = 5; // maximum number of exchanging individuals in one migration process

    double p_c = 0.5; // crossover probability
    double p_m = 0.1; // mutation probability

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    Random rnd = Random(seed, rank);        // prepared instance for Rannyu generator (different primes for each node)

    // create point MPI_Datatype
    MPI_Datatype MPI_point;
    int count = 2;
    int blocklens[] = {1, 1};

    MPI_Aint indices[2];
    indices[0] = (MPI_Aint)offsetof(struct point, x);
    indices[1] = (MPI_Aint)offsetof(struct point, y);

    MPI_Datatype old_types[] = {MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_struct(count, blocklens, indices, old_types, &MPI_point);

    MPI_Type_commit(&MPI_point);

    // create square city map and synchronize it over all nodes (0 is root node)
    std::vector<point> square_map = make_map(rnd, n_cities, "square");
    if (rank == 0)
    {
        write_map("../data/parallel/map.txt", square_map);
    }
    MPI_Bcast(square_map.data(), n_cities, MPI_point, 0, MPI_COMM_WORLD);

    // initialize genetic population (different for each node)
    std::string path;
    if (N_migr == 0)
        path = "../data/non_parallel/";
    else
        path = "../data/parallel/";

    genetic_salesman square_gen = genetic_salesman(path, rnd, pop_size, n_gens, square_map, p_c, p_m, rank);
    square_gen.gen_start_pop();

    std::cout << "Initialized starting population for node " << rank << "/" << size << std::endl;

    std::vector<int> send_individual, recv_individual; // migrating individual
    int send_rank;                                     // rank that sends the migrating individuals
    int recv_rank;                                     // rank that receives a migrating individual
    int n_exchange;                                    // number of exchanging individuals in migration process
    MPI_Status stat;

    for (int i = 1; i <= N_migr; i++) // loop over number of migration processes
    {
        square_gen.evolute_pop(n_steps); // evolute population until next migration process

        if (rank == 0)
        {
            std::cout << i << ". migration process from " << N_migr << std::endl;
        }

        // determine sending and receiving rank for current migration process
        // and share information between nodes
        if (rank == 0)
        {
            send_rank = (int)rnd.Rannyu(0.0, (double)size);
            do
            {
                recv_rank = (int)rnd.Rannyu(0.0, (double)size);
            } while (recv_rank == send_rank);
            n_exchange = (int)rnd.Rannyu(1.0, (double)max_exchange);
        }
        MPI_Bcast(&send_rank, 1, MPI_INTEGER, rank, MPI_COMM_WORLD);
        MPI_Bcast(&recv_rank, 1, MPI_INTEGER, rank, MPI_COMM_WORLD);
        MPI_Bcast(&n_exchange, 1, MPI_INTEGER, rank, MPI_COMM_WORLD);

        for (int j = 0; j < n_exchange; j++) // loop over number of exchange individuals
        {
            // send migration individual from sending rank to receiving rank
            if (rank == send_rank)
            {
                send_individual = square_gen.pop[j];
                MPI_Send(&send_individual[0], n_cities - 1, MPI_INT, recv_rank, 1, MPI_COMM_WORLD);
            }
            // send migration individual from receiving to sending rank
            if (rank == recv_rank)
            {
                send_individual.resize(n_cities - 1);
                MPI_Recv(&send_individual[0], n_cities - 1, MPI_INT, send_rank, 1, MPI_COMM_WORLD, &stat);

                recv_individual = square_gen.pop[j];
                MPI_Send(&recv_individual[0], n_cities - 1, MPI_INT, send_rank, 1, MPI_COMM_WORLD);

                square_gen.pop[j] = send_individual;
            }
            if (rank == send_rank)
            {
                recv_individual.resize(n_cities - 1);
                MPI_Recv(&recv_individual[0], n_cities - 1, MPI_INT, recv_rank, 1, MPI_COMM_WORLD, &stat);

                square_gen.pop[j] = recv_individual;
            }
        }
        if (rank == send_rank || rank == recv_rank)
        {
            square_gen.compute_pop_fitness();
        }
    }

    if (N_migr == 0)
    {
        square_gen.evolute_pop(n_steps);
    }

    MPI_Finalize();

    return 0;
}