/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
    exercise(); // compute everything necessary for exercise 6

    // Input(); //Inizialization

    // for (int iblk = 1; iblk <= nblk; ++iblk) //Simulation
    // {
    //     Reset(iblk); //Reset block averages

    //     for (int istep = 1; istep <= nstep; ++istep)
    //     {
    //         Move(metro);
    //         Measure();
    //         Accumulate(); //Update block averages
    //     }

    //     Averages(iblk); //Print results for current block
    // }

    // ConfFinal(); //Write final configuration

    return 0;
}

// SOLUTION OF EXERCISE 6
//
// Does not use seed.in nor config.final but starts each run of a simulation with a random configuration.
// Seed can be specified inside function.
// Conducts simulation using Metropolis algorithm for h=0 and h=0.02.
// Simulations are done on a range of temp./betas.
// List of final values saved in data folder.
// Same for Gibbs algorithm.
//
// So: choose algorithm --> choose h --> choose temp. --> do simulation (with data blocking) [fuction does all of this automatically]
//
// Equilibration before every simulation is ensured.
void exercise()
{
    cout << "execute everything necessary for 6-th. exercise" << endl;

    exer = 1; // program should know that function exercise is executed

    nrep = 30;     // number of repetitions for each method (metro and gibbs), influences spacing between temp. steps between 0.5 and 2
    nblk = 20;     // total number of blocks
    nstep = 10000; // number of steps per block
    nequi = 10000; // number of steps for equilibration

    nspin = 50; // number of spins
    J = 1;      // coupling constant

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    rnd = Random(seed);                     // prepared instance for Rannyu generator

    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility

    n_props = 4; //Number of observables

    //initial configuration (T = inf)
    for (int i = 0; i < nspin; ++i)
    {
        if (rnd.Rannyu() >= 0.5)
            s[i] = 1;
        else
            s[i] = -1;
    }

    // create vectors for different simulations
    vector<double> beta_list; // list of betas (inverse temperatures)
    for (int i = 0; i < nrep; i++)
    {
        beta_list.push_back(1 / (0.5 + (2.0 - 0.5) / (nrep - 1) * i));
    }

    vector<int> metro_list = {1, 1, 0, 0};
    vector<double> h_list = {0, 0.02, 0, 0.02};
    vector<string> path_list = {"../data/metro_h=0.txt", "../data/metro_h=0.02.txt", "../data/gibbs_h=0.txt", "../data/gibbs_h=0.02.txt"};

    // open files to write equilibration measurements
    ofstream file_metro05("../data/metro_equi_05.txt");
    file_metro05 << "M, u(M), chi(M), m(M)" << endl;
    ofstream file_metro2("../data/metro_equi_2.txt");
    file_metro2 << "M, u(M), chi(M), m(M)" << endl;
    ofstream file_gibbs05("../data/gibbs_equi_05.txt");
    file_gibbs05 << "M, u(M), chi(M), m(M)" << endl;
    ofstream file_gibbs2("../data/gibbs_equi_2.txt");
    file_gibbs2 << "M, u(M), chi(M), m(M)" << endl;

    for (int sim_nr = 0; sim_nr < (int)metro_list.size(); sim_nr++) // loop over different simulations
    {
        ofstream remove_file(path_list[sim_nr]); // clear file

        metro = metro_list[sim_nr];
        h = h_list[sim_nr];

        cout << "metro = " << metro << ", h = " << h << endl;

        for (int irep = 1; irep <= nrep; ++irep) // loop over beta values/temperatures
        {
            beta = beta_list[irep - 1]; // take a irep-th. beta value (that is a certain temperature)

            // equilibration
            for (int iequi = 1; iequi <= nequi; iequi++)
            {
                Move(metro);
                Measure();

                // print equilibration process
                if (irep == 1 && h == 0 && iequi <= 10000)
                {
                    if (metro == 1)
                    {
                        file_metro05 << iequi << ", " << walker[iu] / m_spin << ", " << walker[ix] / m_spin << ", " << walker[im] / m_spin << endl;
                    }
                    if (metro == 0)
                    {
                        file_gibbs05 << iequi << ", " << walker[iu] / m_spin << ", " << walker[ix] / m_spin << ", " << walker[im] / m_spin << endl;
                    }
                }
                if (irep == nrep && h == 0 && iequi <= 10000)
                {
                    if (metro == 1)
                    {
                        file_metro2 << iequi << ", " << walker[iu] / m_spin << ", " << walker[ix] / m_spin << ", " << walker[im] / m_spin << endl;
                    }
                    if (metro == 0)
                    {
                        file_gibbs2 << iequi << ", " << walker[iu] / m_spin << ", " << walker[ix] / m_spin << ", " << walker[im] / m_spin << endl;
                    }
                }
            }

            for (int iblk = 1; iblk <= nblk; ++iblk) // loop over blocks
            {
                Reset(iblk); //Reset block averages

                for (int istep = 1; istep <= nstep; ++istep) // loop over steps in current block
                {
                    Move(metro);
                    Measure();
                    Accumulate(); //Update block averages
                }

                Averages(iblk, irep, path_list[sim_nr], "T, u(T), u_error(T), c(T), c_error(T), chi(T), chi_error(T), m(T), m_error(T)"); //Print results for current block
            }
        }
    }
}

void Input(void)
{
    exer = 0; // program should know that only one run of program and not function exersize is executed

    ifstream ReadInput;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl
         << endl;
    cout << "Nearest neighbour interaction      " << endl
         << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl
         << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;

    //Read input informations
    ReadInput.open("input.dat");

    ReadInput >> temp;
    beta = 1.0 / temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl
         << endl;

    ReadInput >> restart; // if=1 restart from file "config.final"
    if (restart == 1)
        cout << "The program loads start configuration of \"config.final\"" << endl;
    else
        cout << "The program starts from a random configuration" << endl;

    rnd = Random(restart); // prepared instance for Rannyu generator

    ReadInput >> metro; // if=1 Metropolis else Gibbs

    ReadInput >> nblk;

    ReadInput >> nstep;

    if (metro == 1)
        cout << "The program perform Metropolis moves" << endl;
    else
        cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl
         << endl;
    ReadInput.close();

    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility

    n_props = 4; //Number of observables

    //initial configuration
    if (restart == 1)
    {
        ifstream ReadConf("config.final");

        if (ReadConf.is_open() == 0)
        {
            throw std::invalid_argument("config.final does not exist");
        }

        for (int i = 0; i < nspin; ++i)
        {
            ReadConf >> s[i];
        }
        ReadConf.close();
    }
    else
    {
        for (int i = 0; i < nspin; ++i)
        {
            if (rnd.Rannyu() >= 0.5)
                s[i] = 1;
            else
                s[i] = -1;
        }
    }

    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu] / (double)nspin << endl;
}

void Move(int metro)
{
    int o, sm;
    double delta_E;
    double p_ratio;
    double alpha, r;

    for (int i = 0; i < nspin; ++i)
    {
        //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        o = (int)(rnd.Rannyu() * nspin);
        sm = s[o];

        if (metro == 1) //Metropolis
        {
            delta_E = Boltzmann(sm * (-1), o) - Boltzmann(sm, o); // energy difference between proposed and old configuration
            p_ratio = exp(-beta * delta_E);                       // quotient of prob. of proposed and old configuration

            alpha = min(1.0, p_ratio); // acceptance probability (simplified version!)
            r = rnd.Rannyu();          // random number between 0 and 1

            if (r <= alpha) // accept trial
            {
                s[o] *= -1; // do spin flip
                accepted++;
            }

            attempted++;
        }

        else //Gibbs sampling
        {
            delta_E = Boltzmann(1, o) - Boltzmann(-1, o); // energy difference between spin up and spin down configuration
            alpha = 1 / (1 + exp(beta * delta_E));        // conditional property to get spin +1

            r = rnd.Rannyu(); // random number between 0 and 1

            if (r <= alpha)
                s[o] = 1;
            else
                s[o] = -1;

            attempted++;
            accepted++;
        }
    }
}

double Boltzmann(int sm, int ip)
{
    double ene = -J * sm * (s[Pbc(ip - 1)] + s[Pbc(ip + 1)]) - h * sm;
    return ene;
}

void Measure()
{
    double u = 0.0, m = 0.0;

    //cycle over spins
    for (int i = 0; i < nspin; ++i)
    {
        u += -J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)]);
        m += s[i];
    }
    walker[iu] = u;
    walker[ic] = beta * beta * u * u;
    walker[ix] = beta * m * m;
    walker[im] = m;
}

void Reset(int iblk) //Reset block averages
{

    if (iblk == 1)
    {
        for (int i = 0; i < n_props; ++i)
        {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for (int i = 0; i < n_props; ++i)
    {
        blk_av[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}

void Accumulate(void) //Update block averages
{

    for (int i = 0; i < n_props; ++i)
    {
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}

void Averages(int iblk, int irep, string path, string header) //Print results for current block
{
    stima_u = blk_av[iu] / blk_norm / (double)nspin; //Energy
    glob_av[iu] += stima_u;
    glob_av2[iu] += stima_u * stima_u;
    err_u = Error(glob_av[iu], glob_av2[iu], iblk);

    stima_c = blk_av[ic] / blk_norm / (double)nspin - (double)nspin * beta * beta * stima_u * stima_u; //Heat Capacity
    glob_av[ic] += stima_c;
    glob_av2[ic] += stima_c * stima_c;
    err_c = Error(glob_av[ic], glob_av2[ic], iblk);

    stima_m = blk_av[im] / blk_norm / (double)nspin; //Magnetization
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m * stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);

    stima_x = blk_av[ix] / blk_norm / (double)nspin; //susceptibility
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x * stima_x;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);

    if (exer == 0)
    {
        cout << "Block number " << iblk << endl;
        cout << "Acceptance rate " << accepted / attempted << endl
             << endl;

        ofstream Ene, Heat, Mag, Chi;
        const int wd = 12;

        Ene.open("output.ene.0", ios::app);
        Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
        Ene.close();

        Heat.open("output.heat.0", ios::app);
        Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;
        Heat.close();

        Mag.open("output.mag.0", ios::app);
        Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
        Mag.close();

        Chi.open("output.chi.0", ios::app);
        Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;
        Chi.close();

        cout << "----------------------------" << endl
             << endl;
    }

    if ((exer == 1) && (iblk == nblk))
    {
        ofstream file(path, ios::app);

        if (irep == 1)
            file << header << endl;

        file << 1 / beta << ", "
             << glob_av[iu] / (double)iblk << ", " << err_u << ", "
             << glob_av[ic] / (double)iblk << ", " << err_c << ", "
             << glob_av[ix] / (double)iblk << ", " << err_x << ", "
             << glob_av[im] / (double)iblk << ", " << err_m
             << endl;
    }
}

void ConfFinal(void)
{
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl
         << endl;
    WriteConf.open("config.final");
    for (int i = 0; i < nspin; ++i)
    {
        WriteConf << s[i] << endl;
    }
    WriteConf.close();

    rnd.SaveSeed();
}

int Pbc(int i) //Algorithm for periodic boundary conditions
{
    if (i >= nspin)
        i = i - nspin;
    else if (i < 0)
        i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if (iblk == 1)
        return 0.0;
    else
        return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) / (double)(iblk - 1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
