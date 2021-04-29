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
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
    Input(); //Inizialization

    for (int iblk = 1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk); //Reset block averages

        for (int istep = 1; istep <= nstep; ++istep)
        {
            Move(metro);
            Measure();
            Accumulate(); //Update block averages
        }

        Averages(iblk); //Print results for current block
    }

    ConfFinal(); //Write final configuration

    return 0;
}

void exersice()
{
    nspin = 50; // number of spins
    J = 1;      // coupling constant

    int seed[4] = {0000, 0000, 0000, 0001}; // seed for Rannyu generator
    rnd = Random(seed);                     // prepared instance for Rannyu generator

    nblk = 100;
    nstep = 10000;

    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility

    n_props = 4; //Number of observables
}

void Input(void)
{
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
    double energy_up, energy_down;

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
            alpha = 1 / (1 + exp(-beta * delta_E));       // conditional property to get spin +1

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
    int bin;
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

void Averages(int iblk) //Print results for current block
{

    ofstream Ene, Heat, Mag, Chi;
    const int wd = 12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted / attempted << endl
         << endl;

    Ene.open("output.ene.0", ios::app);
    stima_u = blk_av[iu] / blk_norm / (double)nspin; //Energy
    glob_av[iu] += stima_u;
    glob_av2[iu] += stima_u * stima_u;
    err_u = Error(glob_av[iu], glob_av2[iu], iblk);
    Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.heat.0", ios::app);
    stima_x = blk_av[ix] / blk_norm / (double)nspin - (double)nspin * beta * beta * stima_u * stima_u; //Heat Capacity
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x * stima_x;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);
    Heat << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;
    Heat.close();

    Mag.open("output.mag.0", ios::app);
    stima_m = blk_av[im] / blk_norm / (double)nspin; //Magnetization
    glob_av[im] += stima_m;
    glob_av2[im] += stima_m * stima_m;
    err_m = Error(glob_av[im], glob_av2[im], iblk);
    Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    Chi.open("output.chi.0", ios::app);
    stima_x = blk_av[ix] / blk_norm / (double)nspin; //susceptibility
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x * stima_x;
    err_x = Error(glob_av[ix], glob_av2[ix], iblk);
    Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    cout << "----------------------------" << endl
         << endl;
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
