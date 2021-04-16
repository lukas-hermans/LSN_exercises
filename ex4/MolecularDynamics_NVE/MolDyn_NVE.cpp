/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h> // srand, rand: to generate random number
#include <iostream> // cin, cout: Standard Input/Output Streams Library
#include <fstream>  // Stream class to both read and write from/to files.
#include <cmath>    // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main()
{
    Input(); //Inizialization
    int nconf = 1;
    for (int istep = 1; istep <= nstep; ++istep)
    {
        Move(); //Move particles with Verlet algorithm
        if (istep % iprint == 0)
        {
            cout << "Number of time-steps: " << istep << endl;
        }
        if (istep % 10 == 0)
        {
            Measure(); //Properties measurement
            //ConfXYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
            nconf += 1;
        }
        if (istep == nstep - 1)
        {
            ConfBeforeFinal(); // Write penultimate configuration to restart
        }
    }
    ConfFinal(); //Write final configuration to restart

    return 0;
}

//Prepare all stuff for the simulation
void Input(void)
{
    ifstream ReadInput, ReadConf;
    double ep, ek, pr, et, vir;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl
         << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl
         << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator

    ReadInput.open("input.dat"); //Read input

    ReadInput >> restart;
    cout << "Restart using final and penultimate configuration: " << restart << endl;

    ReadInput >> temp;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart / rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol, 1.0 / 3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl
         << endl;
    ReadInput.close();

    //Prepare array for measurements
    iv = 0;      //Potential energy
    ik = 1;      //Kinetic energy
    ie = 2;      //Total energy
    it = 3;      //Temperature
    n_props = 4; //Number of observables

    if (restart == 1) // read final and penultimate configuration and move one time to get velocity
    {
        // read final configuration
        cout << "Read final configuration from file config.final " << endl;
        ReadConf.open("config.final");
        for (int i = 0; i < npart; ++i)
        {
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();

        // read penultimate configuration
        cout << "Read penultimate configuration from file config.final-1 " << endl
             << endl;
        ReadConf.open("config.final-1");
        for (int i = 0; i < npart; ++i)
        {
            ReadConf >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadConf.close();

        Move();
    }
    if (restart == 0) // read final configuration and generate random velocities
    {
        //Read initial configuration
        cout << "Read initial configuration from file config.0 " << endl
             << endl;
        ReadConf.open("config.0");
        for (int i = 0; i < npart; ++i)
        {
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();

        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl
             << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i = 0; i < npart; ++i)
        {
            vx[i] = rand() / double(RAND_MAX) - 0.5;
            vy[i] = rand() / double(RAND_MAX) - 0.5;
            vz[i] = rand() / double(RAND_MAX) - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim = 0; idim < 3; ++idim)
            sumv[idim] /= (double)npart;
        double sumv2 = 0.0, fs;
        for (int i = 0; i < npart; ++i)
        {
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];
        }
    }

    // calculate squared velocity (per particle)
    double sumv2 = 0.0, fs;
    for (int i = 0; i < npart; ++i)
        sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
    sumv2 /= (double)npart;

    // rescale velocities to achieve desired temperature
    fs = sqrt(3 * temp / sumv2); // fs = velocity scale factor
    for (int i = 0; i < npart; ++i)
    {
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;

        xold[i] = Pbc(x[i] - vx[i] * delta);
        yold[i] = Pbc(y[i] - vy[i] * delta);
        zold[i] = Pbc(z[i] - vz[i] * delta);
    }

    return;
}

//Move particles with Verlet algorithm
void Move(void)
{
    double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

    for (int i = 0; i < npart; ++i)
    { //Force acting on particle i
        fx[i] = Force(i, 0);
        fy[i] = Force(i, 1);
        fz[i] = Force(i, 2);
    }

    for (int i = 0; i < npart; ++i)
    { //Verlet integration scheme

        xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
        ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
        znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

        vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
        vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
        vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];

        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
    return;
}

//Compute forces as -Grad_ip V(r)
double Force(int ip, int idir)
{
    double f = 0.0;
    double dvec[3], dr;

    for (int i = 0; i < npart; ++i)
    {
        if (i != ip)
        {
            dvec[0] = Pbc(x[ip] - x[i]); // distance ip-i in pbc
            dvec[1] = Pbc(y[ip] - y[i]);
            dvec[2] = Pbc(z[ip] - z[i]);

            dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
            dr = sqrt(dr);

            if (dr < rcut)
            {
                f += dvec[idir] * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8)); // -Grad_ip V(r)
            }
        }
    }

    return f;
}

//Properties measurement
void Measure()
{
    int bin;
    double v, t, vij;
    double dx, dy, dz, dr;
    ofstream Epot, Ekin, Etot, Temp;

    Epot.open("output_epot.dat", ios::app);
    Ekin.open("output_ekin.dat", ios::app);
    Temp.open("output_temp.dat", ios::app);
    Etot.open("output_etot.dat", ios::app);

    v = 0.0; //reset observables
    t = 0.0;

    //cycle over pairs of particles
    for (int i = 0; i < npart - 1; ++i)
    {
        for (int j = i + 1; j < npart; ++j)
        {

            dx = Pbc(xold[i] - xold[j]); // here I use old configurations [old = r(t)]
            dy = Pbc(yold[i] - yold[j]); // to be compatible with EKin which uses v(t)
            dz = Pbc(zold[i] - zold[j]); // => EPot should be computed with r(t)

            dr = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr);

            if (dr < rcut)
            {
                vij = 4.0 / pow(dr, 12) - 4.0 / pow(dr, 6);

                //Potential energy
                v += vij;
            }
        }
    }

    //Kinetic energy
    for (int i = 0; i < npart; ++i)
        t += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

    stima_pot = v / (double)npart;                //Potential energy per particle
    stima_kin = t / (double)npart;                //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t / (double)npart; //Temperature
    stima_etot = (t + v) / (double)npart;         //Total energy per particle

    Epot << stima_pot << endl;
    Ekin << stima_kin << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

//Write penultimate configuration
void ConfBeforeFinal(void)
{
    ofstream WriteConf;

    cout << "Print penultimate configuration to file config.final-1 " << endl;
    WriteConf.open("config.final-1");

    for (int i = 0; i < npart; ++i)
    {
        WriteConf << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
    }
    WriteConf.close();
    return;
}

//Write final configuration
void ConfFinal(void)
{
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl
         << endl;
    WriteConf.open("config.final");

    for (int i = 0; i < npart; ++i)
    {
        WriteConf << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
    }
    WriteConf.close();
    return;
}

void ConfXYZ(int nconf)
{ //Write configuration in .xyz format
    ofstream WriteXYZ;

    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i = 0; i < npart; ++i)
    {
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
    }
    WriteXYZ.close();
}

//Algorithm for periodic boundary conditions with side L=box
double Pbc(double r)
{
    return r - box * rint(r / box);
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
