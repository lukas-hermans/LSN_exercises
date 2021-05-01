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
#include <cmath>
#include <cstdlib>
#include <vector>
#include "random.h"

using namespace std;

Random ::Random() {}

// Initialize an instance of Random that is prepared for the application of the Rannyu method.
Random::Random(int seed[4])
{
    int p1, p2;

    ifstream Primes("../../tools/Primes");
    Primes >> p1 >> p2;
    Primes.close();

    this->SetRandom(seed, p1, p2);
}

// Initialize an instance of Random that is prepared for the application of the Rannyu method.
// If restart=1, seed from "seed.out" is used, otherwise from "seed.in".
Random::Random(int restart)
{
    int p1, p2;
    int seed[4];

    ifstream Primes("../../tools/Primes");
    Primes >> p1 >> p2;
    Primes.close();

    ifstream Seed;

    if (restart == 1)
    {
        Seed.open("seed.out");
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    }
    else
    {
        Seed.open("seed.in");
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    }

    this->SetRandom(seed, p1, p2);
}

Random ::~Random() {}

void Random ::SaveSeed()
{
    ofstream WriteSeed;
    WriteSeed.open("seed.out");
    if (WriteSeed.is_open())
    {
        WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
        ;
    }
    else
        cerr << "PROBLEM: Unable to open random.out" << endl;
    WriteSeed.close();
    return;
}

double Random ::Gauss(double mean, double sigma)
{
    double s = Rannyu();
    double t = Rannyu();
    double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
    return mean + x * sigma;
}

double Random ::Rannyu(double min, double max)
{
    return min + (max - min) * Rannyu();
}

double Random ::Rannyu(void)
{
    const double twom12 = 0.000244140625;
    int i1, i2, i3, i4;
    double r;

    i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
    i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
    i3 = l3 * m4 + l4 * m3 + n3;
    i4 = l4 * m4 + n4;
    l4 = i4 % 4096;
    i3 = i3 + i4 / 4096;
    l3 = i3 % 4096;
    i2 = i2 + i3 / 4096;
    l2 = i2 % 4096;
    l1 = (i1 + i2 / 4096) % 4096;
    r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

    return r;
}

void Random ::SetRandom(int *s, int p1, int p2)
{
    m1 = 502;
    m2 = 1521;
    m3 = 4071;
    m4 = 2107;
    l1 = s[0];
    l2 = s[1];
    l3 = s[2];
    l4 = s[3];
    n1 = 0;
    n2 = 0;
    n3 = p1;
    n4 = p2;

    return;
}

// Returns an exponentially distributed random variable (computed via the inversion method from a uniformally distributed variable).
double Random::exponential_draw(double lambda)
{
    double y_uniform = this->Rannyu();
    return -1 / lambda * log(1 - y_uniform);
}

// Returns an Lorentzian distributed random variable (computed via the inversion method from a uniformally distributed variable).
double Random::lorentzian_draw(double gamma, double mu)
{
    double y_uniform = this->Rannyu();
    return mu + gamma * tan(M_PI * (y_uniform - 0.5));
}

// Samples cos(phi) with phi between 0 and pi/2 using the rejection method.
double Random::sample_cos()
{
    double x;
    double y;
    double r2;

    bool reject = true;
    while (reject == true) // loop until accepted
    {
        x = this->Rannyu(0, 1);
        y = this->Rannyu(0, 1);
        r2 = x * x + y * y;

        if (r2 < 1)
        {
            return x / sqrt(r2); // sample of cos(phi)
            reject = false;
        }
    }
}

// Returns a random variable that follows the pdf d(x) = 2 * (1 - x) (computed via the inversion method from a uniformally distributed variable).
double Random::importance_draw()
{
    double y_uniform = this->Rannyu();
    return 1 - sqrt(1 - y_uniform);
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
