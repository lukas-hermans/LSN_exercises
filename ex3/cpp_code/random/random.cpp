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
#include <string>
#include "random.h"

using namespace std;

Random ::Random() {}

// Initialize an instance of Random that is prepared for the application of the Rannyu method.
Random::Random(int seed[4])
{
   int p1, p2;

   ifstream Primes("random/Primes");
   Primes >> p1 >> p2;
   Primes.close();

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
   double s = this->Rannyu();
   double t = this->Rannyu();
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

// Returns a vector of N prices calculated from blocks of size L using a direct calculation.
// Type argument specifies the option.
// It can be "put" or any other string, where any other string is considered as "call".
vector<double> Random::blocking_method_direct(string type, int N, int L, double S_0, double T, double K, double r, double sigma)
{
   vector<double> price; // to store option price

   double sum;
   double Z;   // to store standard normal number temporarily
   double S_T; // to store asset price at time T temporarily
   for (int i = 0; i < N; i++)
   {
      sum = 0;
      for (int j = 0; j < L; j++)
      {
         Z = (this->Gauss(0, 1));
         S_T = S_0 * exp((r - pow(sigma, 2) / 2) * T + sigma * Z * sqrt(T));
         (type == "put") ? (sum += exp(-r * T) * max(0.0, K - S_T)) : (sum += exp(-r * T) * max(0.0, S_T - K));
      }

      price.push_back(sum / (double)L); // append price of block i
   }

   return price;
}

// Returns a vector of N prices calculated from blocks of size L using a discretized calculation.
// Type argument specifies the option.
// It can be "put" or any other string, where any other string is considered as "call".
vector<double> Random::blocking_method_discretized(string type, int N, int L, int N_steps, double S_0, double T, double K, double r, double sigma)
{
   vector<double> price; // to store option price

   double sum;
   double Z;      // to store standard normal number temporarily
   double S_disc; // to store discretely calculated asset price (progressively)
   double S_T;    // to store asset price at time T temporarily
   for (int i = 0; i < N; i++)
   {
      sum = 0;
      for (int j = 0; j < L; j++)
      {
         S_disc = S_0;
         for (int step = 0; step < N_steps; step++) // discretized calculation of asset price at time T
         {
            Z = (this->Gauss(0, 1));
            S_disc = S_disc * exp((r - 0.5 * pow(sigma, 2)) * T / N_steps + sigma * Z * sqrt(T / N_steps));
         }
         S_T = S_disc;

         (type == "put") ? (sum += exp(-r * T) * max(0.0, K - S_T)) : (sum += exp(-r * T) * max(0.0, S_T - K));
      }

      price.push_back(sum / (double)L); // append price of block i
   }

   return price;
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
