/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <vector>

using namespace std;

#ifndef __Random__
#define __Random__

class Random
{

private:
    int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;

protected:
public:
    // constructors
    Random();
    Random(int seed[4]);
    Random(int restart);
    // destructor
    ~Random();
    // methods
    void SetRandom(int *, int, int);
    void SaveSeed();
    double Rannyu(void);
    double Rannyu(double min, double max);
    double Gauss(double mean, double sigma);
    double exponential_draw(double lambda);
    double lorentzian_draw(double gamma, double mu);
    double sample_cos();
    double importance_draw();
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
