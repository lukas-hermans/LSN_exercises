/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props = 4;
double walker[m_props];
int n_props;
int iv, ik, it, ie;
double stima_pot, stima_kin, stima_etot, stima_temp;
double err_pot, err_kin, err_etot, err_temp;

// blocking
double blk_av[m_props], blk_norm;
double glob_av[m_props], glob_av2[m_props];

// averages
double acc, att;

//configuration
const int m_part = 108;
double x[m_part], y[m_part], z[m_part], xold[m_part], yold[m_part], zold[m_part];
double vx[m_part], vy[m_part], vz[m_part];

// thermodynamical state
int npart;
double energy, temp, vol, rho, box, rcut;

// simulation
int nblk, nstep, iprint, seed;
double delta;
bool restart;

//functions
void Input(void);
void Move(void);
void Accumulate(void);
void Averages(int);
void Reset(int);
void ConfBeforeFinal(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double sum, double sum2, int iblk);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
