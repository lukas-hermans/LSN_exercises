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
const int m_props = 1000;
double walker[m_props];
int n_props;
int iv, ik, it, ie, iw, igofr;
double vtail, ptail, bin_size, nbins, sd;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres, stima_g;
double err_pot, err_kin, err_etot, err_temp, err_press, err_gdir;

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
int nequi, equi_step, nblk, nstep, iprint, seed;
double delta;
bool restart;

//pigreco
const double pi = 3.1415927;

//functions
void Input(void);
void init_config(int);
void Move(void);
void Accumulate(void);
void Averages(int);
void Reset(int);
void ConfBeforeFinal(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(bool, int iequi = 0, int istep = 0);
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
