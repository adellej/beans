#ifndef __EOS_H__
#define __EOS_H__

class Eos {
public:
  // density, temperature and composition
  double rho, T8;
  // initialization
  void init(int n);
  void tidy(void);
  double *A, *Z, *X;
  // mean molecular weights
  double Ye(void);
  double Yi(void);
  double YZ2(void);
  double Z2;
  // equation of state
  double pe(void);
  double pemod(void);
  double ptot(void);
  double Utot(void);
  double Chabrier_EF(void);
  double Fermi_Inv_1_2(double F);
  double FermiI(int k, double T8, double EF);
  // thermodynamics
  double f(void);
  double CP(void);
  double CV(void);
  double del_ad(void);
  double chi(double *x);
  double Gamma1(void);
  // reaction rates
  double triple_alpha(void);
  double O14AP(void);
  double O14AG(void);
  double O15AG(void);
  double N14PG(void);
  double N14AG(void);
  double C12PG(void);
  double C14AG(void);
  double Ne19PG(void);
  double F17PG(void);
  double F18PA(void);
  double C12C12(void);
  double fC12C12(void);
  double O16O16(void);
  double fO16O16(void);
  double C12O16(void);
  double fC12O16(void);
  double screen(double Z1, double Z2);
  double ogata(int i, int j);
  // opacity
  double eps_nu(void);
  double K_cond(double ef);
  double econd(void);
  double opac(void);
  double gff(double Z, double eta);
  double J(double x,double y);
  double iben(void);

  double lambda1, lambda2;
  double x(void);
  double gamma(void);
  double Uex(void);
  double eta(void);
  double lamei(void);

  double Q;

  double Fep(int flag);
  double Eep(double q);

  double debug;
  double kes, kff, kcond, kgaunt;

  double set_Ye, set_Yi;
  double RHO1,RHO2;

  double tmdrift(double dTdy);
  double Fermi(double n, double eta);
  void Fermi_derivs(double x, double ff[], double dfdx[]);
  static void Wrapper_Fermi_derivs(double x, double ff[], double dfdx[]);

  double find_rho(void);
  static double Wrapper_find_rho_eqn(double r);
  double find_rho_eqn(double r);
  double P;
  
 private:
  int ns;  // number of species
  double Fermi_n, Fermi_alpha;
};

#endif // __EOS_H__
