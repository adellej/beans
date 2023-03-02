#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>

extern "C" {
#include "nrutil.h"
#include "root.h"
}

#define me 510.999
#define F14 0.352
#define F15 0.648

#include "odeint.h"
#include "eos.h"

void* pt2Object;

// ------------------------ initialise ----------------------------------

void Eos::tidy(void)
{
    free_dvector(this->A,1,this->ns);
    free_dvector(this->Z,1,this->ns);
    free_dvector(this->X,1,this->ns);
}

void Eos::init(int n)
{
  this->ns=n;

  this->A=dvector(1,this->ns);
  this->Z=dvector(1,this->ns);
  this->X=dvector(1,this->ns);
  this->Z2=0.0;
  this->set_Ye=0.0;
  this->set_Yi=0.0;
}

// ------------------------ mean molecular weights ------------------------

double Eos::Yi(void)
  // inverse mean molecular weight per ion
{
  //  return this->X+(this->Y/4.0)+(F14*this->Z/14.0)+(F15*this->Z/15.0);
  int i;
  double sum=0.0;

  if (set_Yi > 0.0) sum=set_Yi;
  else {
    for (i=1; i<=this->ns; i++)
      sum+=this->X[i]/this->A[i];
  }
  return sum;
}

double Eos::Ye(void)
  // inverse mean molecular weight per electron
{
  //  return this->X+(this->Y/2.0)+(F14*this->Z*8.0/14.0)+(F15*this->Z*8.0/15.0);
  int i;
  double sum=0.0;

  if (set_Ye > 0.0) sum=set_Ye;
  else {
    for (i=1; i<=this->ns; i++)
      sum+=this->X[i]*this->Z[i]/this->A[i];
  }
  return sum;
}

double Eos::YZ2(void)
  // sum of (X/A)*Z^2
{
  //  return this->X+this->Y+(F14*this->Z*64.0/14.0)+(F15*this->Z*64.0/15.0);
  int i;
  double sum=0.0;
  if (this->Z2 == 0.0) {
    for (i=1; i<=this->ns; i++)
      sum+=this->X[i]*this->Z[i]*this->Z[i]/this->A[i];
  } else {
    sum=this->Z2;  
  }
  return sum;
}

// ------------------------ equation of state ------------------------------

double Eos::pe(void)
  // Calculates the (cgs) electron gas pressure as given by the
  // semi-analytic formula of Paczynski (1983) ApJ 267 315.
{
  double rY, pednr, pedr, pend, ped;
  rY=this->rho*Ye();
  pednr=9.91e-2*pow(rY,5.0/3.0);
  pedr=1.231e1*pow(rY,4.0/3.0);
  ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
  pend=8.254e-7*1e8*this->T8*rY;
  return(1e14*sqrt(pow(ped,2)+pow(pend,2)));
}

double Eos::ptot(void)
  // Calculates the total pressure, that is the sum of electron +
  // ion + radiation pressure. Coulomb corrections included if
  // species 1 has Z>1 i.e. not hydrogen, this is a trick to only
  // apply Coulomb corrections in the ocean
{
  double f;
  if (this->Z[1]>1.0) f=1.0+this->Uex()/3.0; else f=1.0;
  return 8.254e15*this->rho*this->T8*Yi()*f + pe() 
    + 2.521967e17*pow(this->T8,4);
}


double Eos::pemod(void)
  // Modified version of the electron pressure formula
  // to use for heat capacity
  // (Paczynski 1983 ApJ 267 315)
{
  double rY, pednr, pedr, pend, ped;
  rY=this->rho*this->Ye();
  // divide pednr and pedr by appropriate factors
  pednr=9.91e-2*pow(rY,5.0/3.0)/1.32;
  pedr=1.231e1*pow(rY,4.0/3.0)/0.822;
  ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
  pend=8.254e-7*1e8*this->T8*rY;
  return(1e14*sqrt(pow(ped,2)+pow(pend,2)));
}




double Eos::Utot(void)
  // internal energy density
{
  double f, r, tau;
  tau=this->T8/59.4;
  r=this->FermiI(2,this->T8,this->Chabrier_EF())/
    this->FermiI(1,this->T8,this->Chabrier_EF());
  if (this->Z[1]>1.0) f=1.5+this->Uex(); else f=1.5;
  return 8.254e15*this->rho*this->T8*Yi()*f + 
    1.5*this->pe()*(1.0+tau*r)/(1.0+0.5*tau*r) +
    3.0*2.521967e17*pow(this->T8,4);
}


double Eos::FermiI(int k, double T8, double EF)
  // fitting formula for the generalized Fermi integrals
  // from Chabrier and Potekhin 1998
  // nu=k+1/2, with k=0,1, or 2
{
  double c[3][5]={{0.37045057, 0.41258437, 9.777982e-2, 5.3734153e-3, 3.8746281e-5},
		  {0.39603109, 0.69468795, 0.22322760, 1.5262934e-2, 1.3081939e-4},
		  {0.76934619, 1.7891437, 0.70754974, 5.6755672e-2, 5.5571480e-4}};
  double e[3][5]={{0.43139881, 1.7597537, 4.1044654, 7.7467038, 13.457678},
		  {0.81763176, 2.4723339, 5.1160061, 9.0441465, 15.049882},
		  {1.2558461, 3.2070406, 6.1239082, 10.316126, 16.597079}};
  double x[5]={7.265351e-2, 0.2694608, 0.533122, 0.7868801, 0.9569313};
  double z[5]={0.26356032, 1.4134031, 3.5964258, 7.0858100, 12.640801};
  double h[5]={3.818735e-2, 0.1256732, 0.1986308, 0.1976334, 0.1065420};
  double v[5]={0.29505869, 0.32064856, 7.3915570e-2, 3.6087389e-3, 2.3369894e-5};
  int i;
  double I=0.0, F, R, chi, tau;

  tau=T8/59.4;
  chi=EF/(8.625*T8);

  R=sqrt(chi*(1.0+0.5*chi*tau));

  if (chi*tau > 0.01) {
    F=(chi+1.0/tau)*0.5*R-log(1+tau*chi+sqrt(2.0*tau)*R)/pow(2*tau,1.5);
    if (k>0) F=(2*pow(R,3.0)/3.0-F)/tau;
    if (k>1) F=(2*chi*pow(R,3.0)-5.0*F)/(4.0*tau);
  } else {
    F=pow(chi,1.0*k+1.5)/(1.0*k+1.5);
  }

  if (chi <= 0.6) {
    for (i=0; i<5; i++) {
      I+=c[k][i]*sqrt(1+0.5*e[k][i]*tau)/(exp(-e[k][i])+exp(-chi));
    }
  }
  if (chi > 0.6 && chi < 14.0) {
    for (i=0; i<5; i++) {
      I+=h[i]*pow(x[i],1.0*k)*pow(chi, 1.0*k+1.5)*sqrt(1+0.5*chi*x[i]*tau)/
	  (1+exp(chi*x[i]-chi));
      I+=v[i]*pow(z[i]+chi, 1.0*k+0.5)*sqrt(1+0.5*(z[i]+chi)*tau);
    }
  }
  if (chi >= 14.0) {
    I=F+(M_PI*M_PI/6.0)*pow(chi, 1.0*k)*(1.0*k+0.5+0.5*(k+1)*chi*tau)/R;
  }

  return I;
}


double Eos::Fermi_Inv_1_2(double F)
{
  double AN=0.5, RN, DEN, INV, FF;
  int i;
  int M1=2, K1=2, M2=2, K2=2;
  double A1[4]={0.0, 4.4593646e1, 1.1288764e1, 1.0};
  double B1[4]={0.0, 3.9519346e1, -5.7517464, 2.6594291e-1};
  double A2[4]={0.0, 3.4873722e1, -2.6922515e1, 1.0};
  double B2[4]={0.0, 2.6612832e1, -2.0452930e1, 1.1808945e1};
  
  if (F < 4.0) {
    RN=F+A1[M1];
    for (i=M1-1; i>=1; i--) RN=RN*F+A1[i];
    DEN=B1[K1+1];
    for (i=K1; i>=1; i--) DEN=DEN*F+B1[i];
    INV=log(F*RN/DEN);
  } else {
    FF=1.0/pow(F,1.0/(1.0+AN));
    RN=FF+A2[M2];
    for (i=M2-1; i>=1; i--) RN=RN*FF+A2[i];
    DEN=B2[K2+1];
    for (i=K2; i>=1; i--) DEN=DEN*FF+B2[i];
    INV=RN/(DEN*FF);
  }
  
  return INV;
}


double Eos::Chabrier_EF(void)
  // Calculates the Fermi energy in keV including the
  // rest mass using the fit of Chabrier and Potekhin,
  // (1998) Phys Rev E, 58, 4941
  // It uses Antia (1993) to evaluate the inverse Fermi integral
  //
  // rY is (rho/mu_e) in g/cm^3 ;  T is the temperature in K
{
  double EFnr, kT, x, tau, theta;
  double q1,q2,q3, et,etu, corr,F, mc2, rY, T;

  T=this->T8*1e8;
  rY=this->rho*Ye();

  // Electron rest mass in keV
  mc2=510.999;

  // Find kT, x=p_F/m_e c, tau=kT/me and theta=T/T_F
  kT=8.617347*T*1e-8;
  x=1.007e-2*pow(rY,1.0/3.0);
  tau=kT/mc2;
  theta=tau/(sqrt(1+x*x)-1);

  // Calculate the non-relativistic guess
  F=2.0*pow(theta,-1.5)/3.0;
  EFnr=kT*Fermi_Inv_1_2(F);

  // These functions are defined in CP eq. 24
  et=exp(theta); etu=1/et;
  q1=1.5/(et-1);
  q2=12+8/pow(theta,1.5);
  q3=1.366-(etu+1.612*et)/(6.192*pow(theta,0.0944)*etu+
			   5.535*pow(theta,0.698)*et);

  // This is the correction to the non-relativistic EF
  corr=(1+q1*sqrt(tau)+q2*q3*tau)/(1+q2*tau);
  corr*=tau/(1+(tau/(2*theta)));
  corr=1.5*log(1+corr);

  // return E_F including the rest mass

  return mc2+EFnr-kT*corr;
}

// ----------------------- thermodynamics --------------------------------

double Eos::chi(double *x)
{
  double x1, x2, p1, p2;
  x1=*x; p1=ptot();
  x2=*x=1.001*x1; p2=ptot();
  *x=x1;
  return (log(p2)-log(p1))/(log(x2)-log(x1));
}


double Eos::CP(void)
  // Calculates specific heat at constant pressure
{
  double cv, chiT, chirho;

  chirho=chi(&this->rho);
  chiT=chi(&this->T8);
  cv=CV();

  return cv+chiT*chiT*ptot()/(this->rho*this->T8*1e8*chirho);
}

double Eos::Gamma1(void)
{
  double chirho, chiT, gam1, cv;

  chirho=chi(&this->rho);
  chiT=chi(&this->T8);
  cv=CV();
  
  gam1=chirho+chiT*chiT*ptot()/(cv*this->rho*1e8*this->T8);

  return gam1;
}  


double Eos::del_ad(void)
  // calculates dlnT/dlnp at constant entropy
{
  double chirho, chiT, gam1, cv;

  chirho=chi(&this->rho);
  chiT=chi(&this->T8);
  cv=CV();
  
  gam1=chirho+chiT*chiT*ptot()/(cv*this->rho*1e8*this->T8);

  return ptot()*chiT/(cv*gam1*this->rho*1e8*this->T8);
}

double Eos::CV(void)
  // Calculates the specific heat at constant volume (density)
{
  double gg, cvion, cve, cvrad, alpha;

  // first the IONS
  gg=this->gamma();
  
  {  // alpha comes from Chabrier's fit
    double a1,a2,a3,b1,b2,b3,b4;
    a1=-0.9070; a2=0.62954; a3=-0.5*sqrt(3.0)-a1/sqrt(a2);
    b1=4.56e-3; b2=211.6; b3=-1.0e-4; b4=4.62e-3;
    alpha=0.5*pow(gg,1.5)*(a3*(gg-1.0)/pow(gg+1.0,2.0)-a1*a2/pow(gg+a2,1.5))
      +pow(gg,2.0)*(b3*(pow(gg,2.0)-b4)/pow(pow(gg,2.0)+b4,2.0)-
		    b1*b2/pow(gg+b2,2.0));
  }
  cvion=8.3144e7*(1.5+alpha)*this->Yi();
  
  // RADIATION
  cvrad=3.0256e10*pow(this->T8,3)/this->rho;
  
  // ELECTRONS
  { // modified version of Paczynksi's fit for cve
    double dT,temp,p1,p2;
    temp=1.001*this->T8; dT=temp-this->T8;
    p1=this->pemod(); this->T8+=dT; p2=this->pemod(); this->T8-=dT;
    cve=(1/((this->f()-1)*this->rho))*1e-8*(p2-p1)/dT;
  }
  
  // total
  return cve+cvion+cvrad;
}

/*
double Eos::CV(void)
  // Calculates the specific heat at constant volume (density) by 
  // performing a numerical differentiation. The electron specific heat
  // is calculated using the fitting formula given by Paczynski (1983).
{
  double cv, T1, T2, p1, p2, U1, U2;
  T1=this->T8; p1=pe(); U1=this->Uex();
  this->T8=T2=1.001*T1; p2=pe(); U2=this->Uex();
  this->T8=T1;
  cv=0.0;
  // ions
  cv+=8.3144e7*1.5*this->Yi();
  if (this->Z[1] > 1.0) 
    cv+=8.3144e7*this->Yi()*(U2*T2-U1*T1)/(T2-T1);
  // electrons
  cv+=(1/((f()-1)*this->rho))*1e-8*(p2-p1)/(T2-T1);
  // radiation
  cv+=3.0256e10*pow(T1,3)/this->rho;
  return cv;
}
*/

double Eos::f(void)
  // Calculates f = dln ped/dln rho, using the fitting formula given
  // by Paczynski (1983).
{
  double rY, pednr, pedr, ped;
  rY=this->rho*Ye();
  pednr=9.91e-2*pow(rY,5.0/3.0);
  pedr=1.231e1*pow(rY,4.0/3.0);
  ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
  return (5.0*pow(ped/pednr,2) + 4.0*pow(ped/pedr,2))/3.0;
}



// --------------------------- opacity ------------------------------------


double Eos::eps_nu(void)
  // Calculates neutrino emissivity (erg/g/s)
  // by plasma process from Schinder et al. 1987
{
  double a0, a1, a2, b1, b2, b3, c;
  double xi, xi2, xi3, la, la2, la3, g, K;
  double Q1, Q2;

  // variables
  la=this->T8/59.302; la2=la*la; la3=la2*la;
  xi=pow(this->rho*this->Ye()*1e-9, 1.0/3.0)/la;
  xi2=xi*xi; xi3=xi2*xi;
  
  // 1. plasma

  // these coefficients valid for 10^8<T<10^11 K
  a0=2.146e-7; a1=7.814e-8; a2=1.653e-8;
  b1=2.581e-2; b2=1.734e-2; b3=6.990e-4;
  c=0.56457;

  K=pow(this->rho*this->Ye(),3.0);

  // formula from Schinder et al.
  Q1=K*exp(-c*xi)*(a0+a1*xi+a2*xi2)/(xi3+(b1/la)+(b2/la2)+(b3/la3));

  // 2. pair

  // coefficients valid for 10^8 < T < 10^11 K
  a0=5.026e19; a1=1.745e20; a2=1.568e21;
  if (this->T8 < 100.0) {  // 10^8<T<10^10 K
    b1=9.383e-1; b2=-4.141e-1; b3=5.829e-2;
    c=5.5924;
  } else { // 10^10 < T < 10^11 K
    b1=1.2383; b2=-8.141e-1; b3=0.0;
    c=4.9924;
  }

  g=1.0-13.04*la2+133.5*la2*la2+1534*la2*la2*la2+918.6*la2*la2*la2*la2;
  K=g*exp(-2.0/la);

  // formula from Schinder et al.
  Q2=K*exp(-c*xi)*(a0+a1*xi+a2*xi2)/(xi3+(b1/la)+(b2/la2)+(b3/la3));
  
  // return the summed emissivity (divide by rho to get per gram)
  return (Q1+Q2)/this->rho;
}


double Eos::opac(void)
  // Calculates the opacity
{
  double kappa, ef, eta, gaunt;
  int i;

  // This is the fitting formula for the electron scattering
  // opacity from Paczynski
  this->kes=(0.4*Ye())/((1+2.7e11*this->rho*pow(1e8*this->T8,-2.0))*
		  (1+pow(this->T8/4.5,0.86)));

  // Fermi energy
  ef=Chabrier_EF();
  if (ef == 0) {
    eta=Fermi_Inv_1_2(1.105e-4*this->rho*Ye()/pow(this->T8,1.5));
    ef=eta*8.617*this->T8;
  } else {
    eta=(ef-me)/(8.617*this->T8);
  }

  // Free-Free opacity
  this->kff=7.53e-6*this->rho*Ye()/pow(this->T8, 3.5);
  gaunt=0.0;
  
  // this is the general formula
  for (i=1; i<=this->ns; i++)
    gaunt+=this->Z[i]*this->Z[i]*this->X[i]*gff(this->Z[i],eta)/this->A[i];

  // fix for hot CNO
  //gaunt=this->X*gff(1,eta)+this->Y*gff(2,eta)+
  //64.0*((F15/15.0)+(F14/14.0))*this->Z*gff(8,eta);

  // this line is a fix to compare to hendrik/greg stuff (24 Oct 2001)
  //gaunt+=this->Z2*gff(this->Yi()/this->Ye(), eta);

  this->kff*=gaunt;

  // Conduction
  this->kcond=3.024e20*pow(this->T8,3)/(K_cond(ef)*this->rho);
 
  // Add up opacities
  kappa=1.0/((1.0/this->kcond)+(1.0/(this->kff+this->kes)));
  //  printf("%lg %lg %lg\n", this->kes, this->kff, this->kcond);
  
  return kappa;
}


double Eos::gff(double Z1, double eta)
  // Calculates the free-free Gaunt factor for element with
  // charge Z1 using a fitting
  // formula described in the Schatz et al steady state paper
  // Agrees to 10% with Itoh et al. 1991
{
  double gaunt, x, rY, T8_32, gam, ef;

  rY=this->rho*Ye();
  T8_32=pow(this->T8, 1.5);

  if (eta < 100.0) x=log(1.0+exp(eta)); 
  else x=eta; // make sure it doesn't freak out for extremely large eta

  // normalisation and degeneracy piece
  gaunt=1.16*8.02e3*x*T8_32/rY;

  x=pow(1+x,2.0/3.0);
  gam=sqrt(1.58e-3/this->T8)*Z1;

  // Elwert factor
  gaunt*=(1.0-exp(-2*M_PI*gam/sqrt(x+10.0)))/(1.0-exp(-2*M_PI*gam/sqrt(x)));

  // relativistic piece
  gaunt*=1.0+pow(this->T8/7.7, 1.5);

  // send it back
  return gaunt;
}


double Eos::Uex(void)
  // Coulomb correction Uex/kT
{
  double u,g,g2,g3,g14;
  g=this->gamma(); g2=g*g; g3=g2*g; g14=pow(g,0.25);

  if (g < 173.0) {
    u=-0.89813*g+0.98686*g14-0.91095+0.25098/g14;
  } else {
    u=-0.89593*g+1.5+9.65/g+840/g2+1.101e5/g3;
  }

  return u;
}

double Eos::Fep(int flag)
  // "Coulomb log" for electron-phonon scattering 
  // (Baiko & Yakovlev 1995,1996)
  // if flag=0 electrical conductivity; flag>0 thermal conductivity
{
  double R0, R1, R2, G0, G2, t, u2, u1, s, F, K0,K1;
  double alpha, alpha0, a0, a2, x, beta, AA, ZZ;
  double P0,g, K2,P2, c1, c2;

  // constants -- use values for bcc crystal
  a0=0.0174; a2=0.0118; u2=13.0; u1=2.8;

  AA=this->A[1]; ZZ=this->Z[1];
  x=this->x(); beta=x/sqrt(1+x*x);
  t=0.804*this->T8*(0.5*AA/ZZ)/sqrt(1e-9*this->rho);
  s=pow(4*ZZ,-2.0/3.0)+2.323e-3/beta;
  alpha0=1.683*sqrt(x/(AA*ZZ));
  alpha=alpha0*(0.5*u1*exp(-9.1*t)+t*u2);
  //alpha=1e-6;  //  small alpha is the Yakovlev & Urpin limit
  
  G0=u2*t/sqrt(t*t+a0); 
  R0=(exp(-alpha*s)-exp(-alpha))/alpha;
  R1=2*(exp(-alpha*s)*(1+alpha*s)-exp(-alpha)*(1+alpha))/
    (alpha*alpha);
  K0=2*R0-beta*beta*R1;
  
  // correction for finite nuclear size
  if (this->rho < 4e11) g=0.16*pow(this->rho*1e-12,1.0/3.0);
  else g=0.25*pow(this->rho*1e-12*this->Ye(),1.0/3.0);
  //  g=0.0;  switch off finite size effects
  P0=4.787-0.0346*ZZ;
  R2=(exp(-alpha*s)*(alpha*alpha*s*s+2*alpha*s+2)-exp(-alpha)*(alpha*alpha+2*alpha+2))/(alpha*alpha*alpha);
  c1=pow(1.0+pow(18.0*ZZ*M_PI,2.0/3.0)*g*g*(0.5*R1-beta*beta*R2)/(2.5*K0*P0),-P0);
  
  F=G0*K0*c1;
  
  if (flag > 0) { // thermal conductivity so add an extra piece
    P2=2.729-0.0204*ZZ;
    R2=this->Eep(alpha*s)-this->Eep(alpha);
    G2=t/(M_PI*M_PI*pow(t*t+a2,1.5));
    K2=0.5*R2-0.5*beta*beta*R0;
    // correction for finite nuclear size
    c2=pow(1.0+pow(18.0*M_PI*ZZ,2.0/3.0)*g*g*0.5*K0/(10.0*K2*P2),-P2);
    F+=G2*(3*K2-0.5*K0)*c2;
  }
  return F;
}

double Eos::Eep(double q)
  // used by Fep() to calculated thermal conductivity piece
  // Baiko & Yakovlev 1995
{
  double q2,q3,q4,qu;
  q2=q*q; q3=q2*q; q4=q3*q; qu=1.0/q;
  return exp(-q4/(q3+0.1397))*(log(1+qu)-0.5772/(1+2.2757*q2));
}


double Eos::lamei(void)
{
  double x,x1,x2,lam;
  x1=this->x();
  x2=0.22*sqrt(this->T8);
  if (x1>x2) x=x1; else x=x2;
  lam=127*x*sqrt((3.0/this->gamma())+1.5)/
    pow(this->rho*this->Yi(),1.0/3.0);
  return log(lam)-0.5*x*x/(1+x*x);
}  


double Eos::tmdrift(double dTdy)
  // Calculates the thermomagnetic drift velocity
  // dTdy is the local temperature gradient
{
  double u;
  
  // if species 1 is hydrogen, then assume we're in the H/He layer
  if (this->A[1]==1.0) {

    double tau, eta;
    eta=this->eta();
    tau=3.43e-13*pow(this->T8,1.5)/(this->lamei()*this->rho*this->YZ2());
    u=1.5165e11*tau*this->rho*dTdy;
    u*=11.0*this->Fermi(2.0,eta)*this->Fermi(4.5,eta)-
      12.0*this->Fermi(3.0,eta)*this->Fermi(3.5,eta);
    u/=6*pow(this->Fermi(2.0,eta),2.0);

  } else { // ashes
    
    double x, lam, beta;

    x=this->x(); beta=x/sqrt(1+x*x);  
    lam=log(pow(2*M_PI*this->Ye()/(3*this->Yi()),1.0/3.0)*
	    sqrt(1.5+3.0/this->gamma()));
    lam-=0.5*beta*beta;
    
    u=0.5*1.101e-22*this->T8*dTdy*this->econd()/this->Ye();
    u*=3.0-2.0*beta*beta-
      (1.0-beta*beta+beta*beta*beta*beta)/lam;
    u*=sqrt(1+x*x)/(x*x);
    
  }

  return u;
}



double Eos::econd(void)
  // calculates the electrical conductivity
{
  double x1, x2, sig, x, lambda, nu, theta, beta;

  if (this->gamma() < 173.0 || this->Q == 900.0) { // if Q=900 treat as liquid
    
    // This is the method from the WD paper, where I interpolate using x
    // choose appropriate value for x
    x1=this->x();
    x2=0.26*sqrt(this->T8);
    x=sqrt(x1*x1+x2*x2);
    if (x1>x2) x=x1; else x=x2;
    sig=8.48e21*this->Ye()*pow(x,3.0)/(this->YZ2()*(1+x*x));
    sig/=this->lamei();

    // here, write sig directly in terms of Fermi integrals
    //    sig=9.47e23*pow(this->T8,3.0)*this->Fermi(2.0,this->eta())/this->rho;
    //sig/=this->YZ2()*this->lamei();
    
  } else { // solid --- NB assumes A=2Z and single species
    double TU,ka,sm1;

    TU=2.2e8*sqrt(1e-12*this->rho)*this->Ye()*pow(this->Z[1]/60.0,1.0/3.0);

    x=this->x(); beta=x/sqrt(1+x*x);
    nu=9.55e16*this->T8*this->Fep(0)/beta; // phonons

    // add exponential suppression when the Umklapp scatterings freeze out
    //if (this->T8  < 1e-8*TU) 
    nu*=exp(-1e-8*TU/this->T8);

    //nu=9.55e16*this->T8*13.0/beta;
    /* old phonons from Urpin & Yakovlev
    theta=0.56*sqrt(1e-9*this->rho)/this->T8;
    nu=1.24e18*this->T8*(2-beta*beta)/(beta*sqrt(1+pow(theta/3.5,2.0)));
    */

    // Coulomb log from Itoh & Kohyama 1996
    ka=1.92*pow(this->Ye()/this->Yi(),1.0/3.0);
    sm1=0.5*log(1.0+0.4*ka*ka);
    lambda=sm1*(1.0+2.5*beta*beta/(ka*ka))-0.5*beta*beta;

    nu+=1.75e16*this->Q*lambda*sqrt(1+x*x)/this->Z[1]; // impurities

    //sig=1.49e22*x*x*beta*1e16/nu;
    sig=1.52e25*1e17*pow(this->rho*1e-12*Ye(),2.0/3.0)/nu;
  }  
  return sig;
}


double Eos::Fermi(double n, double eta)
{
  double F;
  Ode_Int FERMI;
  
  pt2Object=(void*) this;

  FERMI.init(1);

  //if (eta>30.0) return pow(eta,n+1.0)/(n+1.0);
  this->Fermi_n=n; this->Fermi_alpha=-eta;
  FERMI.set_bc(1,0);
  FERMI.go(0.0, 200.0, 1.0, 1e-8, &Eos::Wrapper_Fermi_derivs);
  F=FERMI.get_y(1,FERMI.kount);

  FERMI.tidy();
  return F;
}

void Eos::Wrapper_Fermi_derivs(double x, double ff[], double dfdx[])
{
  Eos* mySelf = (Eos*) pt2Object;
  // call member
  mySelf->Fermi_derivs(x, ff, dfdx);
}

void Eos::Fermi_derivs(double x, double ff[], double dfdx[])
{
  dfdx[1]=pow(x,Fermi_n)/(1+exp(x+Fermi_alpha));
}


double Eos::find_rho(void)
{
  double old, found, guess, rad, guess1;
  pt2Object=(void*) this;
  old=this->rho;

  //  found=zbrent(Wrapper_find_rho_eqn,1.0,1e12,1e-6);

  // first guess the density
  rad=2.521967e17*pow(this->T8,4);
  //rad=0.0;
  guess=1.49e5*pow(this->P*1e-22,0.75)/this->Ye();  // NR deg electrons
  if (guess < 1e4) guess=(this->P-rad)*1.66e-24/
		     ((this->Ye()+this->Yi())*1.38e-8*this->T8);
  guess1=guess;
  while ((found = zbrent(Wrapper_find_rho_eqn,1e-1*guess,10.0*guess,1e-8))
	 <= 0.12*guess) {
    guess/=9.0;
    //    printf("new guess=%lg\n", guess);
    //    printf("*"); fflush(stdout);
  }

  //rad=2.521967e17*pow(this->T8,4);
  //  printf("%lg %lg %lg %lg %lg\n", guess, found, rad, this->P, guess1);
  
  this->rho=old;
  return found;
}


double Eos::Wrapper_find_rho_eqn(double r)
{
  Eos* mySelf = (Eos*) pt2Object;
  // call member
  return mySelf->find_rho_eqn(r);
}

double Eos::find_rho_eqn(double r)
{
  this->rho=r;
  return this->ptot()-this->P;
}




double Eos::x(void)
{
  double x; 
  x=pow(this->Chabrier_EF()/511.0,2)-1.0; if (x<0.0) x=0.0;   
  x=sqrt(x);
  return x;
}

double Eos::eta(void)
{
  return (this->Chabrier_EF()-511.0)/(8.625*this->T8);
}

double Eos::gamma(void)
{ 
  return 0.11*(this->YZ2()/this->Yi())*pow(this->rho*1e-5*Yi(),1.0/3.0)/
    this->T8;
}

double Eos::K_cond(double ef)
  // Calculates the conductivity due to electron-ion and
  // electron-electron collisions
  // ef is the Fermi energy in keV
{
  double x, lam, f_ei, y, y3, f_ee, rY, x2, K;
  double gam, f_c, theta, beta, corr;

  rY=this->rho*Ye();
  x=this->x(); if (x==0) return 1e-10;
  x2=sqrt(1+x*x); beta=x/x2;

  gam=this->gamma();

  lam=log(pow(2*M_PI*Ye()/(3*Yi()),1.0/3.0)*sqrt(1.5+3.0/gam));
  lam-=0.5*beta*beta;
  this->lambda2=lam;

  if (gam < 173.0 || this->Q == 900.0) { // if Q=900 treat as liquid  

    // The electron-ion collision frequency is calculated as given by
    // Yakovlev & Urpin
    f_ei=1.755e16*this->lambda2*x2*YZ2()/Ye();
    
    // electron-electron collisions
    y=5.771e-3*sqrt(rY/x2)/this->T8;
    f_ee=5.11e15*this->T8*this->T8*pow(x,1.5)*J(x,y)/pow(1+x*x,1.25);
    
    // The collision frequencies just add
    f_c=f_ee+f_ei;

    //printf("%lg %lg\n", f_ee, f_ei);
 
  } else { // solid --- NB assumes A=2Z and single species

    /* old phonons from Yakovlev & Urpin
    theta=0.56*sqrt(1e-9*this->rho)/this->T8;
    lam=(2-beta*beta)/(beta*sqrt(1+pow(theta/3.5,2.0)));
    lam+=pow(theta/5.1,2.0)*(3*this->lambda2-1+0.5*beta*beta)/
      (beta*pow(1+pow(theta/4.2,2.0),1.5));
    f_c=1.24e18*this->T8*lam;
    */
    f_c=9.55e16*this->T8*this->Fep(1)/beta; // phonons

    lam=0.77; // this for <Z>=30
    f_c+=1.75e16*this->Q*lam*x2/this->Z[1];  // impurities
  }
  
  // the conductivity is then as given by Yakovlev & Urpin
  K = 4.116e27*this->T8*rY/(x2*f_c);

  // correction due to thermoelectric field
  //corr=9.34e-4*pow(this->T8,2.0)*(1+x*x)/pow(x,4.0);
  //corr*=pow(6.0-2.0*beta*beta-(1.0-beta*beta+beta*beta*beta*beta)/lam,2.0);
  //K=K*(1.0-corr); if (K<0.0) return 1e-10;

  return K; 
  
}



double Eos::J(double x,double y)
{
  // from Potekhin, Chabrier, & Yakovlev 1997
  double x2=x*x;
  double b2=x2/(1.+x2);
  double y3=y*y*y;
  double y4=y3*y;
  return (1.+0.4*(3.+1./x2)/x2)*(y3*pow(1.+0.07414*y,-3.)*
                                 log((2.810-0.810*b2+y)/y)/3.+
                                 pow(M_PI,5)*y4*pow(13.91+y,-4.)/6);
}


double Eos::iben(void)
  // This is the opacity fit used by Iben 1975
  // (what Taam uses) for free free
{
  double AA,B,C,Ap,mu,T617;
  double lkappa;

  T617=pow(this->T8*100.0,1.7);

  //  mu=G.avZ2A-1;
  //mu=1.0*ZZ*ZZ/AA;
  mu=this->X[1]+this->X[2]+64.0*((F15/15.0)+(F14/14.0))*this->X[3];
  mu-=1;

  AA=(1.25+0.488*sqrt(mu)+0.092*mu)/0.67;

  B=3.86+0.252*sqrt(mu)+0.018*mu;

  C=pow(2.019e-4*this->rho/T617, 2.425);

  Ap=1.0+C*(1.0+C/24.55);

  lkappa=0.67*(log10(this->rho)+AA-B*log10(100.0*this->T8))+log10(Ap);

  return pow(10.0,lkappa);
}



// ----------------------- reaction rates --------------------------------



double Eos::O14AP(void)
  // Energy production rate per gram for 14O(alpha,p)17F reaction
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9;

  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0);
  T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0);
  T932=pow(T9,1.5);
 
  Q=1.191;  //  Q value in MeV

  // f is the rate in (cm^3/mol)/s
  f = 1.68e+13/T923*exp(-39.388/T913-pow(T9/0.717,2.0))*
    (1.0+0.011*T913+13.117*T923+0.971*T9+85.295*T943+16.061*T953)+
    3.31e+04/T932*exp(-11.733/T9)+1.79e+07/T932*exp(-22.609/T9)+9.00e+03
     *pow(T9,11.0/3.0)*exp(-12.517/T9);

  return 1.67e16*(15.0/14.0)*Q*f*this->rho*
    this->X[3]*this->X[2]*exp(this->screen(2,8));
}

double Eos::O15AG(void)
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9;

  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0);
  T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0);
  T932=pow(T9,1.5);
 
  Q=3.529;  //  Q value in MeV

  // f is the rate in (cm^3/mol)/s
  f=3.57e11/T923*exp(-39.584/T913-pow(T9/3.000,2.0))*
    (1.0+0.011*T913-0.273*T923-0.020*T9)+
    5.10e+10/T923*exp(-39.584/T913-pow(T9/1.937,2.0))*
    (1.0+0.011*T913+1.59*T923+0.117*T9+1.81*T943+0.338*T953)+
    3.95e-01/T932*exp(-5.849/T9)+1.90e+01*pow(T9,2.85)*
    exp(-7.356/T9-pow(T9/8.000,2.0));

  return 1.67e16*Q*f*this->rho*this->X[4]*
    this->X[2]*exp(this->screen(2,8));
}

double Eos::C14AG(void)
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9, T92, T945;

  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0);
  T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0);
  T932=pow(T9,1.5);
  T92=pow(T9,2.0);
  T945=pow(T9,0.8);

  Q=6.227;  //  Q value in MeV
  
  // f is the rate in (cm^3/mol)/s
  f=(3.375e8/T92)*exp(-32.513/T913)+(1.528e9/T923)
    *exp(-32.513/T913-pow(T9/2.662,2.0))*
    (1.0+0.0128*T913-0.869*T923-0.0779*T9+0.321*T943+0.0732*T953)
    +9.29e-8/T932*exp(-2.048/T9)+2.77e+03/T945*exp(-9.876/T9);
  
  return (9.64e17/(14.0*4.0))*Q*f*this->rho*(this->X[3]+this->X[4])*
    this->X[2]*exp(this->screen(2,6));
}

double Eos::N14PG(void)
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9;

  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0);
  T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0);
  T932=pow(T9,1.5);
 
  Q=7.297;  //  Q value in MeV

  // f is the rate in (cm^3/mol)/s
  f=4.90e+07/T923*exp(-15.228/T913-pow(T9/3.294,2.0))*
    (1.0+0.027*T913-0.778*T923-0.149*T9+0.261*T943+0.127*T953)+
    2.37e+03/T932*exp(-3.011/T9)+2.19e+04*exp(-12.530/T9);

  return 9.64e17*f*Q*this->rho*(this->X[3]+this->X[4])*this->X[1]*exp(this->screen(1,7))/14.0;
}

double Eos::N14AG(void)
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9;

  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0);
  T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0);
  T932=pow(T9,1.5);
 
  Q=4.415;  //  Q value in MeV

  // f is the rate in (cm^3/mol)/s
  f=7.78e+09/T923*exp(-36.031/T913-pow(T9/0.881,2))*
    (1.00+0.012*T913+1.45*T923+0.117*T9+1.97*T943+0.406*T953)
    +2.36e-10/T932*exp(-2.798/T9)+2.03e+00/T932*exp(-5.054/T9)+
    1.15e+04/T923*exp(-12.310/T9);

  return 9.64e17*f*Q*this->rho*this->X[3]*this->X[2]*exp(this->screen(2,7))/(14.0*4.0);
}

double Eos::C12PG(void)
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9;

  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0);
  T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0);
  T932=pow(T9,1.5);

  Q=1.944;  //  Q value in MeV

  // f is the rate in (cm^3/mol)/s
  f=2.04e+07/T923*exp(-13.690/T913-pow(T9/1.500,2.0))*
    (1.00+0.030*T913+1.19*T923+0.254*T9+2.06*T943+1.12*T953)+
    1.08e+05/T932*exp(-4.925/T9)+2.15e+05/T932*exp(-18.179/T9);

  return 9.64e17*f*Q*this->rho*this->X[3]*this->X[1]/**exp(this->screen(1,6))*//12.0;
}

double Eos::Ne19PG(void)
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9, T954;

  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0);
  T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0);
  T932=pow(T9,1.5);
  T954=pow(T9,1.25);

  Q=2.199;  //  Q value in MeV

  // f is the rate in (cm^3/mol)/s
  f=1.71e+06/T923*exp(-19.431/T913)*(1.0+0.021*T913+0.130*T923
     +1.95e-02*T9+3.86e-02*T943+1.47e-02*T953)+1.89e+05/T923
     *exp(-19.431/T913-pow(T9/1.142,2.0))*(1.0+0.021*T913+2.13*T923
     +0.320*T9+2.80*T943+1.07*T953)+8.45e+03/T954*exp(-7.64/T9);

  return 9.64e17*f*Q*this->rho*this->X[5]*this->X[1]*exp(this->screen(1,10))/19.0;
}

double Eos::F17PG(void)
  // This rate is from Felix --- see ~/reactions/reaclib
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9;
  double a1,a2,a3,a4,a5,a6,a7;
  
  T9=0.1*this->T8; T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0); T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0); T932=pow(T9,1.5);
  
  Q=3.92203;  // Mev
  a1=0.172930e+02; a2=0.434069e-03; a3=-0.180950e+02; a4=0.147122e+00;
  a5=-0.154972e+00; a6=0.167096e-01; a7=-0.741063e+00;

  f=exp(a1+a2/T9+a3/T913+a4*T913+a5*T9+a6*T953+a7*log(T9));

  return 9.64e17*Q*f*this->rho*this->X[3]*this->X[1]*exp(this->screen(1,9))/17.0;
}


double Eos::F18PA(void)
  // This rate is from Felix --- see ~/reactions/reaclib
{
  double Q, f;
  double T923, T913, T943, T953, T932, T9;
  double a1,a2,a3,a4,a5,a6,a7;
  
  T9=0.1*this->T8; T913=pow(T9,1.0/3.0);
  T923=pow(T9,2.0/3.0); T943=pow(T9,4.0/3.0);
  T953=pow(T9,5.0/3.0); T932=pow(T9,1.5);
  
  Q=2.88111;  // Mev
  a1=0.251928e+02; a2=-0.393821e+01; a3=0.104823e+02; a4=-0.244073e+02;
  a5=0.170614e+01; a6=-0.112812e+00; a7=0.877482e+01;     

  f=exp(a1+a2/T9+a3/T913+a4*T913+a5*T9+a6*T953+a7*log(T9));

  return 9.64e17*Q*f*this->rho*this->X[3]*this->X[1]*exp(this->screen(1,9))/18.0;
}


double Eos::screen(double Z1, double Z2)
  // evaluates the screening factor according to 
  // Graboske et al. (1973)
{
  double zbar, zhat, gamma12, b, kb, etab, xib, gamma0;
  double c, c1;

  zbar=Ye(); zhat=YZ2();
  zhat=sqrt(zhat+zbar);
 
  gamma0=1.88e-4*sqrt(this->rho*Yi()/pow(this->T8,3.0));
  gamma12=gamma0*Z1*Z2*zhat;

  if (gamma12 < 0.1) {  // weak screening
    c=0.5*zhat*2*Z1*Z2*gamma0;
  } else {
    // first calculate intermediate screening
    etab=this->X[1]+(2.99*this->X[2]/4.0)+
      (this->X[3]*26.72/14.0)+(this->X[4]*26.72/15.0);
    etab/=pow(zhat,0.58)*pow(zbar,0.28);
    xib=pow(Z1+Z2,1.86)-pow(Z1,1.86)-pow(Z2,1.86);
    c1=0.380*etab*xib*pow(gamma0,0.860);
    // now strong screening
    etab=pow(zbar,1.0/3.0);
    xib=pow(Z1+Z2,5.0/3.0)-pow(Z1,5.0/3.0)-pow(Z2,5.0/3.0);
    xib+=0.316*etab*(pow(Z1+Z2,4.0/3.0)-pow(Z1,4.0/3.0)-pow(Z2,4.0/3.0));
    xib+=(0.737/zbar)*(pow(Z1+Z2,2.0/3.0)-pow(Z1,2.0/3.0)-pow(Z2,2.0/3.0))/
      pow(gamma0,2.0/3.0);
    c=0.624*etab*xib*pow(gamma0,2.0/3.0);

    if (gamma12 < 2) c=c1;  // intermediate
    else {
      if (gamma12 < 5) {  
	if (c > c1) c=c1;  // intermediate if its smaller
      }
    }  // default is strong screening
  }
  //  c=0.0; // no screening
  return c;
}


double Eos::fO16O16(void)
  // rate of 16O+16O->32S
  // from Caughlan & Fowler 1988
  // uses Ogata et al. for screening
{

  double T9, T923, T913, T943, f;
  
  T9=0.1*this->T8;
  T913=pow(T9,1.0/3.0);
  T923=T913*T913;
  T943=T923*T923;

  f=7.10e36/T923*exp(-135.93/T913-0.629*T923-0.445*T943+0.0103*pow(T9,2.0));

  // screening from Ogata et al.
  return f*exp(this->ogata(2,2));
  //  return f*exp(this->screen(8.0,8.0));
}

double Eos::O16O16(void)
  // returns ergs/g for O16O16 NB: assumes the oxygen is species 2
{
  double Q, f;
  
  Q=16.542; // MeV
  f=this->fO16O16(); // n_A <sig v>
  
  return this->X[2]*this->X[2]*(9.64e17*Q)*f*this->rho/(2*16.0*16.0);
}



double Eos::fC12O16(void)
  // rate of 12C+16O -> 28Si
  // from Caughlan & Fowler 1988
  // uses Ogata et al. for screening
{

  double T9A, T9A13, T9A56, T9, T932, T9A23, f;
  
  T9=0.1*this->T8;
  T9A=T9/(1.0+0.055*T9);
  T9A13=pow(T9A,1.0/3.0);
  T9A56=pow(T9A,5.0/6.0);     
  T9A23=pow(T9A,2.0/3.0);
  T932=pow(T9,1.5);
 
  f=1.72e31*T9A56/T932*exp(-106.594/T9A13)/(exp(-0.180*pow(T9A,2.0))+
					    1.06e-3*exp(2.562*T9A23));

  // screening from Ogata et al.
  return f*exp(this->ogata(1,2));
  //return f*exp(this->screen(6.0,6.0));
  //  return f;
}

double Eos::C12O16(void)
  // returns ergs/g for C12C12 NB: assumes the carbon is species 1
{
  double Q, f;
  
  Q=16.755; // MeV
  f=this->fC12O16(); // n_A <sig v>
  
  return this->X[1]*this->X[2]*(9.64e17*Q)*f*this->rho/(12.0*16.0);
}


double Eos::fC12C12(void)
  // rate of 12C+12C -> 24Mg
  // from Caughlan & Fowler 1988
  // uses Ogata et al. for screening
{

  double T9A, T9A13, T9A56, T9, T932, f;
  
  T9=0.1*this->T8;
  T9A=T9/(1.0+0.0396*T9);
  T9A13=pow(T9A,1.0/3.0);
  T9A56=pow(T9A,5.0/6.0);     
  T932=pow(T9,1.5);
  
  f=4.27e26*T9A56/T932*exp(-84.165/T9A13-2.12e-03*pow(T9,3.0));

  // screening from Ogata et al.
  return f*exp(this->ogata(1,1));
  //return f*exp(this->screen(6.0,6.0));
  //  return f;
}

double Eos::C12C12(void)
  // returns ergs/g for C12C12 NB: assumes the carbon is species 1
{
  double Q, f;
  
  Q=13.933; // MeV
  f=this->fC12C12(); // n_A <sig v>
  
  return this->X[1]*this->X[1]*(9.64e17*Q)*f*this->rho/(2*144.0);
}

double Eos::ogata(int i, int j)
  // screening factor from Ogata et al 1993
  // eqs. (19)--(21)
{
  double T8, gam, r6, r613, f, lgam;
  double Z1=this->Z[i], Z2=this->Z[j];
  double A1=this->A[i], A2=this->A[j];
  double hf, A, tau;
  
  A=A1*A2/(A1+A2);

  T8=this->T8;
  r6=this->Ye()*this->rho/1e6;
  r613=pow(r6,1.0/3.0);

  hf=0.25*pow(pow(Z1,1.0/3.0)+pow(Z2,1.0/3.0),3.0)/(Z1+Z2);

  gam=0.23*r613/T8; 
  gam*=2.0*Z1*Z2/(pow(Z1,1.0/3.0)+pow(Z2,1.0/3.0));
  lgam=log(gam);

  tau=9.18*pow(Z1*Z2,2.0/3.0)*pow(A/T8,1.0/3.0);
  f=3.0*gam/tau;
  
  double QQ=(1.148-0.00944*lgam-0.000168*lgam*lgam)*gam-0.15625*gam*f*f*hf+
    (-0.18528+0.03863*lgam+0.01095*f)*gam*f*f*f;
  
  if (i==1 && j==1) { 
    // extra terms for 1+1 reaction; doesn't seem to make much differece
    double DD=0.007*(this->X[2]/(this->Yi()*A2))*pow((Z2/Z1)-1.0,2.0);
    double BB=0.456-0.0130*lgam;
    QQ+=DD*(gam*BB*BB/(0.25*hf*2));
    QQ+=DD*gam*f*f*f*(0.18528-0.03863*lgam);
  }
  
  return QQ;
}

double Eos::triple_alpha(void)
  // Fushiki and Lamb's fit to the triple alpha rate
  // includes screening and is good for pycnonuclear regime
  // the traditional formula is
  // 5.3e11*p->rho*p->rho*pow(p->Y/p->T8,3.0)*exp(-44.0/p->T8);
{
  double r6, T6, r613, r616, T613, T623, T653, T632, T612, u;
  double G1, G2, f1, f2, f3, f4, u32;
  
  // this is the simple rate without screening
  // return 5.3e11*this->rho*this->rho*pow(this->X[2]/this->T8,3.0)*exp(-44.0/this->T8);

  r6=1e-6*this->rho*Ye()*2.0; T6=this->T8*100;

  r613=pow(r6,1.0/3.0); r616=sqrt(r613);
  T613=pow(T6,1.0/3.0); T623=T613*T613; T653=T6*T623;
  T632=pow(T6,1.5); T612=pow(T6,0.5);

  u=1.35*r613/T623; u32=pow(u,1.5);

  f1=exp(60.492*r613/T6);
  f2=exp(106.35*r613/T6);
  if (r6 < 5.458e3) f3=exp(-1065.1/T6)/T632; else f3=0.0;
  if (r6 < 1.836e4) f4=exp(-3336.4/T6)/T632; else f4=0.0;
  
  if (u < 1) {
    G1=f1*(f3+16.16*exp(-134.92/T613)/
	   (T623*(pow(1-4.222e-2*T623,2.0)+2.643e-5*T653)));
    G2=f2*(f4+244.6*pow(1+3.528e-3*T623,5.0)*exp(-235.72/T613)/
	   (T623*(pow(1-2.807e-2*T623,2.0)+2.704e-6*T653)));
  } else {
    G1=f1*f3+1.178*(1+(1.0/u32))*exp(-77.554/r616)/
      (T612*(pow(1-5.680e-2*r613,2.0)+8.815e-7*T6*T6));
    G2=f2*f4+13.48*(1+(1.0/u32))*pow(1+5.070e-3*r613,5.0)*exp(-135.08/r616)/
      (T612*(pow(1-3.791e-2*r613,2.0)+5.162e-8*T6*T6));
  }

  r6=r6/(2.0*Ye());
  return 5.12e29*pow(this->X[2],3.0)*r6*r6*G1*G2;
}

// ----------------------------------------------------------------------------
