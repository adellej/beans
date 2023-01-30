// settle.cc
//
// Settling solutions
// Integrates the thermal profile of a settling atmosphere
// to find ignition conditions
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

extern "C" {
#include "nrutil.h"
#include "root.h"
}

#include "odeint.h"
#include "eos.h"
#include "spline.h"
#include "useful.h"

#define me 510.999      // electron mass in keV
#define KAPPAT 0.22     // opacity at the top
#define F3a 1.9         // enhancement for triple alpha
#define EH 6.03e18       // energy release from H --> He in hot CNO
#define F14 0.352       // mass fractions of CNO elements
#define F15 0.648       // in 14 or 15O in hot CNO cycle
#define epsH 5.94e15    //erg/g/s used to be 5.8e15

// the following avoids name mangling, which will be nice in python:
// instead of _Z6mainerPdS_S_S_PiS_ you can type directly mainer
extern "C" {
  int mainer(double* flu, double* Z, double* X, double* mdo, int* docomp,
	     double* trec, double* alpha, double* fluen, double* radius, double*mass);
}

// -------- Global variables ----------

Eos EOS;

struct {
  int OCEAN; // =0 if atmosphere, =1 if ocean or crust
  double R; // radius of NS in cm
  double M;
  double ZZ;
  double g; // gravity
  double X; // initial hydrogen abundance
  double Z; // initial mass fraction of CNO elements
  double Y; // column -- used for root finder
  double mdot, avmdot; // accretion rate
  double yd; // the depletion column
  double Fb,Qb;  // flux at the base
  double Fermi_n, Fermi_alpha; // for Fermi integrals
  double ycrust, yb, yt, Tt, Ft;
  int COMPRESS;
  int debug;
  double Medd;
  //  int output_flag;
} G;

int SWITCH;

struct {
  FILE *out, *out2, *dat, *yb, *out3, *flux, *ign, *modes, *mag, *datm;
} fp;

Ode_Int ODE;
int VERBOSE;

// --------- Declarations ------------
double find_rho_eqn(double rho);
double find_rho(void);
void derivs(double y, double ff[], double dfdy[]);
double doint(double yb);
double find_yb(void);
void output(int i);
double find_F(void);
double dointF(double F);
void jacobn(double, double *, double *, double **, int) {};

//------------------------------------------------------------------------


int mainer(double* flu, double* Z, double* X, double* mdo, int* docomp,
	   double* trec, double* alpha, double* fluen, double* radius, double*mass)
{
  int flag, n, i;
  double yb, y, dummy, Xbar;


  // ------ set up output files ----------------
  //  fp.out=fopen("out/settle","w");
  //fp.ign=fopen("cube","a");


  // ---------- Initialise parameters --------------------------------------

  //  G.g=2.45; //printf("Enter gravity (g14)..."); scanf("%lg", &G.g);

  G.M = *mass;
  G.M *= 1.989e33; // convert to grams

  G.R = *radius; //printf("Enter radius (km)..."); scanf("%lg", &G.R); G.R*=1e5;
  //printf("Radius = %lg km\n", G.R/1e5);
  G.R *= 1e5; // convert to cm

  G.ZZ  =pow((1.0 - ((2*G.M*6.67428e-8)/(G.R*pow(2.9979e10,2)))),-0.5);

  G.g = G.M*G.ZZ*6.67428e-8/pow(G.R,2);  //printf("Gravity g=%lg\n", G.g);


  //EOS.Q = 900.0; printf("Q=%lg in the crust\n", EOS.Q);

  // Parameters given on command line

  G.Fb = *flu;
  //printf("Setting base flux = %lg\n",G.Fb);

  G.Z = *Z;
  //printf("Setting metallicity = %lg\n",G.Z);

  G.X = *X;
  if (G.X == 0.0) G.X = 1e-10;
  //printf("Setting accreted H fraction X=%lg\n", G.X);

  G.mdot = *mdo;
  if (G.mdot > 10.0) G.mdot /= G.Medd;
  //printf("I get mdot/Edd=%lg\n", G.mdot);
  G.Medd=(1.75*(1.7/(1+G.X))*(1e-8)*(5.01837638e24))/(G.R*G.R);

  G.yd =G.X*G.mdot*G.Medd*EH/(epsH*G.Z);
  //printf("Depletion column=%lg\n", G.yd);

  G.COMPRESS = *docomp;
  //if (G.COMPRESS) printf("Including compressional heating.\n");
  //else printf("No compressional heating!\n");

  //    G.output_flag=atoi(argv[6]);

  for (int ii=1; ii<=1; ii++) {

    // G.avmdot=0.1;
    G.avmdot = G.mdot * 1.0;

    G.Qb = G.Fb;

    if (G.Fb==0.0) {
    	G.Fb=1e-6; printf("Can't have zero, so setting F=1e-6\n");
    }
    if (G.Fb < 1e10) G.Fb*=14.62*G.avmdot*5.8e21;  // convert to cgs
    // factor of 2 converts back to time-average mdot (relevant for crust)
    //printf("in cgs, Fb=%lg\n", G.Fb);

    //    G.Fb=0.1*ii*14.62*G.mdot*5.8e21;  // convert to cgs
    //    if (ii==2) G.debug=1;

    //    printf("\n-----------------------------\nFLUX=%lg, %lg\n",
    //   0.1*ii, G.Fb);


     EOS.init(4);
     EOS.A[1]=1.0; EOS.Z[1]=1.0;  // H
     EOS.A[2]=4.0; EOS.Z[2]=2.0;  // He
     EOS.A[3]=14.0; EOS.Z[3]=8.0;  // 14O
     EOS.A[4]=15.0; EOS.Z[4]=8.0;  // 15O

     // ------ Do integration through the atmosphere ------------------------
     //printf("\nSearching for ignition depth...\n");
     ODE.init(3);
     yb=find_yb();

     // ------ Loop through the results --------------------------------------

     flag=0; // flag used to find base of hydrogen layer

     // loop
     for (i=1; i<=ODE.kount; i++) {

	  // column depth
	  y=ODE.get_x(i);
	  // Calculate composition at this depth
	  EOS.X[1]=G.X*(1.0-y/G.yd);  if (EOS.X[1] < 0.0) EOS.X[1]=0.0;
	  EOS.X[3]=F14*G.Z; EOS.X[4]=F15*G.Z; EOS.X[2]=1.0-G.Z-EOS.X[1];
	  // Output
	  output(i);

	  // Conditions at base of H column
	  if (y > G.yd && flag == 0) {
		flag=1;
		//printf("at base of the H column,");
		//printf(" y=%lg, rho=%lg, T=%lg, rhoYe=%lg (X=%lg)\n",
		//       y, EOS.rho, 1e8*EOS.T8, EOS.rho*EOS.Ye(), EOS.X[1]);
	  }
     }

  // ----------- Summarize ignition conditions --------------------------
  // on screen
  //printf("\n------------- Ignition conditions ----------------------------\n");
  //printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
  //       "Z", "mdot", "T", "y", "P", "Y", "X", "rho", "flux");
  //mprintf("%10.4lg", 9,
  //	 EOS.X[3] + EOS.X[4], G.mdot * 8.8e4, EOS.T8 * 1e8, yb, G.g * yb,
  //	 EOS.X[2], EOS.X[1], EOS.rho, ODE.get_y(2, 1));
  //printf("---------------------------------------------------------------\n");

  // to file
  // mfprintf(fp.ign, "%10lg", -13,
  //  EOS.X[3]+EOS.X[4], G.mdot, 1e8*EOS.T8, y, EOS.X[2],
  //  EOS.X[1], EOS.rho, y/(3600*G.mdot*8.8e4),
  //  EOS.YZ2()/EOS.Ye(), EOS.Ye(), G.X, G.Fb/(14.62*G.mdot*5.8e21),
  //   EOS.C14AG());



  //printf("z=%lg\n", ODE.get_y(3,ODE.kount));
  yb=0.65*yb;
  *trec  = G.ZZ * yb / (3600.0 * G.Medd * G.mdot);

  //printf("Recurrence time =%lg hours\n", *trec);

  if (EOS.X[1]==0.0) Xbar=0.5*G.X*G.yd/y;
  else Xbar=G.X*(1.0-0.5*y/G.yd);
  //printf("Xbar=%lg, Q=%lg, Energy=%lg\n", Xbar, 1.6+4.0*Xbar,
  //       4*PI*G.R*G.R*y*9.64e17*(1.6+4.0*Xbar)/ZZ);

  *alpha = 290. / (1.35 + 6.05 * Xbar);

  *fluen = (4*M_PI*G.R*G.R*y*9.64e17*(1.35+6.05*Xbar)/G.ZZ)/1e39; //units of 1e39 erg/g

  //printf("eps (14C+alpha) is %lg\n", EOS.C14AG());
  //printf("eps (3a) is %lg\n", EOS.triple_alpha());

  //printf("kappa=%lg\n", EOS.opac());

  //printf("t_alpha=%lg hours\n", 5.84e17*EOS.X[2]/(3600.0*EOS.triple_alpha()));

  //mfprintf(fp.ign, "%lg", -12,
  //	    G.Qb, G.Z, G.X, G.mdot, 1e-8*y, EOS.T8, EOS.X[1], EOS.X[2],
  //	    1.6+4.0*Xbar, ZZ*yb/(3600.0*8.8e4*G.mdot),
  //	    1e-39*4*PI*G.R*G.R*y*9.64e17*(1.6+4.0*Xbar)/ZZ,
  //	    290.0/(1.6+4.0*Xbar));
  //fprintf(fp.ign, "\n");


  /*
  // ---------- Now integrate deeper to the ocean --------------------
  printf("\nNow integrating into ocean...\n");

  // composition
  //  EOS.tidy(); EOS.init(1);
  //  EOS.X[1]=1.0; EOS.A[1]=60.0; EOS.Z[1]=30.0; // set composition
  EOS.X[1]=1.0; EOS.A[1]=12.0; EOS.Z[1]=6.0;
  EOS.X[2]=0.0; EOS.X[3]=0.0; EOS.X[4]=0.0;
  printf("Z=%lg, A=%lg, Ye=%lg, YZ2=%lg\n",
	 EOS.Z[1], EOS.A[1], EOS.Ye(), EOS.YZ2());

  // integrate
  ODE.tidy(); ODE.init(4);
  //  EOS.T8=5.0; // set to steady-state temperature
  ODE.set_bc(1,EOS.T8*1e8);  // temperature at top of ocean set by ignition
  ODE.set_bc(2,G.Fb); // constant flux underneath layer
  ODE.set_bc(3,0.0); ODE.set_bc(4,0.0); // height and mom of inertia
  G.OCEAN=1; ODE.go(yb, 1e11, 0.1*yb, 1e-8, derivs);

  // loop through results and output
  G.ycrust=0.0;
  for (i=1; i<=ODE.kount; i++) output(i);
*/


  ODE.tidy();
  EOS.tidy();


  return 0;
  } // loop through flux


  // ---------- tidy up ---------------------------------------------
  //  fclose(fp.out);
  //fclose(fp.ign);

  /// added by MC - to eliminate warning about missing retval
  return 1;
}


// ---------------------------------------------------------------------

void output(int i)
  // Output results from ODE, index i
{
  double y, T, kappa, z, I, deg, eta, dummy, lam_ei, x, gam;
  double beta, taustar, sig, udrift, h, tdiff;
  double Bcrit, Beq, AA, etaK, Bcrit2;

  y=ODE.get_x(i); T=ODE.get_y(1,i); EOS.T8=1e-8*T;
  // the composition has been set already
  G.Y=y; EOS.rho=find_rho();
  //kappa=EOS.opac();
  //z=ODE.get_y(3,i)-ODE.get_y(3,ODE.kount);
  //I=ODE.get_y(4,i)-ODE.get_y(4,ODE.kount)-8*PI*1e24*
  // 4e-6*ODE.get_y(3,ODE.kount)*(ODE.get_x(ODE.kount)-ODE.get_x(i))/3.0;
  //deg=EOS.pe()/(8.254e15*EOS.rho*EOS.T8*EOS.Ye()); // ratio of Pe/Pe(ideal gas)
  //  eta=(EOS.Chabrier_EF()-510.999)/(8.617*EOS.T8);
  //  h=y/EOS.rho;

  //  x=EOS.x(); beta=x/sqrt(1+x*x);
  //gam=0.49*(EOS.YZ2()/EOS.Yi())*pow(EOS.rho*1e-7*EOS.Yi(),1.0/3.0)/EOS.T8;
  //lam_ei=log(pow(2*PI*EOS.Ye()/(3*EOS.Yi()),1.0/3.0)*sqrt(1.5+3.0/gam));
  //lam_ei-=0.5*x/(1+x);

  // Brunt frequency
  // AA=-EOS.chi(&EOS.T8)*EOS.rho/(y*EOS.chi(&EOS.rho));
  //AA*=EOS.del_ad()-1e-8*y*ODE.get_d(1,i)/EOS.T8;

  //  if (EOS.gamma()<173.0) {
  // G.ycrust=y;
  //}

  // Output: y, T, rho, X, z, I, ef, eta, t_accr, flux, lam_ei
  //  mfprintf(fp.out, "%10lg", 13,
  //   y, T, EOS.rho, EOS.X[1], z, I, EOS.Chabrier_EF(), eta,
  //   y/(G.mdot*8.8e4), ODE.get_y(2,i), lam_ei,
  //   EOS.C14AG(), EOS.triple_alpha());

}


// -------------- Find ignition depth by iteration ------------------


double find_F(void)
{
  return zbrent(dointF, 1.001*G.Ft, 1.3*G.Ft, 1e-3);
}

double dointF(double F)
{
  double Tb, x1, x2, heat, cool;
  int i;

  // upper boundary
  //yt=1e3; Tt=pow(1.53e19*(KAPPAT/0.2)*G.Z*yt*(G.yd-0.5*yt),0.25);
  //if (G.yd < yb) Ft=G.Fb+5.8e15*G.Z*(G.yd-yt);
  //else Ft=G.Fb+5.8e15*G.Z*(yb-yt);

  ODE.set_bc(1,G.Tt);
  ODE.set_bc(2,F);
  ODE.set_bc(3,0.0);

  // do integration
  G.OCEAN=0;
  ODE.go(G.yt, G.yb, G.yt, 1e-8, derivs);

  if (G.COMPRESS) printf("Tried F=%lg, base flux = %lg\n", F, ODE.get_y(2,ODE.kount));

  return ODE.get_y(2,ODE.kount)-G.Fb;
}



double find_yb(void)
  // find ignition depth
{
  double yb=zbrent(doint,8.0,10.0,1e-3);
  return pow(10.0,yb);
}

double doint(double yb)
  // integrates the atmosphere
{
  double F, Tb, x1, x2, heat, cool;

  yb=pow(10.0,yb);

  // upper boundary
  G.yt=1e3;
  if (G.yd < yb) G.Ft=G.Fb+epsH*G.Z*G.yd;
  else G.Ft=G.Fb+epsH*G.Z*yb;

  //  G.Tt=pow(1.53e19*(KAPPAT/0.2)*G.Z*G.yt*(G.yd-0.5*G.yt),0.25);
  G.Tt=pow(2650.0*G.Ft*G.yt,0.25);

  //ODE.set_bc(1,Tt);
  //ODE.set_bc(2,Ft);
  //ODE.set_bc(3,0.0);
  //ODE.set_bc(4,0.0);
  // do integration
  //G.OCEAN=0;
  //ODE.go(yt, yb, yt, 1e-8, derivs);
  G.yb=yb;
  if (G.COMPRESS) F=find_F();
  else dointF(G.Ft);
  if (G.COMPRESS) printf("flux/Fb= %lg; flux/Ft=%lg\n", F/G.Fb, F/G.Ft);

  // compare heating with cooling
  Tb=ODE.get_y(1,ODE.kount); EOS.T8=1e-8*Tb;

  //EOS.X[1]=G.X*(1.0-G.Z*yb/(6.9e7*G.mdot)); if (EOS.X[1] < 0.0) EOS.X[1]=0.0;
  EOS.X[1]=G.X*(1.0-yb/G.yd); if (EOS.X[1] < 0.0) EOS.X[1]=0.0;

  EOS.X[3]=F14*G.Z; EOS.X[4]=F15*G.Z; EOS.X[2]=1.0-G.Z-EOS.X[1];
  G.Y=yb; EOS.rho=find_rho();

  // cooling piece
  x1=EOS.opac(); EOS.T8*=1.001; EOS.rho=find_rho();
  x2=EOS.opac(); EOS.T8/=1.001; EOS.rho=find_rho();
  cool=7.564e-5*pow(Tb,3.0)/(EOS.opac()*yb*yb)*
    (4-(log(x2/x1)/log(1.001)));

  // triple alpha
  x1=EOS.triple_alpha(); EOS.T8*=1.001; EOS.rho=find_rho();
  x2=EOS.triple_alpha(); EOS.T8/=1.001; EOS.rho=find_rho();
  heat=EOS.triple_alpha()*(log(x2/x1)/log(1.001))/(1e8*EOS.T8);
  //if (EOS.X[1]>1e-3) heat*=F3a;
  //smoothly turn on F3a between X=0 and X=1:
  if (EOS.X[1] > 0.143)heat*=F3a;
  else{
  heat*=(F3a-1.0)*(EOS.X[1]/0.143) + 1.0;
  }


  //printf("yb=%lg, Tb=%lg, heat-cool/cool=%lg\n", yb, Tb, (heat-cool)/cool);

  return (heat-cool)/cool;
 }



// --------------------------- Derivatives -----------------------------

void derivs(double y, double ff[], double dfdy[])
  // Evaluate derivatives
  // ff[1]=T, ff[2]=F, ff[3]=z
{
  double eps;

  // catch negative temperatures
  if (ff[1]<0.0) ff[1]=1e7;

  // find density
  G.Y=y; EOS.T8=ff[1]*1e-8;
  if (G.OCEAN==0) { // calculate composition if in atmosphere
    EOS.X[1]=G.X*(1.0-y/G.yd);  if (EOS.X[1]<0.0) EOS.X[1]=0.0;
    EOS.X[2]=1.0-EOS.X[1]-G.Z; EOS.X[3]=G.Z*F14; EOS.X[4]=G.Z*F15;
  }
  EOS.rho=find_rho();

  eps=0.0;
  if (G.OCEAN == 0 && EOS.X[1] > 0.0) eps=epsH*G.Z;
  if (G.COMPRESS) eps+=EOS.CP()*G.mdot*G.Medd*(ff[1]*EOS.del_ad()/y-dfdy[1]);

  // heat equation
  dfdy[1]=3305.1*ff[2]*EOS.opac()/pow(ff[1],3.0);
  // flux
  dfdy[2]=-eps;
  // hydrostatic balance
  dfdy[3]=-1.0/EOS.rho;
}




// ---------------------- root finder for density  ---------------------

double find_rho(void)
{
  double old, found;
  old=EOS.rho;
  found = zbrent(find_rho_eqn,10.0,1e11,1e-8);
  EOS.rho=old;
  return found;
}

double find_rho_eqn(double rho)
{
  EOS.rho=rho;
  return (EOS.ptot()-G.g*G.Y);
}
