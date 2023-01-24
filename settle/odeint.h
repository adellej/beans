#ifndef __ODEINT_H__
#define __ODEINT_H__

class Ode_Int {
public:
  int ignore, kount, stiff, verbose, tri;
  double dxsav, minstep;
  void init(int n);
  void tidy(void);
  void go(double x1, double x2, double xstep,
	  double eps, void (*derivs)(double, double[],double[]));
  void go_simple(double x1, double x2, int nstep,
		 void (*derivs)(double, double[],double[]));
  void go_scale(double x1, double x2, double step,
		 void (*derivs)(double, double[],double[]));
  void set_bc(int n, double num);
  double get_x(int i);
  double get_y(int n, int i);
  double get_d(int n, int i);
  double* xa(), *ya(int n);
  double *xp, **yp;

private:
  double **dydxp,*hstr,*ystart;
  int kmax,nok,nbad,nvar;
  void rkck(double y[], double dydx[], int n, double x, double h,
	    double yout[],
	    double yerr[], void (*derivs)(double, double[], double[]));
  void rkqs(double y[], double dydx[], int n, double *x, double htry, 
	    double eps,	double yscal[], double *hdid, double *hnext,
	    void (*derivs)(double, double[], double[]));
  void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	      double h1,double hmin, int *nok, int *nbad,
	      void (*derivs)(double, double [], double []));
  void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
	     void (*derivs)(double, double [], double []));
  void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep,
	void (*derivs)(double, double [], double []));
  void rkscale(double vstart[], int nvar, double x1, double x2, double h1,
	void (*derivs)(double, double [], double []));
  double **d,*x;   // from stifbs.c
  
  void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	     double xs, double htot, int nstep, double yout[],
	     void (*derivs)(double, double [], double []));
  void trisimpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	     double xs, double htot, int nstep, double yout[],
	     void (*derivs)(double, double [], double []));
  void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	      double yscal[], double *hdid, double *hnext,
	      void (*derivs)(double, double [], double []));
  void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  void tridag(double a[], double b[], double c[], double r[], double u[],
	      unsigned long n);  
};

#endif // __ODEINT_H__
