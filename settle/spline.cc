// class Spline
//
// Implements interpolation using splines as in numerical recipes 
// (the routines spline and splint are taken straight from there).
//
// The data should be in the form of (x, y_i) rows, in binary 
// (double precision) with no seperators between the rows. 
// Spline::init reads x (assumed to be in column
// 1) and y (from column ycol). There are ncol columns and n rows.
//
// 25.8.97
// Added 'out_of_bounds_flag'. If the user sets this variable to zero,
// then 'get' returns zero if the x-coord is out of bounds; otherwise (default)
// it returns the first (or last) value in the table
// 2.3.98
// Added 'log_flag'. If this variable is set by the user to a non-zero
// value, the tabulated values are interpreted as log_10 values. Then 
// get returns 10^{table value}.
// (Fixed a bug in this, 29.4.98)
// NB check that the x's are in ascending order
//
// dec 12th 1999 made it reset the file pointer to the beginning of
// the file + one integer since I do that every time anyway!
// ie. includes the fseek(fp,ISIZE,0) automatically
//
// 9th May 2000 
// added getlin(x) which does linear interpolation rather than spline
//
// 30th Sep 2001
// added size() which returns the number of data points in the table
// (ie. it returns the value of this->num)
// 2nd Oct 2001
// added get_x and get_y which allows access to the table
//

#include <stdio.h>
#include <math.h>

extern "C" {
#include "nrutil.h"
}

#include "spline.h"

double Spline::get_x(int i)
{
  return this->xtab[i];
}

double Spline::get_y(int i)
{
  return this->ytab[i];
}

int Spline::size(void)
{
  return this->num;
}

double Spline::get(double x)
{
  double u;

  //printf("num is %d\n", this->num);

  if (x <= this->xtab[1]) {
    //printf("Now in here  1...%lg, %lg\n",x,this->xtab[1]);
    if (this->out_of_bounds_flag == 0) return 0.0;
    else {
      if (log_flag==0.0) return this->ytab[1];
      else return pow(10.0,this->ytab[1]); 
    }
  }

  if (x >= this->xtab[this->num]) {
    //printf("In here  2...%lg, %lg\n",x,this->xtab[this->num]);
    if (this->out_of_bounds_flag == 0) return 0.0;
    else {
      if (log_flag==0.0) return this->ytab[this->num];
      else return pow(10.0,this->ytab[this->num]); 
    }
  }

  //printf("I'm here, flag=%d, %lg, %lg, %lg\n",this->out_of_bounds_flag,
  //this->xtab[1], x, this->xtab[this->num]);

  splint(this->xtab,this->ytab,this->derivs,this->num,x,&u);

  if (log_flag==0.0) return u;
  else return pow(10.0, u);
}

double Spline::getlin(double x)
{
  double u;
  int i;

  if (x <= this->xtab[1]) {
    if (this->out_of_bounds_flag == 0) return 0.0;
    else {
      if (log_flag==0.0) return this->ytab[1];
      else return pow(10.0,this->ytab[1]); 
    }
  }

  if (x >= this->xtab[this->num]) {
    if (this->out_of_bounds_flag == 0) return 0.0;
    else {
      if (log_flag==0.0) return this->ytab[this->num];
      else return pow(10.0,this->ytab[this->num]); 
    }
  }

  i=1;
  while (this->xtab[i]<x) i++;
  u=this->ytab[i-1]+(this->ytab[i]-this->ytab[i-1])*(x-this->xtab[i-1])/
    (this->xtab[i]-this->xtab[i-1]);

  if (log_flag==0.0) return u;
  else return pow(10.0, u);
}


void Spline::minit(double *x, double *y, int n)
  // initialize from memory rather than a file
{
  double d1,d2;
  int i;
  
  this->num=n;
  this->xtab=dvector(1,n);
  this->ytab=dvector(1,n); 
  this->derivs=dvector(1,n); 

  for(i=1; i<=n; i++) {  // copy x and y arrays
    this->xtab[i]=x[i];
    this->ytab[i]=y[i];
  }

  /* These are our best guesses at the first derivatives at the
     start and end points */
  d1=(this->ytab[2]-this->ytab[1])/(this->xtab[2]-this->xtab[1]);
  d2=(this->ytab[n]-this->ytab[n-1])/(this->xtab[n]-this->xtab[n-1]);
  spline(this->xtab,this->ytab,n,d1,d2,this->derivs);

  // default handling of out of bounds
  this->out_of_bounds_flag = 1;
  // default not logs
  this->log_flag=0;

  // 1st x-coordinate
  this->startx=this->xtab[1];
}


void Spline::init(FILE *fp, int ycol, int ncol, int n)
{
  double d1,d2,dummy;
  int i,col;

  this->num=n;
  this->xtab=dvector(1,n);
  this->ytab=dvector(1,n); 
  this->derivs=dvector(1,n); 

  // reset to beginning of file
  fseek(fp,ISIZE,0);

  for(i=1; i<=n; i++) {
    fread(&(this->xtab[i]),1,DSIZE,fp);
    col=1;
    while (++col != ycol) fread(&dummy,1,DSIZE,fp);
    fread(&(this->ytab[i]),1,DSIZE,fp);
    while (col++ != ncol) fread(&dummy,1,DSIZE,fp);
    //    printf("x=%lg, y=%g\n", this->xtab[i], this->ytab[i]);
  }

  /* These are our best guesses at the first derivatives at the
     start and end points */
  d1=(this->ytab[2]-this->ytab[1])/(this->xtab[2]-this->xtab[1]);
  d2=(this->ytab[n]-this->ytab[n-1])/(this->xtab[n]-this->xtab[n-1]);
  spline(this->xtab,this->ytab,n,d1,d2,this->derivs);

  // default handling of out of bounds
  this->out_of_bounds_flag = 1;
  // default not logs
  this->log_flag=0;

  // 1st x-coordinate
  this->startx=this->xtab[1];

}

void Spline::tidy()
{
  free_dvector(this->ytab,1,this->num);
  free_dvector(this->derivs,1,this->num);
}

void Spline::spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
  int i,k;
  double p,qn,sig,un,*u;
  
  u=dvector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_dvector(u,1,n-1);
}

void Spline::splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
  int klo,khi,k;
  double h,b,a;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) nrerror("Bad xa input to routine splint");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}





