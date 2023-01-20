#ifndef __SPLINE_H__
#define __SPLINE_H__

// definitions for fread & fwrite
#define DSIZE sizeof(double)
#define ISIZE sizeof(int)

class Spline {
public:
  double get(double x);
  double getlin(double x);
  double get_x(int i);
  double get_y(int i);
  double startx;
  void init(FILE *fp, int ycol, int ncol, int n);
  void minit(double *x, double *y, int n);
  void tidy();
  int out_of_bounds_flag;
  int log_flag;
  int size(void);
private:
  double *ytab;
  double *derivs;
  double *xtab;
  int num;
  void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
  void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
};

#endif // __SPLINE_H__
