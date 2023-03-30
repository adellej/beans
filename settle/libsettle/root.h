#ifndef __ROOT_H__
#define __ROOT_H__

double zbrent(double (*func)(double), double x1, double x2, double tol);

/* MCU note: is the following function a dead code? never called */
/* There are a few more local functions called from this one, 
   declared nicely inside newt(), see root.c */
void newt(double x[], int n, int *check,
	  void (*vecfunc)(int, double [], double []));

#endif /* __ROOT_H__ */


  
  
