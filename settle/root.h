#ifndef __ROOT_H__
#define __ROOT_H__

#define float double

float zbrent(float (*func)(float), float x1, float x2, float tol);

/* MCU note: is the following function a dead code? never called */
/* There are a few more local functions called from this one, 
   declared nicely inside newt(), see root.c */
void newt(float x[], int n, int *check,
	  void (*vecfunc)(int, float [], float []));

#undef float

#endif /* __ROOT_H__ */


  
  
