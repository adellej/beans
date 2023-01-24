#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include <math.h>

extern "C" {
#include "nrutil.h"
}

#include "useful.h"

double askd(char *text)
{
  double var;
  
  printf("%s...", text);
  scanf("%lg", &var);
  return var;
}

int aski(char *text)
{
  int var;
  
  printf("%s...", text);
  scanf("%d", &var);
  return var;
}

int flines(char *filename)
  // Counts the number of lines in the file 'filename'
{
  int i; 
  char dummy;
  FILE *fp;
  
  fp=fopen(filename,"r");
  
  for (i=0; feof(fp)==0; i++)
    while (fscanf(fp, "%c", &dummy), dummy != '\n');
  
  fclose(fp);
 
  return i-1;
}

void mprintf(char *format, int n, ...)
{
  va_list ap;
  double dummy;
  int i, flag;

  //printf("n=%d\n", n);

  if (n < 0) {
    flag = 0; n=-n;
  } else {
    flag = 1;
  }

  va_start(ap, n);

  //printf("*");
  
  for (i=0; i<n; i++) {
    dummy=va_arg(ap, double);
    printf(format, dummy);
    if (i<n-1) printf(" ");
  }
    
  if (flag == 1) printf("\n");

  va_end(ap);
}

void mfprintf(FILE *fp, char *format, int n, ...)
{
  va_list ap;
  double dummy;
  int i, flag;

  if (n < 0) {
    flag = 0; n=-n;
  } else {
    flag = 1;
  }

  va_start(ap, n);

  for (i=0; i<n; i++) {
    dummy=va_arg(ap, double);
    fprintf(fp, format, dummy);
    if (i<n-1) fprintf(fp, " ");
  }
    
  if (flag == 1) fprintf(fp, "\n");

  va_end(ap);
}


// ---------  sort routine from NR -----------------------------------

#define NRANSI
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort(unsigned long n, double arr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int jstack=0,*istack;
	double a,temp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI


