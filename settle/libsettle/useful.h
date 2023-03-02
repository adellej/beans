#ifndef __USEFUL_H__
#define __USEFUL_H__

int flines(char *filename);
void mprintf(char *format, int n, ...);
void mfprintf(FILE *fp, char *format, int n, ...);
double askd(char *text);
int aski(char *text);

#endif // __USEFUL_H__
