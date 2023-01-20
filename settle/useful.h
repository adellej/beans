#ifndef USEFUL_H
#define USEFUL_H

int flines(char *filename);
void mprintf(char *format, int n, ...);
void mfprintf(FILE *fp, char *format, int n, ...);
double askd(char *text);
int aski(char *text);

#endif // USEFUL_H
