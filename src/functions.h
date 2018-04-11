#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>

#define pi 3.14159265359


double GaussianKernel(double u);
void initialize(double* vect, int size);
double sum(double* vect, int size);
void proport(double* prop, double* vect, int size);
void densite(double* l, double* w, int from, int to, int P, double  bw, double* dens);
void createMatrixFromRVector(double** M, double* v, int n_ind, int n_col);

void colSum_case(double* colsum_vect, double**M, double* Y, int N, int P);
void colSum_ctrl(double* colsum_vect, double**M, double* Y, int N, int P);
