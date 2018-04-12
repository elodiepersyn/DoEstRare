#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>



double sum(double* vect, int size);   //sum of a vector
void initialize(double* vect, int size); //initialize a vector with null values
void createMatrixFromRVector(double** M, double* v, int n_ind, int n_col);