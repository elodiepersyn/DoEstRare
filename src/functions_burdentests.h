#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>
#include "functions_shared.h"

//stats

void rowSum(double* rowsum_vect, double**M, int N, int P); //computes row sums of the matrix M (N x P) into the vector rowsum_vect
void colSum(double* colsum_vect, double**M, int N, int P);  //computes col sums of the matrix M (N x P) into the vector colsum_vect
void colSum_group(double* colsum_vect, double**M, double* Y, int N, int P); //computes col sums of the matrix M (N x P) into the vector colsum_vect for the controls specified into the Y phenotype vector

double mean(double* vect, int size);  //mean of a vector
double max(double* vect, int size);   //max of a vector

double sum_group(double* vect, double* Y, int size);  //sum of the vector for controls

void rowSum_weighted(double* vect, double** M, double* weight, int N, int P); //row sum for the matrix M (N x P), with weights for each colum of M.

void MAF_ctl(double* MAF, double** M, double* Y, int N, int Nu, int autosomal, double* gender, int P); //computes the MAF in controls
void MAF(double* MAF, double** M, int N, int autosomal, double* gender,  int P);  //computes the overall MAF
void MAF2(double* MAF, double** M, int N, int autosomal, double* gender,  int P);

void weights_wSum(double* weight, double** M, double* Y, int N, int Nu, int autosomal, double* gender, int P);  //weight computation for the wSum test
void weights_aSum(double* weight, double** M, double* Y, double* Ymu, int N, int P, double alpha0);  //weight computation for the aSum test
void weights_SKAT(double* weight, double** M, int N, int autosomal, double* gender, int P); //weight computation for the SKAT test
void weights_wSum_tot(double* weight, double** M, int N, int autosomal, double* gender, int P);


double single_score_stat(double* S, double* Y, double* Ymu, int N); //computes the score test statistic (U^2/V) to test the association between S and Y  ( logit(P(Y=1))=a+bS )
double Uscore_univariate(double* S, double* Y, double* Ymu, int N); //computes the score U

//vector

void CopyVector(double* vect, double* vect2, int size);

//matrix

void selectColumn(double* S, double** M, int N, int P, int j);

//algebra
void product_vector_matrix(double* vect_matrix_prod, double* vect, double** M, int size_vect, int size_M);
void product_matrix_vector(double* vect_matrix_prod, double* vect, double** M, int nrow_M, int size_vect);

double scalar_product(double* vect1, double* vect2, int size);
void product_vector_vector1(double* prod, double* vect1, double* vect2, int size);

void MatProduct(double** A, double** B, double** prod, int M, int N, int P);

double Determinant(double **a,int n);
void CoFactor(double **a,int n,double **b);
void Inverse(double** a, int n, double** inv);

void Transpose(double **a,int n);
void Transpose2(double** a, double** transpose, int n, int p);
