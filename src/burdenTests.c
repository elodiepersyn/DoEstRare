#include "burdenTests.h"
#include <stdio.h>
#include <stdlib.h>



//------------ Sum Test

void Sum_stat(double* X_, double* Y, double* Ymu, double* Na_, double* Nu_, double* P_, double* stat_){
  int Na=Na_[0];
  int Nu=Nu_[0];
  int N=Na+Nu;
  int P=P_[0];

  //transform X vector into a matrix
  double** X;
  X=(double**)malloc(N*sizeof(double*));
  int i=0;
  for(i=0; i<N; i++){
    X[i]=(double*)malloc(P*sizeof(double));
  }

  createMatrixFromRVector(X, X_, N, P);

  //computation of allele counts per individual
  double* count_i=malloc(N*sizeof(double));
  rowSum(count_i, X, N, P);

  //computation of the test statistic
  double stat_obs;
  stat_obs=single_score_stat(count_i,  Y,  Ymu, N);
  stat_[0]=stat_obs;

  for(i=0; i<N; i++){
    free(X[i]);
  }
  free(X);
  free(count_i);
}


//------------ CAST Test

void CAST_stat(double* X_, double* Y, double* Ymu, double* Na_, double* Nu_, double* P_, double* stat_){
  int Na=Na_[0];
  int Nu=Nu_[0];
  int N=Na+Nu;
  int P=P_[0];

  //transform X vector into a matrix
  double** X;
  X=(double**)malloc(N*sizeof(double*));
  int i=0;
  for(i=0; i<N; i++){
    X[i]=(double*)malloc(P*sizeof(double));
  }

  createMatrixFromRVector(X, X_, N, P);

  //computation of allele counts per individual and CAST genetic score
  double* count_i=malloc(N*sizeof(double));
  rowSum(count_i, X, N, P);
  for(i=0; i<N; i++){
    if(count_i[i]>0){
      count_i[i]=1;
    }
  }

  //computation of the test statistic
  long double stat_obs;
  stat_obs=single_score_stat(count_i,  Y,  Ymu, N);
  stat_[0]=stat_obs;

  for(i=0; i<N; i++){
    free(X[i]);
  }
  free(X);

  free(count_i);
}


//------------ wSum Test

void wSum_stat(double* X_, double* Y, double* Ymu, double* Na_, double* Nu_, double* P_,  double* autosomal_, double* gender, double* Madsen_, double* wbeta_, double* Madsentot_, double* weights, double* stat_){
  int Na=Na_[0];
  int Nu=Nu_[0];
  int N=Na+Nu;
  int P=P_[0];
  int autosomal=autosomal_[0];
  int Madsen=(int)Madsen_[0];
  int wbeta=(int)wbeta_[0];
  int Madsentot=(int)Madsentot_[0];

  //transform X vector into a matrix
  double** X;
  X=(double**)malloc(N*sizeof(double*));
  int i=0;
  for(i=0; i<N; i++){
    X[i]=(double*)malloc(P*sizeof(double));
  }

  createMatrixFromRVector(X, X_, N, P);

  //computation of weighted allele counts per individual
  double* count_i=malloc(N*sizeof(double));

  if(Madsen==1){
    initialize(weights, P);
    weights_wSum(weights, X, Y, N, Nu, autosomal, gender, P );
  }

  if(wbeta==1){
    initialize(weights, P);
    weights_SKAT(weights, X, N, autosomal, gender, P);
  }

  if(Madsentot==1){
    initialize(weights, P);
    weights_wSum_tot(weights, X, N, autosomal, gender, P );
  }

  //printf("%f \t", weights[0]);
  rowSum_weighted(count_i, X, weights, N, P);

  //computation of the test statistic
  long double stat_obs;
  stat_obs=single_score_stat(count_i,  Y,  Ymu, N);

  //adaptive permutations
   stat_[0]=stat_obs;

  for(i=0; i<N; i++){
    free(X[i]);
  }
  free(X);

  free(count_i);
}


//------------ aSum Test

void aSum_stat(double* X_, double* Y, double* Ymu, double* Na_, double* Nu_, double* P_, double* alpha0_, double* stat_){
  int Na=Na_[0];
  int Nu=Nu_[0];
  int N=Na+Nu;
  int P=P_[0];
  double alpha0=alpha0_[0];

  //transform X vector into a matrix
  double** X;
  X=(double**)malloc(N*sizeof(double*));
  int i=0;
  for(i=0; i<N; i++){
    X[i]=(double*)malloc(P*sizeof(double));
  }

  createMatrixFromRVector(X, X_, N, P);

  //computation of weighted allele counts per individual
  double* count_i=malloc(N*sizeof(double));
  double* weight=malloc(P*sizeof(double));

  weights_aSum(weight, X, Y, Ymu, N, P, alpha0 );
  //printf("%f \t", weight[0]);
  rowSum_weighted(count_i, X, weight, N, P);

  //computation of the test statistic
  long double stat_obs;
  stat_obs=single_score_stat(count_i,  Y,  Ymu, N);
  stat_[0]=stat_obs;

  for(i=0; i<N; i++){
    free(X[i]);
  }
  free(X);
  free(count_i);
  free(weight);
}

