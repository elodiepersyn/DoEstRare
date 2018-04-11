#include "DoEstRare.h"
#include <stdio.h>


void DoEstRare_stat(double* stat, double* X_, double* Y, double* N_, double* P_, double* autosomal_, double* gender, double* l, double* from_, double* to_, double*  bw_){
  int N=N_[0];
  int P=P_[0];
  int from=from_[0];
  int to=to_[0];
  double bw=bw_[0];
  int autosomal=autosomal_[0];


  //TRANSFORM X INTO A MATRIX
  double** X;
  X=(double**)malloc(N*sizeof(double*));
  int i=0;
  for(i=0; i<N; i++){
    X[i]=(double*)malloc(P*sizeof(double));
  }

  createMatrixFromRVector(X, X_, N, P);


  //COMPUTATION OF CHROMOSOME NUMBERS IN CASES AND CONTROLS
  int Na=0;
  for(i=0; i<N; i++){
    if(Y[i]==1){
      Na=Na+1;
    }
  }

  int chra=0, chru=0;
  if(autosomal==1){
    chra=2*Na;
    chru=2*(N-Na);
  }else{
    for(i=0; i<N; i++){
      if(Y[i]==1){
        chra=chra+gender[i];
      }else{
        chru=chru+gender[i];
      }
    }
  }

  //COMPUTATION OF MUTATION COUNTS IN CASES AND CONTROLS
  double* ma=malloc(P*sizeof(double)); //numeric vector of mutation counts in cases for each variant
  double* mu=malloc(P*sizeof(double)); // numeric vector of mutation counts in controls for each variant
  initialize(ma,P);
  initialize(mu,P);
  colSum_case(ma, X, Y, N, P);
  colSum_ctrl(mu, X, Y, N, P);

  //COMPUTATION OF WEIGHTS
  double* weights=malloc(P*sizeof(double));
  initialize(weights, P);
  for(i=0; i<P;i++){
    //R : pbinom(q=x[1], size=2*Na,prob=(x[2]+1)/(2*Nu+2),lower.tail=T)
    // fonction pbinom dans Rmath.h : double pbinom(double x, double n, double p, int lower_tail, int log_p)
    weights[i]=pbinom(ma[i], chra, (mu[i]+1)/(chru+1), 1, 0);
  }

  //COMPUTATION OF BURDEN COMPONENTS
  double Pa=0;
  double Pu=0;
  for(i=0; i<P;i++){
    Pa=Pa+ma[i]*weights[i];
    Pu=Pu+mu[i]*weights[i];
  }
  Pa=Pa/(chra*P)/sum(weights,P);
  Pu=Pu/(chru*P)/sum(weights,P);

  //COMPUTATION OF DENSITIES
  double * dens_cas=malloc((to-from+1)*sizeof(double));
  double * dens_ctrl=malloc((to-from+1)*sizeof(double));
  double * w=malloc(P*sizeof(double));

  initialize(dens_cas, P);
  initialize(dens_ctrl, P);
  initialize(w, P);

  double sum_ma=0;
  double sum_mu=0;
  sum_ma=sum(ma, P);
  sum_mu=sum(mu, P);

  if(sum_ma!=0){
    proport(w, ma, P);
    densite(l, w, from, to, P, bw, dens_cas);
  }

  if(sum_mu!=0){
    proport(w, mu, P);
    densite(l, w, from, to, P, bw, dens_ctrl);
  }

  //COMPUTATION OF STATISTIC
  double* diff=malloc((to-from+1)*sizeof(double));
  initialize(diff, (to-from+1));
  for(i=0; i<(to-from+1); i++){
    diff[i]=Pa*dens_cas[i]-Pu*dens_ctrl[i];
    if(diff[i]<0){
      diff[i]=-diff[i];
    }
  }
  double* vect_int=malloc((to-from)*sizeof(double));
  initialize(vect_int, to-from);
  for(i=0; i<(to-from); i++){
    vect_int[i]=(diff[i]+diff[i+1])/2;
  }


  stat[0]=sum(vect_int,to-from);

  //printf("%f \n", stat[0]);

  free(ma);
  free(mu);
  free(weights);
  free(w);
  free(diff);
  free(dens_cas);
  free(dens_ctrl);
  free(vect_int);

  for(i=0; i<N; i++){
    free(X[i]);
  }
  free(X);

}


