#include "functions.h"

void densite(double* l, double* w, int from, int to, int P, double  bw, double* dens){
  int pos=0;
  int j=0;
  initialize(dens, (to-from+1));
  for(pos=from; pos<=to; pos++){
    dens[(pos-1)]=0;
    for(j=0;j<P; j++){
      dens[(pos-1)]=dens[(pos-1)]+w[j]*GaussianKernel((double)(pos-l[j])/bw)/bw;
    }
  }
}

double GaussianKernel(double u){
  double k;
  k=1/sqrt(2*pi)*exp(-0.5*u*u);
  return k;
}

void initialize(double* vect, int size){
  int i=0;
  for(i=0; i<size; i++){
    vect[i]=0;
  }
}

void proport(double* prop, double* vect, int size){
  initialize(prop, size);
  int i=0;
  for(i=0; i<size; i++){
    prop[i]=vect[i]/sum(vect,size);
  }
}

double sum(double* vect, int size){
  int i=0;
  double somme=0;
  for(i=0; i<size; i++){
    somme=somme+vect[i];
  }
  return somme;
}


//------------create a matrix from a vector
void createMatrixFromRVector(double** M, double* v, int n_ind, int n_col){
  int i=0;
  int j=0;
  for(i=0; i<n_ind; i++){
    for(j=0; j<n_col; j++){
      M[i][j]=v[n_ind*j+i];
    }
  }
}

void colSum_ctrl(double* colsum_vect, double**M, double* Y, int N, int P){
  initialize(colsum_vect, P);
  int i=0;
  int j=0;
  for(i=0; i<N; i++){
    for(j=0; j<P; j++){
      if(Y[i]==0){
        colsum_vect[j]=colsum_vect[j]+M[i][j];
      }
    }
  }
}

void colSum_case(double* colsum_vect, double**M, double* Y, int N, int P){
  initialize(colsum_vect, P);
  int i=0;
  int j=0;
  for(i=0; i<N; i++){
    for(j=0; j<P; j++){
      if(Y[i]==1){
        colsum_vect[j]=colsum_vect[j]+M[i][j];
      }
    }
  }
}
