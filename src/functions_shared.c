#include "functions_shared.h"





double sum(double* vect, int size){
  double sum=0;
  int i=0;
  for(i=0; i<size; i++){
    sum=sum+vect[i];
  }
  return sum;
}


void initialize(double* vect, int size){
  int i=0;
  for(i=0; i<size; i++){
    vect[i]=0;
  }
}


void createMatrixFromRVector(double** M, double* v, int n_ind, int n_col){
  int i=0;
  int j=0;
  for(i=0; i<n_ind; i++){
    for(j=0; j<n_col; j++){
      M[i][j]=v[n_ind*j+i];
    }
  }
}