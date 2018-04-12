#include "functions_burdentests.h"


/*
 * STATS
 */

void rowSum(double* rowsum_vect, double** M, int N, int P){
  int i=0;
  int j=0;
  initialize(rowsum_vect,N);
  for(i=0; i<N; i++){
    for(j=0; j<P; j++){
      rowsum_vect[i]=rowsum_vect[i]+M[i][j];
    }
  }
}

void colSum(double* colsum_vect, double**M, int N, int P){
  initialize(colsum_vect, P);
  int i=0;
  int j=0;
  for(i=0; i<N; i++){
    for(j=0; j<P; j++){
      colsum_vect[j]=colsum_vect[j]+M[i][j];
    }
  }
}

void colSum_group(double* colsum_vect, double**M, double* Y, int N, int P){
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


double mean(double* vect, int size){
  double moy=0;
  moy=sum(vect, size);
  moy=moy/((double)size);
  return moy;
}

double max(double* vect, int size){
  double maximum=vect[0];
  int i=0;
  for(i=1; i<size; i++){
    if(vect[i]>maximum){
      maximum=vect[i];
    }
  }
  return maximum;
}

double sum(double* vect, int size){
  double sum=0;
  int i=0;
  for(i=0; i<size; i++){
    sum=sum+vect[i];
  }
  return sum;
}

double sum_group(double* vect, double* Y, int size){
  double sum=0;
  int i=0;
  for(i=0; i<size; i++){
    if(Y[i]==0){
      sum=sum+vect[i];
    }
  }
  return sum;
}

void rowSum_weighted(double* rowsum_vect, double** M, double* weight, int N, int P){
  int i=0;
  int j=0;

  initialize(rowsum_vect,N);
  for(i=0; i<N; i++){
    for(j=0; j<P; j++){
      rowsum_vect[i]=rowsum_vect[i]+M[i][j]*weight[j];
    }
  }
}

//------------
void MAF_ctl(double* MAF, double** M, double* Y, int N, int Nu, int autosomal, double* gender,  int P){
  initialize(MAF,P);
  colSum_group(MAF, M,  Y, N, P);
  int j=0;
  for(j=0; j<P; j++){
    if(autosomal==1){
      MAF[j]=(MAF[j]+1)/(2*Nu+2);
    }else{
      MAF[j]=(MAF[j]+1)/(sum_group(gender, Y, N)+2);
    }
  }
}

void MAF(double* MAF, double** M, int N, int autosomal, double* gender,  int P){
  initialize(MAF,P);
  colSum(MAF, M, N, P);
  int j=0;
  for(j=0; j<P; j++){
    if(autosomal==1){
      MAF[j]=(MAF[j])/(2*N);
    }else{
      MAF[j]=(MAF[j])/(sum(gender, N));
    }
  }
}

void MAF2(double* MAF, double** M, int N, int autosomal, double* gender,  int P){
  initialize(MAF,P);
  colSum(MAF, M, N, P);
  int j=0;
  for(j=0; j<P; j++){
    if(autosomal==1){
      MAF[j]=(MAF[j]+1)/(2*N+2);
    }else{
      MAF[j]=(MAF[j]+1)/(sum(gender, N)+2);
    }
  }
}



//--------------weights
void weights_wSum(double* weight, double** M, double* Y, int N, int Nu, int autosomal, double* gender, int P){
  initialize(weight,P);
  MAF_ctl(weight, M,  Y, N, Nu, autosomal, gender, P);
  int j=0;
  for(j=0; j<P; j++){
    weight[j]=1/sqrt(N*weight[j]*(1-weight[j]));
  }
}


void weights_wSum_tot(double* weight, double** M, int N, int autosomal, double* gender, int P){
  initialize(weight,P);
  MAF2(weight, M, N, autosomal, gender, P);
  int j=0;
  for(j=0; j<P; j++){
    weight[j]=1/sqrt(N*weight[j]*(1-weight[j]));
  }
}

void weights_SKAT(double* weight, double** M, int N, int autosomal, double* gender, int P){
  initialize(weight,P);
  MAF(weight, M, N, autosomal, gender, P);
  double p=0;
  int j=0;
  for(j=0; j<P; j++){
    p=weight[j];
    weight[j]=(double)dbeta(p, 1, 25, 0);
  }
}


void weights_aSum(double* weight, double** M, double* Y, double* Ymu, int N, int P, double alpha0){
  initialize(weight, P);

  //testing
  int j=0;
  double stat=0;
  double pval;
  double* S=malloc(N*sizeof(double)); // one column of the M matrix
  double U=0;

  for(j=0; j<P; j++){
    selectColumn(S, M, N, P, j);
    U=Uscore_univariate(S, Y, Ymu, N);
    stat=single_score_stat(S, Y, Ymu, N);
    pval=pchisq(stat, 1, 0, 0);
    if(pval<=alpha0 && U<0){
      weight[j]=-1;
    }else{
      weight[j]=1;
    }

  }

  free(S);
}



//score test


double single_score_stat(double* S, double* Y, double* Ymu, int N){
  // computation of Y-Ymu vector
  double* res=malloc(N*sizeof(double));
  int i=0;
  for(i=0; i<N; i++){
    res[i]=Y[i]-Ymu[i];
  }

  //computation of U
  double U=0;
  U= Uscore_univariate(S, Y, Ymu, N);

  //printf("U = %f \n", U);

  //computation of V
  double moyS=mean(S, N);

  long double V1=0;
  long double V2=0;
  long double V=0;

  for(i=0; i<N; i++){
    V1=V1+(res[i]*res[i]);
    V2=V2+((S[i]-moyS)*(S[i]-moyS));
  }
  V=(V1*V2)/(N-1);

  //computation of stat=U^2/V
  long double stat=U*U/V;
  free(res);

  return stat;
}

double Uscore_univariate(double* S, double* Y, double* Ymu, int N){
  // computation of Y-Ymu vector
  double* res=malloc(N*sizeof(double));
  int i=0;
  for(i=0; i<N; i++){
    res[i]=Y[i]-Ymu[i];
  }

  //computation of U
  double U=0;
  for(i=0; i<N; i++){
    U=U+res[i]*S[i];
  }

  free(res);
  return U;
}





//----------initialize a null vector
void initialize(double* vect, int size){
  int i=0;
  for(i=0; i<size; i++){
    vect[i]=0;
  }
}

//---------permutate a vector
void permutate(double* vect, double* vect_perm, int size){
  double random=0;
  int round_random=0;

  int i=0;
  int j=0;
  int ind;

  double* indice=malloc(size*sizeof(double));

  for(i=0; i<size; i++){
    indice[i]=i;
  }

  for(i=0; i<size; i++){
    random=runif(-0.5,(size-1-i)+0.5);
    round_random=round(random);

    ind=indice[round_random];
    vect_perm[i]=vect[ind];

    //update indices
    if(round_random!=(size-1-i)){

      for(j=round_random; j<(size-1-i); j++){
        indice[j]=indice[j+1];
      }
    }
  }

  free(indice);
}

void CopyVector(double* vect, double* vect2, int size){
  int i=0;
  for(i=0; i<size; i++){
    vect2[i]=vect[i];
  }
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

void selectColumn(double* S, double** M, int N, int P, int j){
  int i=0;
  for(i=0; i<N; i++){
    S[i]=M[i][j];
  }
}



//algebra


void product_vector_matrix(double* vect_matrix_prod, double* vect, double** M, int size_vect, int ncol_M){
  initialize(vect_matrix_prod, ncol_M);
  int i=0;
  int j=0;
  for(i=0; i<size_vect; i++){
    for(j=0; j<ncol_M; j++){
      vect_matrix_prod[j]=vect_matrix_prod[j]+vect[i]*M[i][j];
    }
  }
}

void product_matrix_vector(double* vect_matrix_prod, double* vect, double** M, int nrow_M, int size_vect){
  initialize(vect_matrix_prod, nrow_M);
  int i=0;
  int j=0;
  for(i=0; i<nrow_M; i++){
    for(j=0; j<size_vect; j++){
      vect_matrix_prod[i]=vect_matrix_prod[i]+vect[j]*M[i][j];
    }
  }
}


void product_vector_vector1(double* prod, double* vect1, double* vect2, int size){
  initialize(prod, size);
  int i=0;
  for (i=0; i<size; i++){
    prod[i]=vect1[i]*vect2[i];
  }
}



double scalar_product(double* vect1, double* vect2, int size){
  double prod=0;
  int i=0;
  for (i=0; i<size; i++){
    prod=prod+vect1[i]*vect2[i];
  }
  return prod;
}

void MatProduct(double** A, double** B, double** prod, int M, int N, int P){
  //initialization of the matrix
  int m, n, p;
  for(m=0; m<M; m++){
    for(p=0; p<P; p++){
      prod[m][p]=0;
    }
  }

  //multiplication
  for(m=0; m<M; m++){
    for(n=0; n<N; n++){
      for(p=0; p<P; p++){
        prod[m][p]=prod[m][p]+A[m][n]*B[n][p];
      }
    }
  }

}





// code from : https://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html

void Inverse(double** a, int n, double** inv){
  double det=0;
  det=Determinant(a, n);

  CoFactor(a, n, inv);
  Transpose(inv, n);

  int i=0;
  int j=0;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      inv[i][j]=inv[i][j]/det;
    }
  }
}



/*
 Recursive definition of determinate using expansion by minors.
 */
double Determinant(double **a,int n)
{
  int i,j,j1,j2;
  double det = 0;
  double **m = NULL;

  if (n < 1) { /* Error */

  } else if (n == 1) { /* Shouldn't get used */
det = a[0][0];
  } else if (n == 2) {
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  } else {
    det = 0;
    for (j1=0;j1<n;j1++) {
      m = malloc((n-1)*sizeof(double *));
      for (i=0;i<n-1;i++)
        m[i] = malloc((n-1)*sizeof(double));
      for (i=1;i<n;i++) {
        j2 = 0;
        for (j=0;j<n;j++) {
          if (j == j1)
            continue;
          m[i-1][j2] = a[i][j];
          j2++;
        }
      }
      det += pow(-1.0,(double)j1+2.0) * a[0][j1] * Determinant(m,n-1);
      for (i=0;i<n-1;i++)
        free(m[i]);
      free(m);
    }
  }
  return(det);
}

/*
 Find the cofactor matrix of a square matrix
 */
void CoFactor(double **a,int n,double **b)
{
  int i,j,ii,jj,i1,j1;
  double det;
  double **c;

  c = malloc((n-1)*sizeof(double *));
  for (i=0;i<n-1;i++)
    c[i] = malloc((n-1)*sizeof(double));

  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {

      /* Form the adjoint a_ij */
      i1 = 0;
      for (ii=0;ii<n;ii++) {
        if (ii == i)
          continue;
        j1 = 0;
        for (jj=0;jj<n;jj++) {
          if (jj == j)
            continue;
          c[i1][j1] = a[ii][jj];
          j1++;
        }
        i1++;
      }

      /* Calculate the determinate */
      det = Determinant(c,n-1);

      /* Fill in the elements of the cofactor */
      b[i][j] = pow(-1.0,i+j+2.0) * det;
    }
  }
  for (i=0;i<n-1;i++)
    free(c[i]);
  free(c);
}

/*
 Transpose of a square matrix, do it in place
 */
void Transpose(double **a,int n)
{
  int i,j;
  double tmp;

  for (i=1;i<n;i++) {
    for (j=0;j<i;j++) {
      tmp = a[i][j];
      a[i][j] = a[j][i];
      a[j][i] = tmp;
    }
  }
}



/*
 Transpose of a non-square matrix
 */
void Transpose2(double** a, double** transpose, int n, int p)
{
  int i,j;

  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) {
      transpose[j][i]=a[i][j];
    }
  }
}



