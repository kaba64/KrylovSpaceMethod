#include <math.h>

void MVmulti(int nrows, int* ia, int* ja, double* a, double* x, double* q){
  for(int row=0;row<nrows;row++){
    double dot = 0.0;
    for(int j=ia[row];j<ia[row+1];j++){
      dot+=a[j]*x[ja[j]];
    }
    q[row] = dot;
  }
}

void DotProduct(int n, double* xt, double* y,double* q){
  *q = 0.0;
  for(int i=0; i<n;i++){(*q)+=xt[i]*y[i];}
}

void Copy(int n, double* x, double* y){
  for(int i=0; i<n;i++){y[i]=x[i];}
}

void MaxValue(int n, double* x, double* maxval){
  *maxval=fabs(x[0]);
  double temp;
  for(int i=1;i<n;i++){
    temp = fabs(x[i]);
    if(temp>(*maxval)){
      *maxval=temp;
    }
  }
}

/* 
y = sign1*alpha1*x+sign2*alpha2*y 
*/
void Axpy(int n, double alpha1, double alpha2, int sign1, int sign2, double* x, double* y){
  double a1, a2;
  a1 = alpha1*sign1;
  a2 = alpha2*sign2;
  for(int i=0;i<n;i++){y[i]=a1*x[i]+a2*y[i];}
}
