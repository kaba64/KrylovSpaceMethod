#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "sparce_matrix.h"
#include "initialize.h"
#include "MVOperation.h"

const int im=200;
const int jm=200;
const double hx = 2.;
const double hy = 2.;
const double alpha = 1.0;
const double beta = 1.0;
const double eps = 1.0e-12;

int main(int argc, char *argv[])
{
  int *ia, *ja;
  double *a, *rh, *q1, *q2;
  double *x, *y, *dx, *dy;
  double *pn, *pnm, *rn, *rnm;
  double *x_sn, *x_snm;
  double alphan, gaman;
  double maxresidual, temp_dot_1, temp_dot_2;
  int iteration;
  int id, nnz, nrows;
  double elements[5];
  
  id = im;
  nnz = 2*(4*(im-2)+4*(jm-2)+6)+5*(im-2)*(jm-2);
  nrows = jm*jm;
  
  ia     = (int*)malloc((nrows+1)*sizeof(int));
  ja     = (int*)malloc(nnz*sizeof(int));
  a      = (double*)malloc(nnz*sizeof(double));
  x      = (double*)malloc(im*sizeof(double));
  y      = (double*)malloc(jm*sizeof(double));
  dx     = (double*)malloc(im*sizeof(double));
  dy     = (double*)malloc(jm*sizeof(double));
  
  x_sn   = (double*)malloc(nrows*sizeof(double));
  rh     = (double*)malloc(nrows*sizeof(double));
  pn     = (double*)malloc(nrows*sizeof(double));
  rn     = (double*)malloc(nrows*sizeof(double));
  q1     = (double*)malloc(nrows*sizeof(double));
  q2     = (double*)malloc(nrows*sizeof(double));

  grid_generation(im,jm,hx,hy,x,y,dx,dy);
  elements[0] = -1.0/(dy[0]*dy[0]);
  elements[1] = -1.0/(dx[0]*dx[0]);
  elements[2] = 2.0/(dx[0]*dx[0])+2.0/(dy[0]*dy[0]);
  elements[3] = -1.0/(dx[0]*dx[0]);
  elements[4] = -1.0/(dy[0]*dy[0]);
  sparce_matrix(im,jm,id,elements,'D',ia,ja,a);
  fill_rh(im,jm,id,alpha,beta,x,y,rh);
  //writing_to_file(im,jm,id,x,y,rh);
  
  /*Start CGM*/
  iteration = 0;
  initialize_vector(im,jm,id,0.0,x_sn);
  Copy(nrows,rh,rn);
  MaxValue(nrows,rn,&maxresidual);
  writing_residual('C',maxresidual,iteration);
  while(maxresidual>eps){
    if(iteration==0){
      Copy(nrows,rn,pn);
    }else{
      MVmulti(nrows,ia,ja,a,rn,q1);
      MVmulti(nrows,ia,ja,a,pn,q2);
      DotProduct(nrows,pn,q1,&temp_dot_1);
      DotProduct(nrows,pn,q2,&temp_dot_2);
      gaman = -1*(temp_dot_1/temp_dot_2);
      Axpy(nrows,1.0,gaman,1,1,rn,pn);
    }
    DotProduct(nrows,rn,rn,&temp_dot_1);
    MVmulti(nrows,ia,ja,a,pn,q2);
    DotProduct(nrows,pn,q2,&temp_dot_2);
    alphan = temp_dot_1/temp_dot_2;
    Axpy(nrows,alphan,1.0,1,1,pn,x_sn);
    Axpy(nrows,alphan,1.0,-1,1,q2,rn);
    MaxValue(nrows,rn,&maxresidual);
    iteration++;
    writing_residual('C',maxresidual,iteration);
  }
  writing_to_file(im,jm,id,x,y,x_sn);
  free(a); free(ia); free(ja); free(rh);
  free(x); free(y),free(dx); free(dy);
  free(x_sn); free(q1);free(q2);
  free(pn); free(rn);
  
  return 0;
}
