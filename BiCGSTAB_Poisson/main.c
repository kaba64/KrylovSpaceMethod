#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparce_matrix.h"
#include "initialize.h"
#include "MVOperation.h"

const int im=100;
const int jm=200;
const double hx = 1.;
const double hy = 2.;
const double alpha = 1.0;
const double beta = 1.0;
const double eps = 1.0e-7;

int main(int argc, char *argv[])
{
  int *ia, *ja;
  double *a, *rh, *q1, *q2;
  double *x, *y, *dx, *dy;
  double *pn, *rn, *r0, *sn;
  double *x_sn;
  double alphan, gaman, betan, product_co;
  double maxresidual, residual, residual0, temp_dot_1, temp_dot_2, temp_dot_3;
  int iteration;
  int id, nnz, nrows;
  double elements[5], boundary_values[4];
  char name;
  double t1,t2;
  
  if(argc!=2){
    printf("Please enter the name of the iterative method.\n");
  }else{
    name = *argv[1];
  }

  id = im;
  nnz = 2*(4*(im-2)+4*(jm-2)+6)+5*(im-2)*(jm-2);
  nrows = im*jm;
  
  ia     = (int*)malloc((nrows+1)*sizeof(int));
  ja     = (int*)malloc(nnz*sizeof(int));
  a      = (double*)malloc(nnz*sizeof(double));
  x      = (double*)malloc(im*sizeof(double));
  y      = (double*)malloc(jm*sizeof(double));
  dx     = (double*)malloc(im*sizeof(double));
  dy     = (double*)malloc(jm*sizeof(double));
  
  x_sn   = (double*)malloc(nrows*sizeof(double));
  rh     = (double*)malloc(nrows*sizeof(double));
  rn     = (double*)malloc(nrows*sizeof(double));
  
  grid_generation(im,jm,hx,hy,x,y,dx,dy);
  elements[0] = -1.0/(dy[0]*dy[0]);
  elements[1] = -1.0/(dx[0]*dx[0]);
  elements[2] = 2.0/(dx[0]*dx[0])+2.0/(dy[0]*dy[0]);
  elements[3] = -1.0/(dx[0]*dx[0]);
  elements[4] = -1.0/(dy[0]*dy[0]);
  boundary_values[0]= 100.0;
  boundary_values[1]= 0.0;
  boundary_values[2]= 0.0;
  boundary_values[3]= 0.0;
  fill_rh(im,jm,id,alpha,beta,x,y,rh);
  sparce_matrix(im,jm,id,elements,'D',ia,ja,a);
  Updating_rh(im,jm,id,elements,boundary_values,'D',rh);
  /*
    Starting the computation for x
  */
  t1=spendtimeinsecond();
  iteration = 0;
  initialize_vector(im,jm,id,0.0,x_sn);
  Copy(nrows,rh,rn);
  DotProduct(nrows,rn,rn,&residual);
  residual = sqrt(residual);
  residual0 = residual;

  if(name=='G'){
    writing_residual('G',residual/residual0,iteration);
    while(residual>eps){
      for(int row=0;row<nrows;row++){
	double dot = rh[row];
	int colmun;
	for(int j=ia[row];j<ia[row+1];j++){
	  if(ja[j]==row){
	    colmun = j;
	  }else{
	    dot-=a[j]*x_sn[ja[j]];
	  }
	}
	x_sn[row] = dot/a[colmun];
      }
      /*Compute the residual*/
      for(int row=0;row<nrows;row++){
        double dot = rh[row];
        for(int j=ia[row];j<ia[row+1];j++){
	  dot-=a[j]*x_sn[ja[j]];
        }
        rn[row] = dot;
      }
      DotProduct(nrows,rn,rn,&residual);
      residual = sqrt(residual);
      iteration++;
      writing_residual('G',residual/residual0,iteration);
    }
    t2=spendtimeinsecond();
    writing_to_file(im,jm,id,'G',x,y,x_sn);
    printf("Time spent for computing x GSM %lf in second\n",t2-t1);
  }else if(name=='C'){
    pn     = (double*)malloc(nrows*sizeof(double));
    q1     = (double*)malloc(nrows*sizeof(double));
    q2     = (double*)malloc(nrows*sizeof(double));
    /*Start CGM*/
    writing_residual('C',residual/residual0,iteration);
    while(residual>eps){
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
      DotProduct(nrows,rn,rn,&residual);
      residual = sqrt(residual);
      iteration++;
      writing_residual('C',residual/residual0,iteration);
    }
    t2=spendtimeinsecond();
    printf("Time spent for computing x CCGM %lf in second\n",t2-t1);
    writing_to_file(im,jm,id,'C',x,y,x_sn);
    free(q1); free(q2); free(pn);
  }else if(name=='B'){
    /*
      Start BICGSTAB
    */
    r0     = (double*)malloc(nrows*sizeof(double));
    pn     = (double*)malloc(nrows*sizeof(double));
    q1     = (double*)malloc(nrows*sizeof(double));
    q2     = (double*)malloc(nrows*sizeof(double));
    sn     = (double*)malloc(nrows*sizeof(double));
    Copy(nrows,rn,r0);
    while(residual>eps){
      DotProduct(nrows,rn,r0,&temp_dot_2);
      if(temp_dot_2==0) {
	break;
      }
      if(iteration==0){
        Copy(nrows,rn,pn);
      }else{
	betan = (temp_dot_2/temp_dot_1)*(alphan/gaman);
	product_co = betan*gaman;
        Axpy4(nrows,1.0,product_co,betan,1,-1,1,rn,q1,pn);
      }
      MVmulti(nrows,ia,ja,a,pn,q1);
      DotProduct(nrows,r0,q1,&temp_dot_3);
      alphan = temp_dot_2/temp_dot_3;
      Axpy3(nrows,1.0,alphan,1,-1,rn,q1,sn);
      MVmulti(nrows,ia,ja,a,sn,q2);
      temp_dot_1=temp_dot_2;
      DotProduct(nrows,q2,sn,&temp_dot_2);
      DotProduct(nrows,q2,q2,&temp_dot_3);
      gaman = temp_dot_2/temp_dot_3;
      Axpy4(nrows,alphan,gaman,1.0,1,1,1,pn,sn,x_sn);
      Axpy3(nrows,1.0,gaman,1,-1,sn,q2,rn);
      DotProduct(nrows,rn,rn,&residual);
      residual = sqrt(residual);
      iteration++;
      writing_residual('B',residual/residual0,iteration);
    }
    t2=spendtimeinsecond();
    printf("Time spent for computing x BICGSTAB %lf in second\n",t2-t1);
    writing_to_file(im,jm,id,'B',x,y,x_sn);
    free(q1); free(q2); free(pn);
    free(r0); free(sn);
  }
  free(a); free(ia); free(ja); free(rh);
  free(x); free(y),free(dx); free(dy);
  free(x_sn); free(rn);
  
  return 0;
}
