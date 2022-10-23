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
const double eps = 1.0e-6;
const int max_m = 20;

int main(int argc, char *argv[])
{
  int *ia, *ja;
  double *a, *rh, *q1;
  double *x, *y, *dx, *dy;
  double *pn, *rn, *r0, *sn;
  double *x_sn;
  double alphan, gaman, betan, product_co;
  double maxresidual, residual, residual0, temp_dot_1, temp_dot_2, temp_dot_3;
  int iteration;
  int id, nnz, nrows;
  double elements[5], boundary_values[4];
  char name;
  double t1,t2, d_r, temp_h;
  double *v, *rotation, *gn, *ym;
  int position_v_in, position_v_out, position_h, position_h_dot, position_r, count;
  double h[max_m+1][max_m];
    
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
  
  x_sn     = (double*)malloc(nrows*sizeof(double));
  rh       = (double*)malloc(nrows*sizeof(double));
  q1       = (double*)malloc(nrows*sizeof(double));
  v        = (double*)malloc((nrows*(max_m+1))*sizeof(double));
  rotation = (double*)malloc(((2*max_m)*sizeof(double)));
  gn       = (double*)malloc((max_m+1)*sizeof(double));
  ym       = (double*)malloc(max_m*sizeof(double));
  
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
  initialize_h((max_m+1),max_m,0.0,h);
  t1=spendtimeinsecond();
  initialize_vector(im,jm,id,0.0,x_sn);
  Copy(nrows,rh,q1);
  DotProduct(nrows,q1,q1,&residual);
  residual = sqrt(residual);
  residual0 = residual;
  writing_residual('G',residual/residual0,0);
  count=1;
  while(residual>eps){
    
    gn[0] = residual;
    Scale_V(nrows,residual,'Y',q1,&v[0]);
    
    for(iteration=0;iteration<max_m;iteration++){
      /*Compute Avj*/
      position_v_out = nrows*iteration;
      MVmulti(nrows,ia,ja,a,&v[position_v_out],q1);
      position_h=0;
      for(int i=0;i<=iteration;i++){
	/*Compute hij*/
	position_v_in = nrows*i;
	DotProduct(nrows,q1,&v[position_v_in],&temp_dot_1);
	
	h[position_h][iteration] = temp_dot_1;
	position_h++;
	
	/*Compute wj-hijvi*/
	Axpy(nrows,temp_dot_1,1.0,-1,1,&v[position_v_in],q1);
      }
      
      DotProduct(nrows,q1,q1,&temp_dot_1);
      temp_dot_1 = sqrt(temp_dot_1);
   
      h[position_h][iteration]=temp_dot_1;
      Scale_V(nrows,temp_dot_1,'Y',q1,&v[nrows*(iteration+1)]);    
   
      position_h = 0;
      for(int i=0;i<iteration;i++){
	position_r      = i*2;
	temp_h          = rotation[position_r]*h[position_h][iteration]+rotation[position_r+1]*h[position_h+1][iteration];
        h[position_h+1][iteration] = -rotation[position_r+1]*h[position_h][iteration]+rotation[position_r]*h[position_h+1][iteration];
        h[position_h][iteration]   = temp_h;
	position_h++;
      }
      
      position_h = iteration;
      d_r                    = sqrt(h[position_h][iteration]*h[position_h][iteration]+h[position_h+1][iteration]*h[position_h+1][iteration]);
      position_r             = iteration*2;
      rotation[position_r]   = h[position_h][iteration]/d_r;
      rotation[position_r+1] = h[position_h+1][iteration]/d_r;
      h[position_h][iteration]          = d_r;
      h[position_h+1][iteration]        = 0.0;
      gn[iteration+1]        = -rotation[position_r+1]*gn[iteration];
      gn[iteration]          = rotation[position_r]*gn[iteration];
      count++; 
      writing_residual('G',fabs(gn[iteration+1])/residual0,count);
    }
    
    for(int i=max_m-1;i>=0;i--){   
      double dot1=0.0;
      for(int j=i+1;j<max_m;j++){
      	dot1-=(ym[j]*h[i][j]);
      }
      dot1+=gn[i];
      ym[i] = dot1/h[i][i];
      position_v_out = nrows*i;
      Axpy(nrows,ym[i],1.0,1,1,&v[position_v_out],x_sn);
    }
    MVmulti(nrows,ia,ja,a,x_sn,q1);
    Axpy(nrows,1.0,1.0,1,-1,rh,q1);
    DotProduct(nrows,q1,q1,&residual);
    residual = sqrt(residual);
    //printf("end error : %lf\n",residual);
  }
  t2=spendtimeinsecond();
  printf("Time spent for computing x GMRES %lf in second\n",t2-t1);
  writing_to_file(im,jm,id,'G',x,y,x_sn);
  free(a); free(ia); free(ja); free(rh);
  free(x); free(y),free(dx); free(dy);
  free(x_sn);
  free(v); free(rotation); free(gn);
  free(q1);
  return 0;
}
