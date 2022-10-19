#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparce_matrix.h"
#include "initialize.h"
#include "MVOperation.h"
#include "parameter.h"

const int max_m = 4;

int main(int argc, char *argv[])
{
  int *ia, *ja;
  double *a, *rh, *q1, *rn;
  double *x, *y;
  double *dxsp, *dxsm, *dysp, *dysm;
  double *dxspi, *dxsmi, *dyspi, *dysmi;
  double *v, *h, *rotation, *gn, *ym;
  double *x_sn, *x_snp;
  double maxresidual, residual, residual0, temp_dot_1, temp_dot_2;
  double d_r, temp_h;
  int iteration, position_v_in, position_v_out, position_h, position_r;
  int id, nnz, nrows;
  char name;
  double t1,t2;
  
  id       = im;
  nnz      = 4*(5*(im-2)+5*(jm-2)+8)+12*(im-2)*(jm-2);
  nrows    = 2*im*jm;
  
  ia       = (int*)malloc((nrows+1)*sizeof(int));
  ja       = (int*)malloc(nnz*sizeof(int));
  a        = (double*)malloc(nnz*sizeof(double));
  dxsp     = (double*)malloc(im*sizeof(double));
  dxsm     = (double*)malloc(im*sizeof(double));
  dysp     = (double*)malloc(jm*sizeof(double));
  dysm     = (double*)malloc(jm*sizeof(double));
  dxspi    = (double*)malloc(im*sizeof(double));
  dxsmi    = (double*)malloc(im*sizeof(double));
  dyspi    = (double*)malloc(jm*sizeof(double));
  dysmi    = (double*)malloc(jm*sizeof(double));
  x        = (double*)malloc(im*sizeof(double));
  y        = (double*)malloc(jm*sizeof(double));
  
  x_sn     = (double*)malloc(nrows*sizeof(double));
  x_snp    = (double*)malloc(nrows*sizeof(double));
  rh       = (double*)malloc(nrows*sizeof(double));
  rn       = (double*)malloc(nrows*sizeof(double));
  q1       = (double*)malloc(nrows*sizeof(double));
  v        = (double*)malloc((nrows*(max_m+1))*sizeof(double));
  h        = (double*)malloc(((max_m+1)*max_m)*sizeof(double));
  rotation = (double*)malloc(((2*max_m)*sizeof(double)));
  gn       = (double*)malloc((max_m+1)*sizeof(double));
  ym       = (double*)malloc(max_m*sizeof(double));
  /*Generate the grid*/
  grid_generation(im,jm,hx,hy,x,y,dxsm,dxsp,dysm,dysp,dxsmi,dxspi,dysmi,dyspi);
  sparce_matrix(im,jm,id,dxsm,dxsp,dysm,dysp,dxsmi,dxspi,dysmi,dyspi,dt,gama0,mobility,lamda_ch,s_ch,eps_ch,'N',ia,ja,a);
  /*Initialize the phi*/
  initialize_phi(im,jm,id,hx,hy,x,y,l_ref,eps_ch,x_sn);
  Copy(nrows,x_sn,x_snp);
  /*Fill the right-hand side*/
  fill_rh(im,jm,id,eps_ch,s_ch,dt,x_sn,rh);
  /*Computing the error*/
  MVmulti(nrows,ia,ja,a,x_snp,rn);
  Axpy(nrows,1.0,1.0,1,-1,rh,rn);
  
  DotProduct(nrows,rn,rn,&residual);
  residual = sqrt(residual);
  residual0 = residual;
  /*Computing v0*/
  //Copy(nrows,rn,&v[0]);
  Scale_V(nrows,residual0,'Y',rn,&v[0]);
  /*initialize g0 for computing error*/
  initialize_VC(max_m+1,0.0,gn);
  gn[0] = residual0;
  
  for(iteration=0;iteration<max_m;iteration++){
    /*Compute Avj*/
    position_v_out = nrows*iteration;
    MVmulti(nrows,ia,ja,a,&v[position_v_out],q1);
    
    for(int i=0;i<=iteration;i++){
      /*Compute hij*/
      position_v_in = nrows*i;
      DotProduct(nrows,q1,&v[position_v_in],&temp_dot_1);
      position_h = iteration*(max_m+1)+i;
      //printf("at %d : h(%d,%d) = %.16f \n",position_h,iteration,i,temp_dot_1);
      h[position_h]=temp_dot_1;
      //printf("h(%d) : %lf\n", position_h,temp_dot_1);
      /*Compute wj-hijvi*/
      Axpy(nrows,temp_dot_1,1.0,-1,1,&v[position_v_in],q1);
    }
    
    DotProduct(nrows,q1,q1,&temp_dot_1);
    temp_dot_1 = sqrt(temp_dot_1);
    position_h = iteration*(max_m+1)+iteration+1;
    h[position_h]=temp_dot_1;
    
    temp_dot_2 =h[position_h];
    
    position_h = iteration*(max_m+1);
    for(int i=0;i<iteration;i++){
      position_r      = i*2;
      h[position_h]   = rotation[position_r]*h[position_h]+rotation[position_r+1]*h[position_h+1];
      h[position_h+1] = -rotation[position_r+1]*h[position_h]+rotation[position_r]*h[position_h+1];
      gn[i+1]         = -rotation[position_r+1]*gn[i];
      gn[i]           = rotation[position_r]*gn[i];
      position_h++;
    }
    
    d_r                    = sqrt(h[position_h]*h[position_h]+h[position_h+1]*h[position_h+1]);
    position_r             = iteration*2;
    rotation[position_r]   = h[position_h]/d_r;
    rotation[position_r+1] = h[position_h+1]/d_r;
    temp_h                 = rotation[position_r]*h[position_h]+rotation[position_r+1]*h[position_h+1];
    h[position_h+1]        = -rotation[position_r+1]*h[position_h]+rotation[position_r]*h[position_h+1];
    h[position_h]          = temp_h;
    gn[iteration+1]        = -rotation[position_r+1]*gn[iteration];
    gn[iteration]          = rotation[position_r]*gn[iteration];
    printf("error = %lf\n",gn[iteration]);
    
    if(temp_dot_1<1.0e-16){
      break;
    }
    Scale_V(nrows,temp_dot_1,'Y',q1,&v[nrows*(iteration+1)]);
    //printf("h+1_{%d} = %.30f\n",iteration,temp_dot_1);
    //for(int i=0;i<=(iteration+1);i++){
    //  position_h = iteration*(max_m+1)+i;
    //  printf("h_{%d,%d} = %.16f\n",i,iteration,h[position_h]);
    // }
  }
  //for(int i=iteration;i>=0;i--){
  //  double dot=0.0;
  //  position_h = (iteration-i)*(max_m+1)+i-1;
  //  for(int j=i;j>=iteration;j--){
  //    position_h = ;
  //    dot=+(ym[j]*h[]);
  //  }
  //  ym[i] = ;
  //}
    
  free(a); free(ia); free(ja);
  free(x); free(y); free(dxsm); free(dxsp); free(dysm); free(dysp);
  free(dxsmi); free(dxspi); free(dysmi); free(dyspi);
  free(x_sn); free(x_snp); free(rn); free(rh); free(q1);
  free(v); free(h); free(rotation); free(gn);
  
  return 0;
}
