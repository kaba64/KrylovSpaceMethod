#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparce_matrix.h"
#include "initialize.h"
#include "MVOperation.h"
#include "parameter.h"

const int max_m = 10;

void print_c(int n, double* x){
  for (int i=0; i<n; i++)
    {
      printf("a(%d) = %lf \n",i,x[i]);
    }
}

int main(int argc, char *argv[])
{
  int *ia, *ja;
  double *a, *rh, *q1, *q2;
  double *x, *y;
  double *dxsp, *dxsm, *dysp, *dysm;
  double *dxspi, *dxsmi, *dyspi, *dysmi;
  double *pn, *rn, *r0, *sn;
  double *x_sn, *x_snp;
  double alphan, gaman, betan, product_co;
  double maxresidual, residual, residual0, temp_dot_1, temp_dot_2, temp_dot_3, scale_g;
  int iteration, k;
  int id, nnz, nrows;
  char name;
  double t1,t2;
  double v[max_m][2*im*jm], h[max_m][2*im*jm];
  
  id     = im;
  nnz    = 4*(5*(im-2)+5*(jm-2)+8)+12*(im-2)*(jm-2);
  nrows  = 2*im*jm;
  
  ia     = (int*)malloc((nrows+1)*sizeof(int));
  ja     = (int*)malloc(nnz*sizeof(int));
  a      = (double*)malloc(nnz*sizeof(double));
  dxsp   = (double*)malloc(im*sizeof(double));
  dxsm   = (double*)malloc(im*sizeof(double));
  dysp   = (double*)malloc(jm*sizeof(double));
  dysm   = (double*)malloc(jm*sizeof(double));
  dxspi  = (double*)malloc(im*sizeof(double));
  dxsmi  = (double*)malloc(im*sizeof(double));
  dyspi  = (double*)malloc(jm*sizeof(double));
  dysmi  = (double*)malloc(jm*sizeof(double));
  x      = (double*)malloc(im*sizeof(double));
  y      = (double*)malloc(jm*sizeof(double));
  
  x_sn   = (double*)malloc(nrows*sizeof(double));
  x_snp  = (double*)malloc(nrows*sizeof(double));
  rh     = (double*)malloc(nrows*sizeof(double));
  rn     = (double*)malloc(nrows*sizeof(double));
  
  grid_generation(im,jm,hx,hy,x,y,dxsm,dxsp,dysm,dysp,dxsmi,dxspi,dysmi,dyspi);
  sparce_matrix(im,jm,id,dxsm,dxsp,dysm,dysp,dxsmi,dxspi,dysmi,dyspi,dt,gama0,mobility,lamda_ch,s_ch,eps_ch,'N',ia,ja,a);
  initialize_phi(im,jm,id,hx,hy,x,y,l_ref,eps_ch,x_sn);
  initialize_phi(im,jm,id,hx,hy,x,y,l_ref,eps_ch,x_snp);
  //fill_rh(im,jm,id,eps_ch,s_ch,dt,x_sn,rh);
  //Updating_rh(im,jm,id,elements,boundary_values,'D',rh);
  /*
    Starting the computation for x
  */
  //t1=spendtimeinsecond();
  //iteration = 0;
  //Copy(nrows,x_sn,x_snp);
  //MVmulti(nrows,ia,ja,a,x_snp,rn);
  //Axpy(nrows,1.0,1.0,1,-1,rh,rn);
  //DotProduct(nrows,rn,rn,&residual);
  //residual = sqrt(residual);
  //residual0 = residual;
  
  /*
    Start BICGSTAB
  */
  r0     = (double*)malloc(nrows*sizeof(double));
  pn     = (double*)malloc(nrows*sizeof(double));
  q1     = (double*)malloc(nrows*sizeof(double));
  q2     = (double*)malloc(nrows*sizeof(double));
  sn     = (double*)malloc(nrows*sizeof(double));
  //while((k*dt)<0.25){
  //  printf("Outer loop : %d\n",k);
  //iteration = 0;
  Copy(nrows,x_snp,x_sn);
  fill_rh(im,jm,id,eps_ch,s_ch,dt,x_sn,rh);
  MVmulti(nrows,ia,ja,a,x_snp,rn);
  Axpy(nrows,1.0,1.0,1,-1,rh,rn);
  DotProduct(nrows,rn,rn,&residual);
  residual = sqrt(residual);
  residual0 = residual;
  //Copy(nrows,rn,v[0]);
  scale_g = residual0;
  //Scale_V(nrows,scale_g,'Y',v[0]);
  //  while(residual>eps){
  //for(iteration=0;iteration<max_m;iteration++){
    //MVmulti(nrows,ia,ja,a,v[iteration],q1);
    //for(int i=0;i<=iteration;i++){
    //  DotProduct(nrows,q1,v[i],&temp_dot_1);
    //  h[iteration][i]=temp_dot_1;
    //}
    //Axpy3(nrows,1.0,h[iteration][0],1,-1,q1,v[0],v[iteration+1]);
    //for(int i=1;i<=iteration;i++){
    //  Axpy(nrows,h[iteration][i],1.0,-1,1,v[i],v[iteration+1]);
    //}
    //DotProduct(nrows,v[iteration+1],v[iteration+1],&temp_dot_1);
    //h[iteration][iteration+1]=temp_dot_1;
    //Scale_V(nrows,h[iteration][iteration+1],'Y',v[iteration+1]);
    
  //}
  
  //t2=spendtimeinsecond();
  //printf("Time spent for computing x BICGSTAB %lf in second\n",t2-t1);
  //writing_to_file(im,jm,id,'G',x,y,x_snp);
  
  free(q1); free(q2); free(pn);
  free(r0); free(sn);
  free(a); free(ia); free(ja); free(rh);
  free(x); free(y); free(dxsm); free(dxsp); free(dysm); free(dysp);
  free(dxsmi); free(dxspi); free(dysmi); free(dyspi);
  free(x_sn); free(x_snp); free(rn);
  
  return 0;
}
