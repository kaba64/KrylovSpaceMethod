#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparce_matrix.h"
#include "initialize.h"
#include "MVOperation.h"
#include "parameter.h"

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
  double maxresidual, residual, residual0, temp_dot_1, temp_dot_2, temp_dot_3;
  int iteration, k, n_w;
  int id, nnz, nrows;
  char name;
  double t1,t2;
  
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
  k=0;
  n_w = 0;
  while((k*dt)<200.5){
    printf("Outer loop : %d\n",k);
    iteration = 0;
    Copy(nrows,x_snp,x_sn);
    fill_rh(im,jm,id,eps_ch,s_ch,dt,x_sn,rh);
    MVmulti(nrows,ia,ja,a,x_snp,rn);
    Axpy(nrows,1.0,1.0,1,-1,rh,rn);
    DotProduct(nrows,rn,rn,&residual);
    residual = sqrt(residual);
    residual0 = residual;
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
      Axpy4(nrows,alphan,gaman,1.0,1,1,1,pn,sn,x_snp);
      Axpy3(nrows,1.0,gaman,1,-1,sn,q2,rn);
      DotProduct(nrows,rn,rn,&residual);
      residual = sqrt(residual);
      iteration++;
      writing_residual('B',residual/residual0,iteration);
    }
    if((k*dt)>=0.25 && n_w ==0){
      writing_to_file(im,jm,id,'a',x,y,x_snp);
      n_w++;
    }else if((k*dt)>=2.0 && n_w ==1){
      writing_to_file(im,jm,id,'b',x,y,x_snp);
      n_w++;
    }else if((k*dt)>=8.0 && n_w ==2){
      writing_to_file(im,jm,id,'c',x,y,x_snp);
      n_w++;
    }else if((k*dt)>=20.0 && n_w ==3){
      writing_to_file(im,jm,id,'d',x,y,x_snp);
      n_w++;
    }else if((k*dt)>=40.0 && n_w ==4){
      writing_to_file(im,jm,id,'e',x,y,x_snp);
      n_w++;
    }else if((k*dt)>=70.0 && n_w ==5){
      writing_to_file(im,jm,id,'f',x,y,x_snp);
      n_w++;
    }else if((k*dt)>=100.0 && n_w ==6){
      writing_to_file(im,jm,id,'g',x,y,x_snp);
      n_w++;
    }else if((k*dt)>=200.0 && n_w ==7){
      writing_to_file(im,jm,id,'h',x,y,x_snp);
      n_w++;
    }
    k++;
  }
  //t2=spendtimeinsecond();
  //printf("Time spent for computing x BICGSTAB %lf in second\n",t2-t1);
  writing_to_file(im,jm,id,'B',x,y,x_snp);
  
  free(q1); free(q2); free(pn);
  free(r0); free(sn);
  free(a); free(ia); free(ja); free(rh);
  free(x); free(y); free(dxsm); free(dxsp); free(dysm); free(dysp);
  free(dxsmi); free(dxspi); free(dysmi); free(dyspi);
  free(x_sn); free(x_snp); free(rn);
  
  return 0;
}
