#include <stdio.h>
#include <math.h>
#include <string.h>

const double pi = acos(-1.);

void grid_generation(int im, int jm, double hx, double hy, double* x, double* y, double* dx, double* dy){
  double dx_t, dy_t;
  dx_t = hx/im;
  dy_t = hy/jm;
  for(int i=0;i<im;i++){
    if(i==0){
      x[i]  = dx_t*0.5;
      dx[i] = dx_t;
    }else{
      x[i]  = x[i-1]+dx_t;
      dx[i] = dx_t;
    }
  }
  for(int j=0;j<jm;j++){
    if(j==0){
      y[j]  = dy_t*0.5;
      dy[j] = dy_t;
    }else{
      y[j]  = y[j-1]+dy_t;
      dy[j] = dy_t;
    }
   }
}

void initialize_vector(int im, int jm, int id, double value, double * f){
  for(int j=0;j<jm;j++){
    for(int i=0;i<im;i++){
      f[j*id+i]=value;
    }
  }
}

void fill_rh(int im, int jm, int id, double alpha, double beta, double* x, double* y, double* f){
  for(int j=0;j<jm;j++){
    for(int i=0;i<im;i++){
      f[j*id+i]= 0.0;
    }
  }
}

void writing_to_file(int im, int jm, int id, double* x, double* y, double* f){
  
  FILE *ofile;
  char filename[] = "solution.txt";
  ofile = fopen(filename, "w");
  for(int j=0;j<jm;j++){
    for(int i=0;i<im;i++){
      fprintf(ofile, "%d\t%d\t%lf\t%lf\t%lf\n",i,j,x[i],y[j],f[j*id+i]);
    }
  }
  fclose(ofile);
}

void writing_to_file_exact(int im, int jm, int id, double* x, double* y, double hx, double hy, double t){
  
  FILE *ofile;
  char filename[] = "exaxt_solution.txt";
  double sum;
  ofile = fopen(filename, "w");
  for(int j=0;j<jm;j++){
    for(int i=0;i<im;i++){
      sum=0.0;
      for(int k=1;k<100;k++){
	sum+=((1-((int)pow(-1,k)))*sinh((k*pi*(hy-y[j]))/hx)*sin((k*pi*x[i])/hx))/(k*sinh((k*pi*(hy))/hx));
      }
      sum = ((2*t)/pi)*sum;
      fprintf(ofile, "%d\t%d\t%lf\t%lf\t%lf\n",i,j,x[i],y[j],sum);
    }
  }
  fclose(ofile);
}

void writing_residual(char name, double res, int num){

  FILE *ofile;
  char filename[] = "residual.txt";
  char name1[3] = {name ,'_','\0'};
  ofile = fopen(strcat(name1,filename), "a");
  fprintf(ofile, "%d\t%.15f\n",num,res);
  fclose(ofile);
}
