#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

const double pi = acos(-1.);
void grid_generation(int im, int jm, double hx, double hy, double* x, double* y,double* dxsm, double* dxsp, double* dysm, double* dysp,
		     double* dxsmi, double* dxspi, double* dysmi, double* dyspi){
  double dx_t, dy_t;
  dx_t = hx/im;
  dy_t = hy/jm;
  
  /* grids in x-direction */
  for(int i=0;i<im;i++){
    if(i==0){
      x[i]      = dx_t*0.5;
    }else{
      x[i]      = x[i-1]+dx_t;
    }
  }
  for(int i=0;i<im;i++){
    if(i==0){
      dxsm[i]   = x[i]*2;
      dxsp[i]   = x[i+1]-x[i];
    }else if(i>0 && i<im-1){
      dxsm[i]   = x[i]-x[i-1];      
      dxsp[i]   = x[i+1]-x[i];
    }else{
      dxsm[i]   = x[i]-x[i-1];
      dxsp[i]   = dxsm[i];
    }
    dxsmi[i] = 1.0/dxsm[i];
    dxspi[i] = 1.0/dxsp[i];
  }
  /* grids in y-direction */
  for(int j=0;j<jm;j++){
    if(j==0){
      y[j]     = dy_t*0.5;
    }else{
      y[j]     = y[j-1]+dy_t;
    }
   }
  for(int j=0;j<jm;j++){
    if(j==0){
      dysm[j]   = y[j]*2;
      dysp[j]   = y[j+1]-y[j];
    }else if(j>0 && j<jm-1){
      dysm[j]   = y[j]-y[j-1];
      dysp[j]   = y[j+1]-y[j];
    }else{
      dysm[j]   = y[j]-y[j-1];
      dysp[j]   = dysm[j];
    }
    dysmi[j] = 1.0/dysm[j];
    dyspi[j] = 1.0/dysp[j];
   }
}

void initialize_vector(int im, int jm, int id, double value, double * f){
  int jid, ig, igp;
  for(int j=0;j<jm;j++){
    jid = j*id;
    for(int i=0;i<im;i++){
      ig = 2*(jid+i);
      igp = ig+1;
      f[ig]=value;
      f[igp]=value;
    }
  }
}

void initialize_VC(int n, double value, double * f){
  for(int i=0;i<n;i++){
    f[i]=value;
  }
}

void initialize_phi(int im, int jm, int id, double hx, double hy, double* x, double* y, double radius, double eps_ch,double * f){
  double cx0, cx1, cx2, cy1, cy2;
  double distance;
  int jid, ig, igp;
  
  cx0 =	hx*0.5;
  cx1 = hx*0.5-radius;
  cy1 =	hy*0.5;
  cx2 =	hx*0.5+radius;
  cy2 = hy*0.5;
  
  for(int j=0;j<jm;j++){
    jid = j*id;
    for(int i=0;i<im;i++){
      ig = 2*(jid+i);
      igp = ig+1;
      f[ig] = 0.;
      if(x[i]<=cx0){
	distance = sqrt((x[i]-cx1)*(x[i]-cx1)+(y[j]-cy1)*(y[j]-cy1));
	f[igp]=tanh((radius-distance)/(eps_ch*sqrt(2.0)));
      }else{
	distance = sqrt((x[i]-cx2)*(x[i]-cx2)+(y[j]-cy2)*(y[j]-cy2));
	f[igp]=tanh((radius-distance)/(eps_ch*sqrt(2.0)));
      }
    }
  }
}

void fill_rh(int im, int jm, int id, double eps_ch, double s_ch, double dt, double* x_s, double* f){
  int idj, ig, igp;
  double epsdi, dti;
  
  epsdi = 1.0/(eps_ch*eps_ch);
  dti = 1.0/dt;
  
  for(int j=0;j<jm;j++){
    idj = id*j;
    for(int i=0;i<im;i++){
      ig     = 2*(idj+i);
      igp    = ig+1;
      f[ig]  = epsdi*x_s[igp]*(1.0-x_s[igp]*x_s[igp])+s_ch*epsdi*x_s[igp];
      f[igp] = dti*x_s[igp];
    }
  }
}

void writing_to_file(int im, int jm, int id, char name, double* x, double* y, double* f){
  
  FILE *ofile;
  char filename[] = "solution.txt";
  char name1[3] = {name ,'_','\0'};
  int jid, ng;
  
  ofile = fopen(strcat(name1,filename), "w");
  for(int j=0;j<jm;j++){
    jid = j*id;
    for(int i=0;i<im;i++){
      ng = 2*(jid+i)+1;
      fprintf(ofile, "%d\t%d\t%lf\t%lf\t%lf\n",i,j,x[i],y[j],f[ng]);
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

double spendtimeinsecond()
{
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
