#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "sparce_matrix.h"

const int im=4;
const int jm=4;
const double hx = 1.0;
const double hy = 1.0;

int main(int argc, char *argv[])
{
  int * ia;
  int * ja;
  double * a;
  int  num, id, nnz, nrows, row, column;
  int i1, i2, i3, i4, i5;
  double elements[5];
  
  elements[0] = -1.0;
  elements[1] = -1.0;
  elements[2] = 4.0;
  elements[3] = -1.0;
  elements[4] = -1.0;
  id = im;
  nnz = 2*(4*(im-2)+4*(jm-2)+6)+5*(im-2)*(jm-2);
  nrows = jm*jm+1;
  ia = (int*)malloc(nrows*sizeof(int));
  ja = (int*)malloc(nnz*sizeof(int));
  a  = (double*)malloc(nnz*sizeof(double));
  sparce_matrix(im,jm,id,elements,ia,ja,a);
  for(int i =0; i<nrows;i++){
    printf("ia(%d) = %d\n",i,ia[i]);
  }
  free(a);
  free(ia);
  free(ja);
  return 0;
}
