void sparce_matrix(int im, int jm, int id, double* elements , char boundary, int* ia, int* ja, double* a){
  int row, num, sign;
  int i1, i2, i3, i4, i5;
  double a1, a2, a3, a4, a5;
  if(boundary=='D'){
    sign = -1;
  }else if(boundary=='N'){
    sign = 1;
  }
  a1 = elements[0];
  a2 = elements[1];
  a3 = elements[2];
  a4 = elements[3];
  a5 = elements[4];
  row = 0;
  num = 0;
  for(int j=0; j<jm; j++) {
    if(j==0){
      for (int i=0; i<im; i++) {
  	i2 = i-1;
  	i3 = i;
  	i4 = i+1;
  	i5 = i+id;
  	ia[row]=num;
  	if(i==0){
	  ja[num]= i3;
  	  a[num] = sign*(a1+a2)+a3;
  	  num++;
  	  ja[num]= i4;
	  a[num] = a4;
	  num++;
	  ja[num]= i5;
	  a[num] = a5;
	  num++;
  	}else if(i>0 && i<im-1){
  	  ja[num]= i2;
	  a[num] = a2;
  	  num++;
  	  ja[num]= i3;
	  a[num] = sign*a1+a3;
	  num++;
	  ja[num]= i4;
	  a[num] = a4;
	  num++;
	  ja[num]= i5;
	  a[num] = a5;
  	  num++;
  	}else{
	  ja[num]= i2;
	  a[num] = a2;
	  num++;
	  ja[num]= i3;
	  a[num] = sign*(a1+a3)+a4;
	  num++;
	  ja[num]= i5;
	  a[num] = a5;
	  num++;
  	}
	row++;
      }
    }else if(j>0 && j<jm-1){
      for (int i=0; i<im; i++) {
	i1 = (j-1)*id+i;
	i2 = j*id+i-1;
	i3 = j*id+i;
	i4 = j*id+i+1;
	i5 = (j+1)*id+i;
	ia[row]=num;
	if(i==0){
	  ja[num]= i1;
	  a[num] = a1;
	  num++;
	  ja[num]= i3;
	  a[num] = sign*a2+a3;
	  num++;
	  ja[num]= i4;
	  a[num] = a4;
	  num++;
	  ja[num]= i5;
	  a[num] = a5;
	  num++;
	}else if(i>0 && i<im-1){
	  ja[num]= i1;
	  a[num] = a1;
	  num++;
	  ja[num]= i2;
	  a[num] = a2;
	  num++;
	  ja[num]= i3;
	  a[num] = a3;
	  num++;
	  ja[num]= i4;
	  a[num] = a4;
	  num++;
	  ja[num]= i5;
	  a[num] = a5;
	  num++;
	}else{
	  ja[num]= i1;
	  a[num] = a1;
	  num++;
	  ja[num]= i2;
	  a[num] = a2;
	  num++;
	  ja[num]= i3;
	  a[num] = a3+a4*sign;
	  num++;
	  ja[num]= i5;
	  a[num] = a5;
	  num++;
	}
	row++;
      }
    }else{
      for (int i=0; i<im; i++) {
  	i1 = (j-1)*id+i;
  	i2 = j*id+i-1;
  	i3 = j*id+i;
  	i4 = j*id+i+1;
  	i5 = (j+1)*id+i;
  	ia[row]=num;
  	if(i==0){
  	  ja[num]= i1;
  	  a[num] = a1;
  	  num++;
  	  ja[num]= i3;
          a[num] = sign*(a2+a5)+a3;
  	  num++;
          ja[num]= i4;
	  a[num] = a4;
  	  num++;
  	}else if(i>0 && i<im-1){
	  ja[num]= i1;
	  a[num] = a1;
  	  num++;
          ja[num]= i2;
          a[num] = a2;
          num++;
          ja[num]= i3;
          a[num] = a3+a5*sign;
          num++;
          ja[num]= i4;
          a[num] = a4;
          num++;
  	}else{
          ja[num]= i1;
          a[num] = a1;
          num++;
          ja[num]= i2;
          a[num] = a2;
          num++;
          ja[num]= i3;
          a[num] = a3+(a4+a5)*sign;
          num++;
  	}
  	row++;
      }
    }
  }
  ia[row]=num;
}
