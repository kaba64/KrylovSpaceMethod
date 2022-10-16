void sparce_matrix(int im, int jm, int id, double* dxsm, double* dxsp, double* dysm, double* dysp, double* dxsmi, double* dxspi,
		   double* dysmi, double* dyspi, double dt, double gama0, double mobility, double lamda_ch, double s_ch, double eps_ch,
		   char boundary, int* ia, int* ja, double* a){
  int row, num, sign;
  int ai0, ai1, ai2, ai3, ai4, ai5;
  double a0, a1, a2, a3, a4, a5;
  int bi0, bi1, bi2, bi3, bi4, bi5;
  double b0, b1, b2, b3, b4, b5;
  int is, js, jid;
  double value_b, value_l, value_r, value_t, eps_2i;
  
  if(boundary=='D'){
    sign = -1;
  }else if(boundary=='N'){
    sign = 1;
  }
  a0 = gama0/dt;
  b0 = -1.0;
  eps_2i = 1./(eps_ch*eps_ch);
  row = 0;
  num = 0;
  for(int j=0; j<jm; j++) {
    jid = j*id;
    b1 = (-2.0*dysmi[j])/(dysp[j]+dysm[j]);
    b5 = (-2.0*dyspi[j])/(dysp[j]+dysm[j]);
    a1 = mobility*lamda_ch*b1;
    a5 = mobility*lamda_ch*b5;
    if(j==0){
      for (int i=0; i<im; i++) {
	b2 = (-2.0*dxsmi[i])/(dxsp[i]+dxsm[i]);
	b3 = 2.0*(dxspi[i]*dxsmi[i]+dyspi[j]*dysmi[j]);
        b4 = (-2.0*dxspi[i])/(dxsp[i]+dxsm[i]);
	a2 = mobility*lamda_ch*b2;
        a3 = mobility*lamda_ch*b3;
        a4 = mobility*lamda_ch*b4;
	b3 +=(s_ch*eps_2i); 
	bi0 = 2*(jid+i);
	ai0 = bi0+1;
        ai1 = 2*(jid-id+i);
        ai2 = 2*(jid+i-1);
        ai3 = 2*(jid+i);
        ai4 = 2*(jid+i+1);
        ai5 = 2*(jid+id+i);
	bi1 = ai1+1;
	bi2 = ai2+1;
	bi3 = ai3+1;
	bi4 = ai4+1;
	bi5 = ai5+1;
	
  	ia[row]=num;
  	if(i==0){
	  ja[num]= bi0;
	  a[num] = b0;
	  num++;
	  ja[num]= bi3;
  	  a[num] = sign*(b1+b2)+b3;
  	  num++;
  	  ja[num]= bi4;
	  a[num] = b4;
	  num++;
	  ja[num]= bi5;
	  a[num] = b5;
	  num++;
	  /*A new row*/
	  row++;
	  ia[row]=num;
	  ja[num]= ai3;
          a[num] = sign*(a1+a2)+a3;
          num++;
	  ja[num]= ai0;
	  a[num] = a0;
	  num++;
	  ja[num]= ai4;
          a[num] = a4;
          num++;
          ja[num]= ai5;
          a[num] = a5;
          num++;
  	}else if(i>0 && i<im-1){
	  ja[num]= bi2;
          a[num] = b2;
          num++;
	  ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = sign*b1+b3;
          num++;
          ja[num]= bi4;
          a[num] = b4;
          num++;
          ja[num]= bi5;
          a[num] = b5;
          num++;
	  /*A new row*/
          row++;
	  ia[row]=num;
          ja[num]= ai2;
          a[num] = a2;
          num++;
          ja[num]= ai3;
          a[num] = sign*a1+a3;
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
          ja[num]= ai4;
          a[num] = a4;
          num++;
          ja[num]= ai5;
          a[num] = a5;
          num++;
  	}else{
	  ja[num]= bi2;
          a[num] = b2;
          num++;
          ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = sign*(b1+b4)+b3;
          num++;
          ja[num]= bi5;
          a[num] = b5;
          num++;
          /*A new row*/
	  row++;
          ia[row]=num;
          ja[num]= ai2;
          a[num] = a2;
          num++;
          ja[num]= ai3;
          a[num] = sign*(a1+a4)+a3;
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
          ja[num]= ai5;
          a[num] = a5;
          num++;
  	}
	row++;
      }
    }else if(j>0 && j<jm-1){
      for (int i=0; i<im; i++) {
	b2  = (-2.0*dxsmi[i])/(dxsp[i]+dxsm[i]);
        b3  = 2.0*(dxspi[i]*dxsmi[i]+dyspi[j]*dysmi[j]);
	b4  = (-2.0*dxspi[i])/(dxsp[i]+dxsm[i]);
        a2  = mobility*lamda_ch*b2;
        a3  = mobility*lamda_ch*b3;
        a4  = mobility*lamda_ch*b4;
        b3 += (s_ch*eps_2i);
        bi0 = 2*(jid+i);
        ai0 = bi0+1;
        ai1 = 2*(jid-id+i);
        ai2 = 2*(jid+i-1);
        ai3 = 2*(jid+i);
        ai4 = 2*(jid+i+1);
        ai5 = 2*(jid+id+i);
        bi1 = ai1+1;
        bi2 = ai2+1;
        bi3 = ai3+1;
        bi4 = ai4+1;
        bi5 = ai5+1;
	
	ia[row]=num;
	if(i==0){
	  ja[num]= bi1;
          a[num] = b1;
          num++;
          ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = sign*b2+b3;
          num++;
          ja[num]= bi4;
          a[num] = b4;
          num++;
          ja[num]= bi5;
          a[num] = b5;
          num++;
          /*A new row start*/
          row++;
          ia[row]=num;
	  /**/
          ja[num]= ai1;
          a[num] = a1;
          num++;
	  ja[num]= ai3;
          a[num] = sign*a2+a3;
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
          ja[num]= ai4;
          a[num] = a4;
          num++;
          ja[num]= ai5;
          a[num] = a5;
          num++;
	}else if(i>0 && i<im-1){
	  ja[num]= bi1;
          a[num] = b1;
          num++;
	  ja[num]= bi2;
          a[num] = b2;
          num++;
          ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = b3;
          num++;
          ja[num]= bi4;
          a[num] = b4;
          num++;
          ja[num]= bi5;
          a[num] = b5;
          num++;
          /*A new row start*/
	  row++;
          ia[row]=num;
	  /**/
	  ja[num]= ai1;
          a[num] = a1;
          num++;
          ja[num]= ai2;
          a[num] = a2;
          num++;
          ja[num]= ai3;
          a[num] = a3;
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
          ja[num]= ai4;
          a[num] = a4;
          num++;
          ja[num]= ai5;
          a[num] = a5;
          num++;
	}else{
	  ja[num]= bi1;
          a[num] = b1;
          num++;
          ja[num]= bi2;
          a[num] = b2;
          num++;
          ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = b3+sign*b4;
          num++;
          ja[num]= bi5;
          a[num] = b5;
          num++;
          /*A new row start*/
          row++;
          ia[row]=num;
	  /**/
          ja[num]= ai1;
          a[num] = a1;
          num++;
          ja[num]= ai2;
          a[num] = a2;
          num++;
          ja[num]= ai3;
          a[num] = a3+sign*a4;
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
          ja[num]= ai5;
          a[num] = a5;
          num++;
	}
	row++;
      }
    }else{
      for (int i=0; i<im; i++) {
	b2  = (-2.0*dxsmi[i])/(dxsp[i]+dxsm[i]);
	b3  = 2.0*(dxspi[i]*dxsmi[i]+dyspi[j]*dysmi[j]);
	b4  = (-2.0*dxspi[i])/(dxsp[i]+dxsm[i]);
        a2  = mobility*lamda_ch*b2;
        a3  = mobility*lamda_ch*b3;
        a4  = mobility*lamda_ch*b4;
	b3 += (s_ch*eps_2i);
        bi0 = 2*(jid+i);
        ai0 = bi0+1;
        ai1 = 2*(jid-id+i);
        ai2 = 2*(jid+i-1);
        ai3 = 2*(jid+i);
        ai4 = 2*(jid+i+1);
        ai5 = 2*(jid+id+i);
        bi1 = ai1+1;
        bi2 = ai2+1;
        bi3 = ai3+1;
        bi4 = ai4+1;
        bi5 = ai5+1;
	
  	ia[row]=num;
  	if(i==0){
  	  ja[num]= bi1;
          a[num] = b1;
          num++;
          ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = b3+sign*(b2+b5);
          num++;
          ja[num]= bi4;
          a[num] = b4;
          num++;
          /*A new row start*/
          row++;
          ia[row]=num;
          /**/
	  ja[num]= ai1;
          a[num] = a1;
          num++;
          ja[num]= ai3;
          a[num] = a3+sign*(a2+a5);
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
          ja[num]= ai4;
          a[num] = a4;
          num++;
  	}else if(i>0 && i<im-1){
	  ja[num]= bi1;
          a[num] = b1;
          num++;
          ja[num]= bi2;
          a[num] = b2;
          num++;
          ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = b3+sign*b5;
          num++;
          ja[num]= bi4;
          a[num] = b4;
          num++;
          /*A new row start*/
          row++;
          ia[row]=num;
          /**/
	  ja[num]= ai1;
          a[num] = a1;
          num++;
          ja[num]= ai2;
          a[num] = a2;
          num++;
          ja[num]= ai3;
          a[num] = a3+sign*a5;
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
          ja[num]= ai4;
          a[num] = a4;
          num++;
  	}else{
          ja[num]= bi1;
          a[num] = b1;
          num++;
          ja[num]= bi2;
          a[num] = b2;
          num++;
          ja[num]= bi0;
          a[num] = b0;
          num++;
          ja[num]= bi3;
          a[num] = b3+sign*(b4+b5);
          num++;
          /*A new row start*/
          row++;
          ia[row]=num;
          /**/
	  ja[num]= ai1;
          a[num] = a1;
          num++;
          ja[num]= ai2;
          a[num] = a2;
          num++;
          ja[num]= ai3;
          a[num] = a3+sign*(a4+a5);
          num++;
          ja[num]= ai0;
          a[num] = a0;
          num++;
  	}
  	row++;
      }
    }
  }
  ia[row]=num;
}
/*
  Updating the right hand side
*/
void Updating_rh(int im, int jm, int id, double* elements, double* values , char boundary, double* rh){
  double a1, a2, a3, a4, a5;
  double value_b, value_l, value_r, value_t;
  
  if(boundary=='D'){
    value_b = values[0];
    value_l = values[1];
    value_r = values[2];
    value_t = values[3];
  }else if(boundary=='N'){
    value_b = 0.0;
    value_l = 0.0;
    value_r = 0.0;
    value_t = 0.0;
  }
  a1 = elements[0];
  a2 = elements[1];
  a3 = elements[2];
  a4 = elements[3];
  a5 = elements[4];
  for(int j=0; j<jm; j++) {
    if(j==0){
      for (int i=0; i<im; i++) {
  	if(i==0){
	  rh[j*id+i] -= 2*(a1*value_b+a2*value_l); 
  	}else if(i>0 && i<im-1){
	  rh[j*id+i] -= 2*(a1*value_b);
  	}else{
	  rh[j*id+i] -= 2*(a1*value_b+a4*value_r);
  	}
      }
    }else if(j>0 && j<jm-1){
      rh[j*id+0] -= 2*(a2*value_l);
      rh[j*id+im-1] -= 2*(a4*value_r);
    }else{
      for (int i=0; i<im; i++) {
	if(i==0){
	  rh[j*id+i] -= 2*(a2*value_l+a5*value_t);
	}else if(i>0 && i<im-1){
	  rh[j*id+i] -= 2*(a5*value_t);
	}else{
	  rh[j*id+i] -= 2*(a4*value_r+a5*value_t);
	}
      }
    }
  }
}
