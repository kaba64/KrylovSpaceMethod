void MVmulti(int nrows, int* ia, int* ja, double* a, double* x, double* q);
void DotProduct(int n, double* xt, double* y,double* q);
void Copy(int n, double* x, double* y);
void MaxValue(int n, double* x, double* maxval);
void Axpy(int n, double alpha1, double alpha2, int sign1, int sign2, double* x, double* y);
void Axpy3(int n, double alpha1, double alpha2, int sign1, int sign2, double* x, double* y, double* z);
void Axpy4(int n, double alpha1, double alpha2, double alpha3, int sign1, int sign2, int sign3, double* x, double* y, double* z);
void Scale_V(int nrows, double value, char inverse, double* x);
