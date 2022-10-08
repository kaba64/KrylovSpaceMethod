void grid_generation(int im, int jm, double hx, double hy, double* x, double* y, double* dx, double* dy);
void initialize_vector(int im, int jm, int id, double value, double * f);
void fill_rh(int im, int jm, int id, double alpha, double beta, double* x, double* y, double* f);
void writing_to_file(int im, int jm, int id, double* x, double* y, double* f);
void writing_to_file_exact(int im, int jm, int id, double* x, double* y, double hx, double hy, double t);
void writing_residual(char name, double res, int num);
