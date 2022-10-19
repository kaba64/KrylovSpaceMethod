void grid_generation(int im, int jm, double hx, double hy, double* x, double* y,double* dxsm, double* dxsp, double* dysm, double* dysp,
                     double* dxsmi, double* dxspi, double* dysmi, double* dyspi);
void initialize_vector(int im, int jm, int id, double value, double * f);
void initialize_VC(int i, double value, double * f);
void initialize_phi(int im, int jm, int id, double hx, double hy, double* x, double* y, double radius, double eps_ch,double * f);
void fill_rh(int im, int jm, int id, double eps_ch, double s_ch, double dt, double* x_s, double* f);
void writing_to_file(int im, int jm, int id, char name, double* x, double* y, double* f);
void writing_to_file_exact(int im, int jm, int id, double* x, double* y, double hx, double hy, double t);
void writing_residual(char name, double res, int num);
double spendtimeinsecond();
