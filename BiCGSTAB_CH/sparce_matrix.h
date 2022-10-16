void sparce_matrix(int im, int jm, int id, double* dxsm, double* dxsp, double* dysm, double* dysp, double* dxsmi, double* dxspi,
                   double* dysmi, double* dyspi, double dt, double gama0, double mobility, double lamda_ch, double s_ch, double eps_ch,
                   char boundary, int* ia, int* ja, double* a);
void Updating_rh(int im, int jm, int id, double* elements, double* values , char boundary, double* rh);
