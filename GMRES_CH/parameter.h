#include <math.h>

const double hx = 1.;
const double hy = 1.;
const double sigma =0.072;
const double viscosity = 0.001;
const double uref = sigma/viscosity;
const double Cn = 0.025;
const double l_ref = 0.19;
const double eps_ch = 0.01;//Cn*l_ref;
const double Pe = 100.;
const double mobility = 0.001;//(2.0*sqrt(2.0)*uref*l_ref*eps_ch)/(3.0*Pe*sigma);
const double lamda_ch = 0.001;//(3.0/(2.0*sqrt(2.0)))*sigma*eps_ch;
const double dt = 0.002;
const double s_ch = 1.2;
const double gama0 = 1.0;
const int im=200;
const int jm=200;
const double eps = 1.0e-10;
