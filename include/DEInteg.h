#ifndef PROYECTO_DEINTEG_H
#define PROYECTO_DEINTEG_H

#include "Matrix.h"

Matrix DEInteg(Matrix (*func)(double, Matrix),double t,double tout,double relerr,double abserr,double n_eqn,Matrix y);
#endif //PROYECTO_DEINTEG_H
