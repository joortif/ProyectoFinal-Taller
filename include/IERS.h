#ifndef PROYECTO_IERS_H
#define PROYECTO_IERS_H

#include "Matrix.h"

void IERS(double **eop, double Mjd_TT, double &UTC1_UTC, double &TAI_UTC, double &x_pole, double &y_pole);

#endif //PROYECTO_IERS_H
