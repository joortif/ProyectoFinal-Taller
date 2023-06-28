#ifndef PROYECTO_MEASUPDATE_H
#define PROYECTO_MEASUPDATE_H

#include "Matrix.h"

Matrix MeasUpdate(Matrix& x, Matrix& z, Matrix& g, Matrix& s, Matrix& G, Matrix& P, int n);

#endif //PROYECTO_MEASUPDATE_H
