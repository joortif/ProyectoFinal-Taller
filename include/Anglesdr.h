#ifndef PROYECTO_ANGLESDR_H
#define PROYECTO_ANGLESDR_H

#include "Matrix.h"

void anglesdr(double az1, double az2, double az3, double el1,double el2,double el3,double Mjd1,double Mjd2,double Mjd3,
              Matrix rsite1,Matrix rsite2,Matrix rsite3,Matrix & r2, Matrix & v2);

#endif //PROYECTO_ANGLESDR_H
