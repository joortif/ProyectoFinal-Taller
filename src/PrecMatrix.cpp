#include "../include/Matrix.h"
#include "../include/SAT_Const.h"
#include "../include/R_z.h"
#include "../include/R_y.h"


/**
 * Precession transformation of equatorial coordinates
 * @param Mjd_1 Epoch given (Modified Julian Date TT)
 * @param Mjd_2 Epoch to precess to (Modified Julian Date TT)
 * @return Precession transformation matrix
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2){

    double zeta;
    double z;
    double theta;
    double T;
    double dT;

    T  = (Mjd_1-MJD_J2000)/36525.0;
    dT = (Mjd_2-Mjd_1)/36525.0;

    //Precession angles
    zeta = ( (2306.2181+(1.39656-0.000139*T)*T)+ ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
    z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs;
    theta =  ( (2004.3109-(0.85330+0.000217*T)*T)- ((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;

    return R_z(-z) * R_y(theta) * R_z(-zeta);
}