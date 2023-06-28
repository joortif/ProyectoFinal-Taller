#include <cmath>

#include "../include/Matrix.h"
#include "../include/SAT_Const.h"

/**
 *
 * @brief Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 * @param lon The longitude in radians.
 * @param lat The latitude in radians.
 * @param h The altitude in meters.
 * @return A 3D vector (as a Vector object) representing the position in meters.
 */
Matrix Position(double lon, double lat, double h){

    double R_equ;
    double f;
    double e2;
    double CosLat;
    double SinLat;
    double N;
    Matrix r(3,1);

    R_equ = R_Earth;
    f = f_Earth;
    e2= f*(2.0-f);
    CosLat = cos(lat);
    SinLat = sin(lat);


    N = R_equ / sqrt(1.0-e2*SinLat*SinLat);
    r(0,0) = (N+h)*CosLat*cos(lon);
    r(1,0) = (N+h)*CosLat*sin(lon);
    r(2,0) = ((1.0-e2)*N+h)*SinLat;

    return r;
}