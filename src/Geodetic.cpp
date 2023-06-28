#include <valarray>
#include "../include/Matrix.h"
#include "../include/SAT_Const.h"
#include "../include/Norm.h"

/**
 * Geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])  from given position vector (r [m])
 * @param r Position vector r[m]
 * @param lon Longitude geodetic coordinate [rad]
 * @param rad Latitude geodetic coordinate [rad]
 * @param m Latitude altitude coordinate [m]
 */
void Geodetic(Matrix r, double & lon, double & lat, double & h){
    double R_equ = R_Earth;
    double f = f_Earth;

    double epsRequ;
    double e2;
    double X;
    double Y;
    double Z;
    double rho2;
    double dZ;
    double ZdZ;
    double Nh;
    double SinPhi;
    double N;
    double dZ_new;

    epsRequ = eps*R_equ;        // Convergence criterion
    e2      = f*(2.0-f);        // Square of eccentricity
    X = r(0,0);            // Cartesian coordinates
    Y = r(1,0);
    Z = r(2,0);
    rho2 = X*X + Y*Y;           // Square of distance from z-axis

//  Check validity of input data
    if (norm(r) == -1) {
        throw "Invalid input in Geodetic constructor";
        lon = 0.0;
        lat = 0.0;
        h = -R_Earth;
    }

    // Iteration
    dZ = e2 * Z;

    while (true) {
        ZdZ = Z + dZ;
        Nh = sqrt(rho2 + ZdZ * ZdZ);
        SinPhi = ZdZ / Nh; // Sine of geodetic latitude
        N = R_equ / sqrt(1.0 - e2 * SinPhi * SinPhi);
        dZ_new = N * e2 * SinPhi;
        if (fabs(dZ - dZ_new) < epsRequ) {
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2(Y, X);
    lat = atan2(ZdZ, sqrt(rho2));
    h = Nh - N;

}