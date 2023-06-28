#include "../include/Matrix.h"
#include "../include/Globals.h"
#include "../include/IERS.h"
#include "../include/Timediff.h"
#include "../include/PrecMatrix.h"
#include "../include/SAT_Const.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/AccelHarmonic.h"

/**
 * Computes the acceleration of an Earth orbiting satellite due to
 *    - the Earth's harmonic gravity field,
 *    - the gravitational perturbations of the Sun and Moon
 *    - the solar radiation pressure and
 *    - the atmospheric drag
 * @param x Terrestrial Time (Modified Julian Date)
 * @param Y Satellite state vector in the ICRF/EME2000 system
 * @return Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
 */
Matrix Accel(double x, Matrix Y){
    double UT1_UTC, TAI_UTC, x_pole, y_pole, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1;

    IERS(eopdata, auxParam.Mjd_TT + x/86400.0, UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    Mjd_UT1 = auxParam.Mjd_TT + x/86400.0 + (UT1_UTC-TT_UTC)/86400.0;

    Matrix P = PrecMatrix(MJD_J2000,auxParam.Mjd_TT+x/86400.0);
    Matrix N = NutMatrix(auxParam.Mjd_TT + x/86400.0);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    double dataSubY[] = {Y(0,0),Y(1,0),Y(2,0)};
    Matrix subY = Matrix(3,1,dataSubY,3);
    Matrix a = AccelHarmonic(subY, E, auxParam.n, auxParam.m);

    double datadY[] = {Y(3,0), Y(4,0), Y(5,0), a(0,0),a(1,0),a(2,0)};
    Matrix dY = Matrix(6,1,datadY,6);
    return dY;


}