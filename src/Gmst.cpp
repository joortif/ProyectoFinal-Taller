#include <cmath>
#include "../include/Frac.h"
#include "../include/SAT_Const.h"

/**
 * Greenwich Mean Sidereal Time
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return GMST in [rad]
 */
double gmst(double Mjd_UT1){
    double Secs = 86400.0;              //Seconds per day
    double MJD_J2000 = 51544.5;

    double Mjd_0;
    double UT1;
    double T_0;
    double T;
    double gmst;
    double gmstime;

    Mjd_0 = floor(Mjd_UT1);
    UT1 = Secs*(Mjd_UT1-Mjd_0);

    T_0   = (Mjd_0  -MJD_J2000)/36525.0;
    T     = (Mjd_UT1-MJD_J2000)/36525.0;

    gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1 + (0.093104-6.2e-6*T)*T*T;    // [s]

    gmstime = 2*pi*Frac(gmst/Secs);       // [rad], 0..2pi

    return gmstime;
}