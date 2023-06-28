#include <valarray>
#include "../include/Gmst.h"
#include "../include/SAT_Const.h"
#include "../include/EqnEquinox.h"

/**
 * Greenwich Apparent Sidereal Time
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return GAST in [rad]
 */
double gast(double Mjd_UT1){
    double gstime;
    double a = gmst(Mjd_UT1);
    double b = EqnEquinox(Mjd_UT1);
    double res = fmod(a+b, 2*pi);
    gstime = res;
    return gstime;
}