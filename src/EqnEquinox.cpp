#include <cmath>
#include "../include/NutAngles.h"
#include "../include/MeanObliquity.h"

/**
 * Computation of the equation of equinoxes.
 * The equation of the equinoxes dpsi*cos(eps) is the right ascension of
 * the mean equinox referred to the true equator and equinox and is equal
 * to the difference between apparent and mean sidereal time.
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Equation of equinoxes
 */
double EqnEquinox(double Mjd_TT){
    double dpsi;
    double deps;
    double EqE;

    NutAngles(Mjd_TT, dpsi, deps);
    EqE = dpsi * cos(MeanObliquity(Mjd_TT));
    return EqE;
}