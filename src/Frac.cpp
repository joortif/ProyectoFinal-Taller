#include <cmath>

/**
* @brief Calculates the fractional part of a number (y=x-[x]). Last modified:   2015/08/12   M. Mahooti
*
* @param x The input number for which the fractional part needs to be calculated.
* @return The fractional part of x.
*/

double Frac(double x){
    return x-floor(x);
}