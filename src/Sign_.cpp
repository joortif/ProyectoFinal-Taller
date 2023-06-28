#include <cmath>
#include "../include/Sign_.h"

/**
 * Returns the absolute value of 'a' with the sign of 'b'
 * @param a Value for which the absolute value with the sign of 'b' is desired
 * @param b Value that determines the sign of the result
 * @return Result, which is the absolute value of 'a' with the sign of 'b'
 */
double sign_(double a, double b){
    if (b>=0.0){
        return fabs(a);
    }
    return -fabs(a);
}