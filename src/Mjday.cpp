#include <cmath>
#include "../include/Mjday.h"

/**
 * @brief Calculates the Modified Julian Date from a given date and time.
 *
 * @param year The year of the date.
 * @param mon The month of the date.
 * @param day The day of the date.
 * @param hr The hour of the time in Universal Time.
 * @param min The minute of the time in Universal Time.
 * @param sec The second of the time in Universal Time.
 * @return double The Modified Julian Date.
 */
double Mjday(int yr, int mon, int day, int hr=0, int min=0, int sec=0){

    double jd = 367.0 *yr - floor((7*(yr + floor((mon+9)/12.0))) *0.25) + floor(275 * mon / 9.0) + day + 1721013.5 + ((sec/60.0 + min) / 60.0 + hr) / 24.0;

    return jd-2400000.5;
}

