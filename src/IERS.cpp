#include <cmath>
#include "../include/IERS.h"
#include "../include/SAT_Const.h"

/**
 * Management of IERS time and polar motion data
 * @param eop The matrix containing the Earth Orientation Parameters.
 * @param Mjd_TT The Modified Julian Date, Terrestrial Time.
 * @param UT1_UTC [out] The UT1-UTC time difference in seconds.
 * @param TAI_UTC [out] The TAI-UTC time difference in seconds.
 * @param x_pole [out] The pole coordinate in arcseconds.
 * @param y_pole [out] The pole coordinate in arcseconds.
 */
void IERS(double ** eop, double Mjd_TT, double &UT1_UTC, double &TAI_UTC, double &x_pole, double &y_pole){
    double Arcs;
    double mj;
    int nop;

    Arcs = 3600.0*180.0/pi;
    mj = round(Mjd_TT);
    nop = 19716;

    for (int i=0; i<nop; i++){
        if (mj == round(eop[i][3])){
            // Setting of IERS Earth rotation parameters
            // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])

            UT1_UTC = eop[i][6];      // UT1-UTC time difference [s]
            TAI_UTC = eop[i][12];     // TAI-UTC time difference [s]
            x_pole  = eop[i][4]/Arcs; // Pole coordinate [rad]
            y_pole  = eop[i][5]/Arcs; // Pole coordinate [rad]
        }
    }

}