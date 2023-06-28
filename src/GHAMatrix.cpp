#include "../include/Matrix.h"
#include "../include/R_z.h"
#include "../include/Gast.h"

/**
 * Transformation from true equator and equinox to Earth equator and Greenwich meridian system
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return Greenwich Hour Angle matrix
 */
Matrix GHAMatrix(double Mjd_UT1){
    Matrix GHAMatrix(3,3);

    GHAMatrix = R_z(gast(Mjd_UT1));
    return GHAMatrix;
}