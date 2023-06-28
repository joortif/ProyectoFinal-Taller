#include "../include/Matrix.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"
#include "../include/R_x.h"
#include "../include/R_z.h"

/**
 * Transformation from mean to true equinator and equinox
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Nutation matrix
 */
Matrix NutMatrix(double Mjd_TT){
    double eps;
    double dpsi;
    double deps;
    Matrix NutMat(3,3);

    eps = MeanObliquity(Mjd_TT);
    NutAngles(Mjd_TT, dpsi, deps);

    NutMat = R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);
    return NutMat;

}