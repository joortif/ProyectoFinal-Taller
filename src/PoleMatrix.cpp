#include "../include/Matrix.h"
#include "../include/R_y.h"
#include "../include/R_x.h"

/**
 * Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
 * @param xp Pole coordinate X
 * @param yp Pole coordinate Y
 * @return Pole matrix
 */
Matrix PoleMatrix(double xp, double yp){
    Matrix PoleMat(3,3);

    PoleMat = R_y(-xp) * R_x(-yp);

    return PoleMat;
}