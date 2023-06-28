#include "../include/Matrix.h"
#include <cmath>

/**
 * @brief Calculates the 3x3 rotation matrix around the Y-axis.
 * @param angle Rotation angle in radians
 * @return The rotation matrix
 */
Matrix R_y(double angle){
    double C;
    double S;
    Matrix rotmat = Matrix(3,3);

    C=cos(angle);
    S=sin(angle);

    rotmat(0,0) =   C;  rotmat(0,1) = 0.0;   rotmat(0,2) = -1.0*S;
    rotmat(1,0) = 0.0;  rotmat(1,1) = 1.0;   rotmat(1,2) =    0.0;
    rotmat(2,0) =   S;  rotmat(2,1) = 0.0;   rotmat(2,2) =      C;

    return rotmat;

}
