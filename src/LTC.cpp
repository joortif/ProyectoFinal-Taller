#include "../include/Matrix.h"
#include "../include/R_y.h"
#include "../include/R_z.h"

/**
 * @brief Transformation from Greenwich meridian system to local tangent coordinate.
 * @param lon Geodetic East longitude in radians
 * @param lat Geodetic latitude in radians
 * @return The rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
 */
Matrix LTC(double lon, double lat){

    Matrix M=R_y(-1.0*lat)*R_z(lon);

    for(int j = 0; j < 3; j++){
        double Aux = M(0,j); M(0,j) = M(1,j); M(1,j) = M(2,j); M(2,j) = Aux; //Transpose M
    }

    return M;
}