#include <valarray>
#include "../include/Matrix.h"
#include "../include/SAT_Const.h"

/**
 * Computes azimuth, elevation and partials from local tangent coordinates
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame)
 * @param Az Azimuth [rad]
 * @param El Elevation [rad]
 * @param dAds Partials of azimuth w.r.t. s
 * @param dEds Partials of elevation w.r.t. s
 */
void AzElPa(Matrix& s, double &Az, double &El, Matrix &dAds, Matrix &dEds){

    double pi2 = 2.0*pi;
    double rho = sqrt(s(0,0)*s(0,0)+s(1,0)*s(1,0));

    Az = atan2(s(0,0),s(1,0));

    if (Az<0.0){
        Az+=pi2;
    }

    El = atan(s(2,0)/rho);

    double datadAds[] = {s(1,0)/(rho*rho),-s(0,0)/(rho*rho),0.0};
    double datadEds[] = {(-s(0,0)*s(2,0)/rho)/s.dot(s), (-s(1,0)*s(2,0)/rho)/s.dot(s), rho/s.dot(s)};

    dAds = Matrix(3,1,datadAds,3);
    dEds = Matrix(3,1,datadEds, 3);
}