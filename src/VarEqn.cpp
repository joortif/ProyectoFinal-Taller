#include "../include/Matrix.h"
#include "../include/IERS.h"
#include "../include/Globals.h"
#include "../include/Timediff.h"
#include "../include/PrecMatrix.h"
#include "../include/SAT_Const.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/AccelHarmonic.h"
#include "../include/G_AccelHarmonic.h"
#include "../include/NutMatrix.h"

/**
 * Computes the variational equations, i.e. the derivative of the state vector and the state transition matrix
 * @param x Time since epoch in [s]
 * @param yPhi (6+36)-dim vector comprising the state vector (y) and the state transition matrix (Phi) in column wise storage order
 * @return Derivative of yPhi
 */
Matrix VarEqn(double x, Matrix yPhi){
   double UT1_UTC, TAI_UTC, x_pole, y_pole, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1;

    IERS(eopdata, auxParam.Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_UT1 = auxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    //Transformation matrix
    Matrix P = PrecMatrix(MJD_J2000,auxParam.Mjd_TT + x/86400);
    Matrix N = NutMatrix(auxParam.Mjd_TT + x/86400);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    //State vector components
    double dataR[] = {yPhi(0,0),yPhi(1,0),yPhi(2,0)};
    double dataV[] = {yPhi(3,0), yPhi(4,0), yPhi(5,0)};
    Matrix r = Matrix(3,1,dataR,3);
    Matrix v = Matrix(3,1,dataV,3);
    Matrix Phi = Matrix(6,6);

    int k=0;
    //State transition matrix
    for (int i=0; i<6; i++){
        for (int j=1; j<=6; j++){
            Phi(i,j-1)=yPhi(6*j+k,0);
        }
        k++;
    }

    // Acceleration and gradient
    Matrix a = AccelHarmonic(r,E,auxParam.n, auxParam.m);
    Matrix G = G_AccelHarmonic(r,E,auxParam.n,auxParam.m);


    // Time derivative of state transition matrix
    Matrix yPhip = Matrix(42,1);
    Matrix dfdy = Matrix(6,6);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dfdy(i, j) = 0.0;                   // dv/dr(i,j)
            dfdy(i + 3, j) = G(i, j);           // da/dr(i,j)
            if (i == j) {
                dfdy(i, j + 3) = 1;
            }
            else {
                dfdy(i, j + 3) = 0.0;             // dv/dv(i,j)
            }
            dfdy(i + 3, j + 3) = 0.0;           // da/dv(i,j)
        }
    }

    Matrix Phip = dfdy*Phi;

    // Derivative of combined state vector and state transition matrix

    for (int i = 0; i < 3; i++) {
        yPhip(i, 0) = v(i, 0);                 // dr/dt(i)
        yPhip(i + 3, 0) = a(i, 0);             // dv/dt(i)
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 1; j <= 6; j++) {
            yPhip(6 * j + i, 0) = Phip(i, j-1);   // dPhi/dt(i,j)
        }
    }

    return yPhip;

}