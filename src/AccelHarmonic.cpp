#include <valarray>
#include <iostream>
#include "../include/Matrix.h"
#include "../include/Norm.h"
#include "../include/Legendre.h"
#include "../include/Globals.h"
#include <iostream>

using namespace std;

/**
 * @brief Computes the acceleration due to the harmonic gravity field of the central body.
 *
 * This function calculates the acceleration due to the harmonic gravity field of the central body
 * using the provided inputs.
 *
 * @param Mjd_TT Modified Julian Date of TT.
 * @param r Satellite position vector in the inertial system.
 * @param E Transformation matrix to body-fixed system.
 * @param n_max Maximum degree.
 * @param m_max Maximum order (m_max <= n_max; m_max = 0 for zonals only).
 * @return Acceleration (a = d^2r/dt^2).
 */
Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max){

    double d, latgc, lon, dUdr, dUdlatgc, dUdlon, q3, q2, q1, b1, b2, b3, r2xy, ax, ay, az;
    Matrix pnm = Matrix(1,1);
    Matrix dpnm = Matrix(1,1);

    double gm=398600.4415e9;              // [m^3/s^2]; JGM3/EGM96
    double r_ref = 6378.1363e3;           // Radius Earth [m]; JGM3/EGM96

    // Body-fixed position
    Matrix r_bf = E * r;
    d = norm(r_bf);
    latgc = asin(r_bf(2,0)/d);
    lon = atan2(r_bf(1,0),r_bf(0,0));

    Legendre(n_max, m_max, latgc, pnm, dpnm);

    dUdr = 0;
    dUdlatgc = 0;
    dUdlon = 0;
    q3 = 0; q2 = q3; q1 = q2;
    for (int n = 0; n<=n_max; n++){
        b1 = (-gm/(d*d))*pow(r_ref/d,n)*(n+1);
        b2 =  (gm/d)*pow(r_ref/d,n);
        b3 =  (gm/d)*pow(r_ref/d,n);
        for (int m = 0; m<=n; m++){
            q1 += pnm(n,m)*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
            q2 += dpnm(n,m)*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
            q3 += m*pnm(n,m)*(Snm[n][m]*cos(m*lon)-Cnm[n][m]*sin(m*lon));
        }
        dUdr+=q1*b1;
        dUdlatgc+=q2*b2;
        dUdlon+=q3*b3;
        q3=0;q2=q3;q1=q2;
    }

    //Body-fixed acceleration
    r2xy = r_bf(0,0)*r_bf(0,0)+r_bf(1,0)*r_bf(1,0);
    ax = (1/d*dUdr-r_bf(2,0)/((d*d)*sqrt(r2xy))*dUdlatgc)*r_bf(0,0)-(1/r2xy*dUdlon)*r_bf(1,0);
    ay = (1/d*dUdr-r_bf(2,0)/((d*d)*sqrt(r2xy))*dUdlatgc)*r_bf(1,0)+(1/r2xy*dUdlon)*r_bf(0,0);
    az =  1/d*dUdr*r_bf(2,0)+sqrt(r2xy)/(d*d)*dUdlatgc;

    double datos[] = {ax,ay,az};
    Matrix a_bf(3,1,datos,3);

    return E.transpose()*a_bf;



}