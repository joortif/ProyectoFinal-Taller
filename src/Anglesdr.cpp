#include <valarray>
#include "../include/Matrix.h"
#include "../include/SAT_Const.h"
#include "../include/Geodetic.h"
#include "../include/LTC.h"
#include "../include/IERS.h"
#include "../include/Timediff.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/Norm.h"
#include "../include/Doubler.h"
#include "../include/Globals.h"


/**
 * Solves the problem of orbit determination using three optical sightings.
 *
 * @param[in] az1 Azimuth at t1 in radians.
 * @param[in] az2 Azimuth at t2 in radians.
 * @param[in] az3 Azimuth at t3 in radians.
 * @param[in] el1 Elevation at t1 in radians.
 * @param[in] el2 Elevation at t2 in radians.
 * @param[in] el3 Elevation at t3 in radians.
 * @param[in] Mjd1 Modified Julian Date of t1.
 * @param[in] Mjd2 Modified Julian Date of t2.
 * @param[in] Mjd3 Modified Julian Date of t3.
 * @param[in] rsite1 ijk site1 position vector in meters.
 * @param[in] rsite2 ijk site2 position vector in meters.
 * @param[in] rsite3 ijk site3 position vector in meters.
 * @param[out] r2 ijk position vector at t2 in meters.
 * @param[out] v2 ijk velocity vector at t2 in meters per second.
 */
void anglesdr(double az1, double az2, double az3, double el1,double el2,double el3,double Mjd1,double Mjd2,double Mjd3,
              Matrix rsite1,Matrix rsite2,Matrix rsite3,Matrix & r2, Matrix & v2){
    double magr1in, magr2in, tol, pctchg, t1, t3, lon1, lat1, h1, lon2, lat2, h2, lon3, lat3, h3, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole,
            UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_TT, Mjd_UT1, magr1old, magr2old, magrsite1, magrsite, magrsite3,
            cc1, cc2, ktr, magrsite2,f1,f2,q1,magr1,magr2,a,deltae32, f, g, magr1o, deltar1, f1delr1,f2delr1,q2,pf1pr1,pf2pr1, magr2o, deltar2,
            f1delr2,f2delr2,q3, pf1pr2, pf2pr2, delta, delta1, delta2;
    Matrix M1(3,3);
    Matrix M2(3,3);
    Matrix M3(3,3);
    Matrix P(3,1);
    Matrix N(3,1);
    Matrix T(3,1);
    Matrix E(3,3);
    Matrix los1(3,1);
    Matrix los2(3,1);
    Matrix los3(3,1);
    Matrix r3(3,1);

    char direct='y';

    magr1in = 1.1*R_Earth;
    magr2in = 1.11*R_Earth;
    tol = 1e-8*R_Earth;
    pctchg = 0.005;

    t1 = (Mjd1 - Mjd2)*86400.0;
    t3 = (Mjd3 - Mjd2)*86400.0;

    los1(0,0) = cos(el1)*sin(az1); los1(1,0) = cos(el1)*cos(az1); los1(2,0)=sin(el1);
    los2(0,0) = cos(el2)*sin(az2); los2(1,0) = cos(el2)*cos(az2); los2(2,0)=sin(el2);
    los3(0,0) = cos(el3)*sin(az3); los3(1,0) = cos(el3)*cos(az3); los3(2,0)=sin(el3);


    Geodetic(rsite1, lon1, lat1, h1);
    Geodetic(rsite2, lon2, lat2, h2);
    Geodetic(rsite3, lon3, lat3, h3);

    M1 = LTC(lon1, lat1);
    M2 = LTC(lon2, lat2);
    M3 = LTC(lon3, lat3);

    //body-fixed system
    los1 = M1.transpose()*los1;
    los2 = M1.transpose()*los2;
    los3 = M1.transpose()*los3;

    //mean of date system (JD2000)
    Mjd_UTC = Mjd1;
    IERS(eopdata, Mjd_UTC,UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los1 = E.transpose()*los1;
    rsite1 = E.transpose()*rsite1;

    Mjd_UTC = Mjd2;
    IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(MJD_J2000, Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    los2 = E.transpose() * los2;
    rsite2 = E.transpose() * rsite2;

    Mjd_UTC = Mjd3;
    IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    los3 = E.transpose()*los3;
    rsite3 = E.transpose()*rsite3;

    magr1old  = 99999999.9;
    magr2old  = 99999999.9;
    magrsite1 = norm(rsite1);
    magrsite2 = norm(rsite2);
    magrsite3 = norm(rsite3);

    cc1 = 2.0*los1.dot(rsite1);
    cc2 = 2.0*los2.dot(rsite2);
    ktr = 0;

    while (fabs(magr1in-magr1old)>tol || fabs(magr2in-magr2old)>tol){
        ktr++;
        doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,
                direct,r2,r3,f1,f2,q1,magr1,magr2,a,deltae32);
        f = 1.0 - a/magr2*(1.0-cos(deltae32));
        g  = t3 - sqrt(pow(a,3)/GM_Earth)*(deltae32-sin(deltae32));
        v2 = (r3 - r2*f)/g;

        magr1o = magr1in;
        magr1in = (1.0+pctchg)*magr1in;
        deltar1 = pctchg*magr1in;
        doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,
                 direct,r2,r3,f1delr1,f2delr1,q2,magr1,magr2,a,deltae32);
        pf1pr1 = (f1delr1-f1)/deltar1;
        pf2pr1 = (f2delr1-f2)/deltar1;

        magr1in = magr1o;
        deltar1 = pctchg*magr1in;
        magr2o = magr2in;
        magr2in = (1.0+pctchg)*magr2in;
        deltar2 = pctchg*magr2in;
        doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3
                 ,direct,r2,r3,f1delr2,f2delr2,q3,magr1,magr2,a,deltae32);
        pf1pr2 = (f1delr2-f1)/deltar2;
        pf2pr2 = (f2delr2-f2)/deltar2;

        magr2in = magr2o;
        deltar2 = pctchg*magr2in;

        delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
        delta1 = pf2pr2*f1 - pf1pr2*f2;
        delta2 = pf1pr1*f2 - pf2pr1*f1;

        deltar1 = -delta1/delta;
        deltar2 = -delta2/delta;

        magr1old = magr1in;
        magr2old = magr2in;

        magr1in = magr1in + deltar1;
        magr2in = magr2in + deltar2;

    }

    doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,rsite3,t1,t3,
             direct,r2,r3,f1,f2,q1,magr1,magr2,a,deltae32);

    f  = 1.0 - a/magr2*(1.0-cos(deltae32));
    g  = t3 - sqrt(pow(a,3)/GM_Earth)*(deltae32-sin(deltae32));
    v2 = (r3 - r2*f)/g;


}