#include <fstream>
#include <iostream>

#include "include/Globals.h"
#include "include/SAT_Const.h"
#include "include/Mjday.h"
#include "include/Position.h"
#include "include/Anglesdr.h"
#include "include/Accel.h"
#include "include/IERS.h"
#include "include/Timediff.h"
#include "include/LTC.h"
#include "include/Gmst.h"
#include "include/R_z.h"
#include "include/TimeUpdate.h"
#include "include/AzElPa.h"
#include "include/MeasUpdate.h"
#include "include/Norm.h"
#include "include/VarEqn.h"
#include "include/DEInteg.h"

using namespace std;

double **eopdata = nullptr;
double **Cnm = nullptr;
double **Snm = nullptr;
AuxParam auxParam;

//--------------------------------------------------------------------------
//
//  Initial Orbit Determination using Double-R-Iteration and Extended Kalman
//  Filter methods
//
// Last modified:   2015/08/12   M. Mahooti
//
// Refrences:
//
//   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
//   Applications", Springer Verlag, Heidelberg, 2000
//
//   D. Vallado, "Fundamentals of Astrodynamics and Applications",
//   3rd Edition, 2007
//
//   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003
//
//--------------------------------------------------------------------------
int main(){

    double sigma_range, sigma_az, sigma_el, lat, lon, alt, Mjd1, Mjd2, Mjd3, Mjd0,
            UT1_UTC, TAI_UTC, x_pole, y_pole, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC,
            Mjd_UTC, t, t_old, Mjd_TT, Mjd_UT1, theta, Azim, Elev, Dist;
    Matrix r2 = Matrix(3,1);
    Matrix v2 = Matrix(3,1);
    Matrix YM = Matrix(1,1);
    Cnm = new double*[361];
    for (int i=0; i<361; i++){
        Cnm[i] = new double[361];
    }
    Snm = new double*[361];
    for (int i=0; i<361; i++){
        Snm[i] = new double[361];
    }
    double temp[6];

    ifstream fid("data/egm.txt");
    if (!fid.is_open()){
        cerr << "ERROR: eegm.txt NOT FOUND" << endl;
    }

    for (int n = 0; n< 360; n++) {
        for (int j = 0; j <= n; j++) {
            for (int i = 0; i < 6; i++) {
                fid >> temp[i];
            }
            Cnm[n][j] = temp[2];
            Snm[n][j] = temp[3];
        }
    }

    fid.close();

    eopdata = new double*[19716];
    for (int i = 0; i < 19716; i++) {
        eopdata[i] = new double[13];
    }

    ifstream fid2("data/eop19620101.txt");
    if (!fid2.is_open()){
        cerr << "ERROR: eop19620101.txt NOT FOUND" << endl;
    }

    /*
     ----------------------------------------------------------------------------------------------------
     | Date MJD x y UT1-UTC LOD dPsi dEpsilon dX dY DAT
     |(0h UTC) " " s s " " " " s
     ----------------------------------------------------------------------------------------------------
    */
    for (int i = 0; i < 19716; i++) {
        for (int j = 0; j < 13; j++) {
            fid2 >> eopdata[i][j];
        }
    }

    fid2.close();

    int nobs=18;
    string tline;
    int Y, M, D, h, m;
    float s, az, el;

    double **obs = new double*[nobs];
    for (int i=0; i<nobs; i++){
        obs[i] = new double[4];
    }


    ifstream fid3("data/GEOS3.txt");
    if (!fid3.is_open()){
        cerr << "ERROR: GEOS3.txt NOT FOUND" << endl;
    }

    for (int i=0; i<nobs; i++){
        getline(fid3, tline);
        if (tline.empty()){
            break;
        }

        Y = stoi(tline.substr(0,4));
        M = stoi(tline.substr(5,2));
        D = stoi(tline.substr(8,2));
        h = stoi(tline.substr(12,2));
        m = stoi(tline.substr(15,2));
        s = stof(tline.substr(18,5));
        az = stof(tline.substr(25,8));
        el = stof(tline.substr(35,8));
        Dist = stof(tline.substr(44,9));
        obs[i][0] = Mjday(Y,M,D,h,m,s);
        obs[i][1] = Rad*az;
        obs[i][2] = Rad*el;
        obs[i][3] = 1e3*Dist;
    }

    fid3.close();

    sigma_range = 92.5;         // [m]
    sigma_az = 0.0224*Rad;      // [rad]
    sigma_el = 0.0139*Rad;      // [rad]

    //Kaena Point station
    lat = Rad*21.5748;          // [rad]
    lon = Rad*(-158.2706);      // [rad]
    alt = 300.20;               // [m]

    Matrix Rs = Position(lon, lat, alt);

    Mjd1 = obs[0][0];
    Mjd2 = obs[8][0];
    Mjd3 = obs[17][0];

    anglesdr(obs[0][1],obs[8][1],obs[17][1],obs[0][2],obs[8][2],obs[17][2],Mjd1,Mjd2,Mjd3,
             Rs,Rs,Rs, r2, v2);

    double dataY0[] = {r2(0,0),r2(1,0),r2(2,0),v2(0,0),v2(1,0),v2(2,0)};
    Matrix Y0_apr = Matrix(6,1,dataY0,6);

    Y0_apr.print();
    //Final de la primera parte del proyecto

    /*
    Mjd0 = Mjday(1995,1,29,2,38,0.0);

    Mjd_UTC = obs[8][0];
    IERS(eopdata, Mjd_UTC,UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC,UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    auxParam.Mjd_TT = Mjd_UTC + TT_UTC/86400;
    auxParam.n = 10;
    auxParam.m = 10;

    n_eqn = 3;

    Matrix YDE = Matrix(1,1);
    YDE = DEInteg (Accel,0,-(obs[9][1]-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

    Matrix P = Matrix(6,6);

    for (int i=0; i<3; i++){
        P(i,i) = 1e8;
    }
    for (int i=3; i<6; i++){
        P(i,i) = 1e3;
    }

    Matrix LT = LTC(lon, lat);
    Matrix yPhi = Matrix(42,1);
    Matrix Phi = Matrix(6,6);

    // Measurement loop
    t=0;

    Matrix Y_old = Matrix(1,1);
    for (int i=0; i<nobs; i++){
        //Previous step
        t_old = t;
        Y_old = YDE;

        // Time increment and propagation
        Mjd_UTC = obs[i][0];            // Modified Julian Date
        t = (Mjd_UTC-Mjd0)*86400.0;     // Time since epoch [s]

        IERS(eopdata, Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole);
        timediff(UT1_UTC, TAI_UTC,UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        auxParam.Mjd_UTC = Mjd_UTC;
        auxParam.Mjd_TT = Mjd_TT;

        for (int ii=0;ii<6;ii++){
            yPhi(ii,0) = Y_old(ii,0);
            for (int j=0;j<6;j++){
                if (ii==j){
                    yPhi(6*j+ii,0) = 1; //POSIBLE ERROR
                } else {
                    yPhi(6*j+ii,0) = 0;
                }
            }
        }

        yPhi = DEInteg(VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);

        //Extract state transition matrices
        for (int j=0; j<6; j++){
            for (int k = 0; k < 6; k++) {
                Phi(k,j) = yPhi(6 * j + k + 1, 0); //POSIBLE ERROR
            }
        }

        YM = DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old);

        // Topocentric coordinates
        theta = gmst(Mjd_UT1);
        Matrix U = R_z(theta);
        double dataR[] = {YDE(0,0),YDE(1,0), YDE(2,0)};
        Matrix r = Matrix(3,1,dataR,3);
        Matrix s = LT*(U*r-Rs);                          // Topocentric position [m]

        Matrix none = Matrix(0,0);
        // Time update
        P = TimeUpdate(P, Phi, none);
        Matrix dAds = Matrix(1,1);
        Matrix dEds = Matrix(1,1);

        // Azimuth and partials
        AzElPa(s, Azim, Elev, dAds, dEds);
        Matrix dAdY = Matrix(6,1);
        dAdY = dAds*LT*U;

        Matrix K = Matrix(1,1);
        double dataAz[] = {Azim};
        Matrix AzimD = Matrix(1,1,dataAz,1);
        double datasigma[] = {sigma_az};
        Matrix sigma_azD = Matrix(1,1,datasigma,1);
        // Measurement update
        Matrix m3 = Matrix(1,1);
        m3(0,0) = obs[i][1];
        K = MeasUpdate(YDE, m3, AzimD, sigma_azD, dAdY, P, 6);

        // Elevation and partials

    }

    IERS(eopdata, obs[17][0], UT1_UTC, TAI_UTC, x_pole, y_pole);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    auxParam.Mjd_UTC = Mjd_UTC;
    auxParam.Mjd_TT = Mjd_TT;

    YDE = DEInteg(Accel,0,-(obs[18][1]-obs[1][1])*86400.0,1e-13,1e-6,6,YM);
    double dataY_true[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
    Matrix Y_true = Matrix(6,1,dataY_true,6);

    cout << "Error of Position Estimation" << endl;
    cout << "[m] " << YDE(0,0)-Y_true(0,0) << endl;
    cout << "[m] " << YDE(1,0)-Y_true(1,0) << endl;
    cout << "[m] " << YDE(2,0)-Y_true(2,0) << endl;
    cout << "Error of Position Estimation" << endl;
    cout << "[m/s] " << YDE(3,0)-Y_true(3,0) << endl;
    cout << "[m/s] " << YDE(4,0)-Y_true(4,0) << endl;
    cout << "[m/s] " << YDE(5,0)-Y_true(5,0) << endl;
    */

    return 0;
}