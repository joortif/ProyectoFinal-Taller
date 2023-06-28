#include <valarray>
#include "../include/Matrix.h"
#include "../include/Norm.h"
#include "../include/SAT_Const.h"

/**
 * Performs calculations for orbit propagation and orbital mechanics.
 *
 * @param cc1 Constant value used in the calculations.
 * @param cc2 Constant value used in the calculations.
 * @param magrsite1 Magnitude of site position 1.
 * @param magrsite2 Magnitude of site position 2.
 * @param magr1in Magnitude of initial position 1.
 * @param magr2in Magnitude of initial position 2.
 * @param los1 Line of sight position vector 1.
 * @param los2 Line of sight position vector 2.
 * @param los3 Line of sight position vector 3.
 * @param rsite1 Site position vector 1.
 * @param rsite2 Site position vector 2.
 * @param rsite3 Site position vector 3.
 * @param t1 Time 1.
 * @param t3 Time 3.
 * @param direct Direction indicator.
 * @param r2 Resultant position vector 2 (modified output).
 * @param r3 Resultant position vector 3 (modified output).
 * @param f1 Resultant value f1 (modified output).
 * @param f2 Resultant value f2 (modified output).
 * @param q1 Resultant value q1 (modified output).
 * @param magr1 Magnitude of position 1 (modified output).
 * @param magr2 Magnitude of position 2 (modified output).
 * @param a Resultant value a (modified output).
 * @param deltae32 Resultant value deltae32 (modified output).
 */
void doubler(double cc1, double cc2, double magrsite1, double magrsite2, double magr1in, double magr2in,
             Matrix los1, Matrix los2, Matrix los3, Matrix rsite1, Matrix rsite2, Matrix rsite3,
             double t1, double t3, char direct, Matrix &r2, Matrix &r3, double &f1, double &f2,
             double &q1, double &magr1, double &magr2, double &a, double &deltae32){


    double rho1, rho2, rho3, magr3, cosdv21, sindv21, dv21, cosdv31, sindv31, dv31,
            cosdv32, sindv32, dv32, c1, c3, p, ecosv1, ecosv2, ecosv3, esinv2, e, n, s, c, sinde32, cosde32,
            sinde21, cosde21, deltae21, deltam32, deltam12, sindh32, sindh21, deltah32, deltah21;
    Matrix r1(3,1);
    Matrix w(3,1);

    rho1 = (-cc1 + sqrt(cc1*cc1-4*(magrsite1*magrsite1-magr1in*magr1in))) / 2.0;
    rho2 = (-cc2 + sqrt(cc2*cc2-4*(magrsite2*magrsite2-magr2in*magr2in))) / 2.0;

    r1 = los1 * rho1 + rsite1;
    r2 = los2 * rho2 + rsite2;

    magr1 = norm(r1);
    magr2 = norm(r2);

    if (direct == 'y'){
        w = r1.cross(r2) / (magr1 * magr2);

    } else {
        w = r1.cross(r2)*-1 / (magr1*magr2);
    }

    rho3 = -1*rsite3.dot(w) / los3.dot(w);
    r3 = los3 * rho3 + rsite3;
    magr3 = norm(r3);

    cosdv21 = r2.dot(r1) / (magr2 * magr1);
    sindv21 = norm(r2.cross(r1)) / (magr2 * magr1);
    dv21 = atan2(sindv21, cosdv21);

    cosdv31 = r3.dot(r1) / (magr3 * magr1);
    sindv31 = sqrt(1.0 - cosdv31 * cosdv31);
    dv31 = atan2(sindv31, cosdv31);

    cosdv32 = r3.dot(r2) / (magr3 * magr2);
    sindv32 = norm(r3.cross(r2)) / (magr3 * magr2);
    dv32 = atan2(sindv32, cosdv32);

    if (dv31 > pi){
        c1 = (magr2 * sindv32) / (magr1 * sindv31);
        c3 = (magr2 * sindv21) / (magr3 * sindv31);
        p = (c1 * magr1 + c3 * magr3 - magr2) / (c1 + c3 - 1);
    } else {
        c1 = (magr1 * sindv31) / (magr2 * sindv32);
        c3 = (magr1 * sindv21) / (magr3 * sindv32);
        p = (c3 * magr3 - c1 * magr2 + magr1) / (-c1 + c3 + 1);
    }

    ecosv1 = p /magr1 - 1;
    ecosv2 = p /magr2 - 1;
    ecosv3 = p /magr3 - 1;

    if (dv21 != pi)
        esinv2 = (-cosdv21 * ecosv2 + ecosv1) / sindv21;
    else
        esinv2 = (cosdv32 * ecosv2 - ecosv3) / sindv31;

    e = sqrt(ecosv2 * ecosv2 + esinv2 * esinv2);
    a = p / (1 - e * e);


    if (e * e < 0.99)
    {
        n = sqrt(GM_Earth / (a * a * a));

        s = magr2 / p * sqrt(1 - e * e) * esinv2;
        c = magr2 / p * (e * e + ecosv2);

        sinde32 = magr3 / sqrt(a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s;
        cosde32 = 1 - magr2 * magr3 / (a * p) * (1 - cosdv32);
        deltae32 = atan2(sinde32, cosde32);

        sinde21 = magr1 / sqrt(a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s;
        cosde21 = 1 - magr2 * magr1 / (a * p) * (1 - cosdv21);
        deltae21 = atan2(sinde21, cosde21);

        deltam32 = deltae32 + 2 * s * pow(sin(deltae32 / 2), 2) - c * sin(deltae32);
        deltam12 = -deltae21 + 2 * s * pow(sin(deltae21 / 2), 2) + c * sin(deltae21);

    } else {
        n = sqrt(GM_Earth / pow(-a,3));

        s = magr2 / p * sqrt(e * e - 1) * esinv2;
        c = magr2 / p * (e * e + ecosv2);

        sindh32 = magr3 / sqrt(-a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s;
        sindh21 = magr1 / sqrt(-a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s;

        deltah32 = log(sindh32 + sqrt(sindh32 * sindh32 + 1));
        deltah21 = log(sindh21 + sqrt(sindh21 * sindh21 + 1));

        deltam32 = -deltah32 + 2 * s * pow(sinh(deltah32 / 2), 2) + c * sinh(deltah32);
        deltam12 = deltah21 + 2 * s * pow(sinh(deltah21 / 2), 2) - c * sinh(deltah21);
        deltae32 = deltah32;

    }


    f1 = t1 - deltam12 / n;
    f2 = t3 - deltam32 / n;
    q1 = sqrt(f1*f1+f2*f2);
}
