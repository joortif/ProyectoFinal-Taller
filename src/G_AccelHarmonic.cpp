#include "../include/Matrix.h"
#include "../include/AccelHarmonic.h"

/**
 * Computes the gradient of the Earth's harmonic gravity field.
 *
 * @param r       Satellite position vector in the true-of-date system.
 * @param U       Matrix U (modified Julian Date).
 * @param n_max   Gravity model degree.
 * @param m_max   Gravity model order.
 *
 * @return        The gradient (G=da/dr) in the true-of-date system as a 3x3 matrix.
 */
Matrix G_AccelHarmonic(Matrix &r, Matrix &U, int n_max, int m_max){
    double d = 1.0;
    Matrix G(3,3);
    Matrix dr(3,1);

    for (int i=0;i<3;i++){
        dr(0,0) = 0;
        dr(1,0) = 0;
        dr(2,0) = 0;
        dr(i,0)=d;
        Matrix rPlusdr = r+dr/2;
        Matrix rMinusdr = r-dr/2;
        Matrix da = AccelHarmonic( rPlusdr,U, n_max, m_max ) -
                AccelHarmonic(rMinusdr,U, n_max, m_max );

        for (int j = 0; j < 3; j++) {
            G(j, i) = da(j, 0) / d;
        }
    }
    return G;
}