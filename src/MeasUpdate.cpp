#include "../include/Matrix.h"
#include <iostream>

/**
 * Performs the measurement update step in a Kalman filter by updating the state and covariance matrix.
 *
 * @param x   The state vector to be updated.
 * @param z   The measurement vector.
 * @param g   The measurement function.
 * @param s   The measurement noise standard deviation vector.
 * @param G   The measurement function Jacobian matrix.
 * @param P   The covariance matrix to be updated.
 * @param n   The dimension of the state vector.
 *
 * @return The Kalman gain matrix, the updated state vector, and the updated covariance matrix.
 */

Matrix MeasUpdate(Matrix& x, Matrix& z, Matrix& g, Matrix& s, Matrix& G, Matrix& P, int n){
    int m = z.fils();
    Matrix Inv_W(m,m);

    for (int i=0; i<m; i++) {
        Inv_W(i, i) = s(i, 0) * s(i, 0); //Inverse weight (measurement covariance)
    }
    Matrix K = P * G.transpose() * (Inv_W + G * P *G.transpose()).inverse();

    x = x+ K * (z-g);
    P = (Matrix::identity(n)-K*G)*P;
    return K;

}