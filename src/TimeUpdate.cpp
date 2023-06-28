#include "../include/Matrix.h"

/**
 * Performs the time update step in a Kalman filter by updating the covariance matrix.
 *
 * @param P    The covariance matrix to be updated.
 * @param Phi  The state transition matrix.
 * @param Qdt  The process noise covariance matrix. (Optional: default value is zero matrix)
 *
 * @return The updated covariance matrix.
 */
Matrix TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt){
    Matrix Qdt2 = Qdt;
    if (Qdt.cols() == 0 || Qdt.fils() == 0){
        Qdt2 = Matrix(P.fils(), P.cols());
    }

    Matrix P2 = Phi * P * Phi.transpose() + Qdt2;
    return P2;
}