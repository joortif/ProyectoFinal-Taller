#include <cmath>
#include "../include/Matrix.h"

/**
 * Calculates the vector norm of a vector matrix.
 * @param v The vector matrix.
 * @return The vector norm of the matrix. If the matrix is not a valid vector, -1 is returned.
 */
double norm(Matrix v){

    if (v.cols() == 1){
        int tam = v.fils();
        double suma = 0.0;

        for (int i=0; i<tam; i++){
            suma += v(i,0) * v(i,0);
        }

        return sqrt(suma);

    }
    return -1;
}