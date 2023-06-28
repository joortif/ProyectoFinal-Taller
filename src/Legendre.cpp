#include <valarray>
#include "../include/Matrix.h"

/**
 * Computes the Legendre polynomials and their derivatives up to order (n, m) at a given angle (fi).
 *
 * @param n    The maximum degree of the Legendre polynomial.
 * @param m    The maximum order of the Legendre polynomial.
 * @param fi   The angle in radians.
 * @param pnm  Output matrix containing the Legendre polynomials.
 * @param dpnm Output matrix containing the derivatives of the Legendre polynomials.
 */
void Legendre(int n, int m, double fi, Matrix &pnm, Matrix &dpnm){
    int j=0;
    int k=2;
    if (n>m){
        pnm = Matrix(n+1, n+1);
        dpnm = Matrix(n+1,n+1);
    } else {
        pnm = Matrix(n+1, m+1);
        dpnm = Matrix(n+1,m+1);
    }

    pnm(0,0) = 1.0;
    dpnm(0,0) = 0.0;
    pnm(1,1)=sqrt(3)*cos(fi);
    dpnm(1,1)=-sqrt(3)*sin(fi);

    // diagonal coefficients
    for (double i = 1; i < n; i++) {
        pnm(i+1, i+1) = sqrt((2*(i+1)+1.0)/(2*(i+1))) * cos(fi) * pnm(i, i);
    }

    for (double i = 1; i < n; i++) {
        dpnm(i+1, i+1) = sqrt((2*(i+1)+1.0)/(2*(i+1))) * ((cos(fi) * dpnm(i, i)) - (sin(fi) * pnm(i, i)));
    }

    // horizontal first step coefficients
    for (double i = 0; i < n; i++) {
        pnm(i+1, i) = sqrt(2*(i+1)+1.0) * sin(fi) * pnm(i, i);
    }

    for (double i = 0; i < n; i++) {
        dpnm(i+1, i) = sqrt(2*(i+1)+1.0) * ((cos(fi) * pnm(i, i)) + (sin(fi) * dpnm(i, i)));
    }

    // horizontal second step coefficients
    while(true){
        for (double i=k;i<=n; i++){
            pnm(i,j)=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*pnm(i-1,j))
            -(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*pnm(i-2,j)));
        }
        j++;
        k++;
        if (j>m){
            break;
        }
    }
    j=0;
    k=2;

    while (true){
        for (double i=k;i<=n;i++){
            dpnm(i,j)=sqrt((2.0*i+1.0)/((i-j)*(i+j)))*((sqrt(2.0*i-1.0)*sin(fi)*dpnm(i-1,j))
            +(sqrt(2*i-1)*cos(fi)*pnm(i-1,j))-(sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3.0))*dpnm(i-2,j)));
        }
        j++;
        k++;
        if (j>m){
            break;
        }
    }
}