#include "include/Matrix.h"

int main2(){
    double data[] = {1,2,3,4,-5,6,7,8,9};
    Matrix m = Matrix(3,3,data,9);
    m.inverse().print();
    Matrix::identity(3).print();

    return 0;
}