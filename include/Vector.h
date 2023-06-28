#ifndef PROYECTO_VECTOR_H
#define PROYECTO_VECTOR_H

#include <iostream>

using namespace std;

class Vector {

private:
    double *data;
    int tam;

public:
    Vector(int n);

    Vector & operator=(Vector v);
    double& operator()(int i) const;
    void print();

};


#endif //PROYECTO_VECTOR_H
