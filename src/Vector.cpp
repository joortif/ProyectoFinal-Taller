#include "../include/Vector.h"

using namespace std;

/**
 * @brief Constructor of the Vector class.
 *
 * @param n The size of the vector.
 */
Vector::Vector(int n) {
    this->tam = n;
    data = new double[tam];
}

/**
 * @brief Overloaded assignment operator.
 *
 * @param v The vector to assign to.
 * @return Vector& The reference to the current vector.
 */
Vector & Vector::operator=(Vector v) {
    this->tam = v.tam;
    this->data = new double[tam];
    for (int i=0; i<tam; i++){
        data[i] = v.data[i];
    }
    return *this;
}

/**
 * @brief Overloaded operator that returns a reference to an element in the vector.
 *
 * @param i The index of the element.
 * @return double& The reference to the element.
 */
double &Vector::operator()(int i) const{
    return data[i-1];
}

/**
 * @brief Prints the vector to the console.
 */
void Vector::print() {
    cout << "[";
    for (int i=0; i<this->tam; i++){
        cout << this->data[i];
        if (i != tam -1) {
            cout << ", ";
        }
    }
    cout << "]" << endl;

}

