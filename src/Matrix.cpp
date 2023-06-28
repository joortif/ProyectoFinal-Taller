#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>

/**
    * @brief Constructor that creates a matrix with the given number of rows and columns.
    * @param fil The number of rows.
    * @param col The number of columns.
*/
Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}

/**
    * @brief Constructor that creates a matrix with the given number of rows and columns, and fills it with values from a given array.
    * @param fil The number of rows.
    * @param col The number of columns.
    * @param v An array of values to fill the matrix with.
    * @param n The number of values in the array.
*/
Matrix::Matrix(int fil, int col, const double v[], int n): fil(fil), col(col)
{
    initMatrix();

    int k = 0;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

/**
    * @brief Copy constructor.
    * @param m The matrix to be copied.
*/
Matrix::Matrix(const Matrix& m)
{
    *this = m;
}

/**
    * @brief Destructor that frees the memory used by the matrix.
*/
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];

    delete[] matrix;
}


/**
    * @brief Initializes the matrix with all elements set to 0.
*/
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}

/**
    * @brief Overloaded assignment operator.
    * @param matrix2 The matrix to be assigned.
    * @return A reference to the newly assigned matrix.
*/
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    col = matrix2.col;
    fil = matrix2.fil;
    this->initMatrix();

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];

    return *this;
}

/**
    * @brief Overloaded addition operator that adds two matrices element-wise.
    * @param matrix2 The matrix to be added.
    * @return A new matrix that is the sum of the two input matrices.
*/
Matrix Matrix::operator+(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];

    return result;
}

/**
    * @brief Overloaded subtraction operator that subtracts two matrices element-wise.
    * @param matrix2 The matrix to be subtracted.
    * @return A new matrix that is the difference between the two input matrices.
*/
Matrix Matrix::operator-(const Matrix& matrix2)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];

    return result;
}

/**
    * @brief Overloaded multiplication operator that multiplies two matrices together.
    * @param matrix2 The matrix to be multiplied.
    * @return A new matrix that is the product of the two input matrices.
*/
Matrix Matrix::operator*(const Matrix& matrix2)
{
    Matrix result(fil, matrix2.col);

    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }

    return result;
}

/**
 * Multiplies a Matrix object by a scalar value.
 *
 * @param n The scalar value to multiply the matrix by.
 * @return The resulting matrix after the multiplication.
 */
Matrix Matrix::operator*(const double n){
    Matrix res(fil, col);

    for (int i=0; i< this->fil; i++){
        for (int j = 0; j < col; j++){
            res.matrix[i][j] = this->matrix[i][j]*n;
        }
    }
    return res;
}

/**
 * Performs element-wise division of the matrix by a scalar value.
 *
 * @param n The scalar value to divide the matrix by.
 * @return A new Matrix object resulting from the division.
 */
Matrix Matrix::operator/(const double n){
    Matrix res(fil, col);


    for (int i = 0; i < this->fil; i++){
        for (int j = 0; j < col; j++){
            res.matrix[i][j] = this->matrix[i][j] / n;
        }
    }
    return res;
}

/**
 * Negates all elements of the matrix.
 *
 * @return A new Matrix object with negated elements.
 */
Matrix Matrix::operator-() {
    Matrix res(fil, col);

    for (int i = 0; i < this->fil; i++){
        for (int j = 0; j < col; j++){
            res.matrix[i][j] = -this->matrix[i][j];
        }
    }
    return res;
}

/**
    * @brief Overloaded parentheses operator that allows access to individual elements of the matrix.
    * @param i The row index.
    * @param j The column index.
    * @return A reference to the element at the given indices.
*/
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i][j];
}

/**
    * @brief Prints the matrix to standard output.
*/
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::setprecision(30) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/**
 * Retrieves the number of rows in the matrix.
 * @return The number of rows in the matrix.
 */
int Matrix::fils() {
    return this->fil;
}

/**
 * Retrieves the number of columns in the matrix.
 * @return The number of columns in the matrix.
 */
int Matrix::cols() {
    return this->col;
}

/**
 * Calculates the dot product between two vector matrices.
 * @param matrix2 The second matrix for the dot product.
 * @return The dot product result as a scalar value.
 */
double Matrix::dot(const Matrix &matrix2) {
    if (col == 1 && matrix2.col == 1 && fil == matrix2.fil) {
        double result = 0.0;

        for (int i = 0; i < fil; i++) {
            result += matrix[i][0] * matrix2.matrix[i][0];
        }

        return result;
    } else {
        throw "Invalid matrices";
    }
}

/**
 * Calculates the cross product between two matrices.
 * @param matrix2 The second matrix for the cross product.
 * @return The cross product result as a new matrix.
 */
Matrix Matrix::cross(const Matrix& matrix2) {
    Matrix result(3, 1);
    if (col == 1 && matrix2.col == 1 && fil == 3 && matrix2.fil == 3) {
        result(0, 0) = this->matrix[1][0] * matrix2.matrix[2][0] - this->matrix[2][0] * matrix2.matrix[1][0];
        result(1, 0) = this->matrix[2][0] * matrix2.matrix[0][0] - this->matrix[0][0] * matrix2.matrix[2][0];
        result(2, 0) = this->matrix[0][0] * matrix2.matrix[1][0] - this->matrix[1][0] * matrix2.matrix[0][0];
        return result;
    } else {
        throw "Invalid values";
    }
}

/**
 * Transposes the matrix. The current matrix remains unchanged.
 * @return A new matrix that is the transpose of the current matrix.
 */
Matrix Matrix::transpose(){
    Matrix result(col, fil);

    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            result(j, i) = matrix[i][j];
        }
    }

    return result;
}
/**
 * Creates an identity matrix of the specified size.
 *
 * @param n The size of the identity matrix (number of rows and columns).
 * @return The identity matrix of size n x n.
 */
Matrix Matrix::identity(int n) {
    Matrix result(n,n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                result(i, j) = 1.0;
            } else {
                result(i, j) = 0.0;
            }
        }
    }

    return result;
}

/**
 * Calculates the inverse of a square matrix.
 *
 * @return The inverse of the matrix.
 * @throws An exception if the matrix is not square or is not invertible.
 */
Matrix Matrix::inverse() {
    if (fil != col) {
        throw "Matrix is not square, inverse does not exist.";
    }

    int n = fil;
    Matrix augmented(n, 2 * n);


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented(i, j) = matrix[i][j];
        }
        augmented(i, i + n) = 1.0;
    }

    for (int i = 0; i < n; i++) {
        if (augmented(i, i) == 0.0) {
            int swap_row = i + 1;
            while (swap_row < n && augmented(swap_row, i) == 0.0) {
                swap_row++;
            }
            if (swap_row == n) {
                throw "Matrix is not invertible.";
            }
            for (int j = 0; j < 2 * n; j++) {
                double temp = augmented(i, j);
                augmented(i, j) = augmented(swap_row, j);
                augmented(swap_row, j) = temp;
            }
        }

        double pivot = augmented(i, i);

        for (int j = 0; j < 2 * n; j++) {
            augmented(i, j) /= pivot;
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmented(k, i);
                for (int j = 0; j < 2 * n; j++) {
                    augmented(k, j) -= factor * augmented(i, j);
                }
            }
        }
    }

    Matrix inverse(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inverse(i, j) = augmented(i, j + n);
        }
    }

    return inverse;
}



