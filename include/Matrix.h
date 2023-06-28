#ifndef _MATRIX_
#define _MATRIX_

class Matrix
{
    public:
        Matrix(int fil, int col);
        Matrix(int fil, int col, const double v[], int n);
        Matrix(const Matrix& m);
        ~Matrix();

        int fils();
        int cols();
 
        Matrix& operator=(const Matrix& matrix2);
        Matrix  operator+(const Matrix& matrix2);
        Matrix  operator-(const Matrix& matrix2);
        Matrix  operator*(const Matrix& matrix2);
        Matrix  operator*(const double n);
        Matrix operator/(const double n);
        Matrix operator-();
        double& operator()(const int i, const int j) const;
 
        void print();
        double dot(const Matrix& matrix2);
        Matrix cross(const Matrix& matrix2);
        Matrix transpose();
        static Matrix identity(int n);
        Matrix inverse();

 
    private:
        void initMatrix();
        int fil;
        int col;
        double **matrix;
};

#endif
