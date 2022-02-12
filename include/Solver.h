#pragma once
#include "Matrix.h"
#include "CSRMatrix.h"

template<class T>
class Solver
{

public:
    // -------- CLASS MEMBERS------------------------------------------------------------------------

    // max number of iterations set to default
    int iterations = 100000;


    // count the number of iterations was done durig execution fo a solver
    // this number will be overwritten with every new solver run
    // consideration: consider making this as an output of solver execution
    int numiterationsDone = -1;

    // -------- CONSTRUCTORS--------------------------------------------------------------------------

    Solver();
    Solver(int maxIterations);

    //------------ACCESSORY FUNCTIONS----------------------------------------------------------------

    // Forward substitution
    void forwardSubstitution(Matrix<T>& lower_triangular_mat, Matrix<T>& output, Matrix<T>& vector_mat);

    // Backward Substitution
    void backSubstitution(Matrix<T>& upper_triangular_mat, Matrix<T>& output, Matrix<T>& vector_mat);

    // Error calculation
    double errorRelativeCalc(Matrix<T>& X, Matrix<T>& Xprev);
    double errorAbsoluteCalc(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& vector_x);
    double conditionNumber(Matrix<T>& b_1, Matrix<T>& b_2, Matrix<T>& x_1, Matrix<T>& x_2);

    // Check if matrix is diagonally dominant
    bool diagonal_dominant(Matrix<T>& A);

    //------------SOLVERS------------------------------------------

    // Direct solvers: 

    void gaussElimination(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& output);

    void LUDecomposition_Solver(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T> output);

    void choleskySolver(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& output);

    // Iterative solvers:
    // apply default tolerance, so the user is able to specify tolerance themselves or use default.
    // Did not make it a class member, since it is only applicable to some solvers(iterative)
    void gaussSeidel(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& output, bool relativeError = false, double tolerance = 0.000001);

    void conjugateGradient(Matrix<T>& A, Matrix<T>& B, Matrix<T>& X, double tolerance = 0.000001);

    void jacobiDense(Matrix<T>& A, Matrix<T>& B, Matrix<T>& X, bool relativeError = false, double tolerance = 0.000001);


    //------------SOLVERS------------------------------------------
    // apply default tolerance, so the user is able to specify tolerance themselves or use default
    void jacobiSparse(CSRMatrix<T>& A, Matrix<T>& B, Matrix<T>& X, double tolerance = 0.000001);
    void conjugateGradientSparse(CSRMatrix<T>& A, Matrix<T>& vector_b, Matrix<T>& X, double tolerance = 0.000001);

private:

};
