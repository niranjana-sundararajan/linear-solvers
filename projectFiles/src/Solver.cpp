#pragma once
#include <iostream>
#include <algorithm>
#include <cmath>
#include "CSRMatrix.h"
#include "Solver.h"

using namespace std;

// --------------------------------------------------------------------------
// ----------------- CONSTRUCTOR FOR SOLVER CLASS---------------------------
// --------------------------------------------------------------------------

template <class T>
Solver<T>::Solver(int maxIterations)
{
    // -- Constructor I --------------------------------------------------------
    // 
    // Constructor to define  the maximum number of iterations.
    // This is specific to the iterative solvers.
    // 
    //-------------------------------------------------------------------------


    this->iterations = maxIterations;

    ;
}

template <class T>
Solver<T>::Solver()
{
    // -- Constructor II --------------------------------------------------------
    // 
    // Constructor for direct methods and use of default values for tolerance and 
    // maximum number of iterations
    // 
    //---------------------------------------------------------------------------

}


// -------------------------------------------------------------------------------
// ---- ACCESSORY FUNCTIONS REQUIRED IN SOLVER CLASS------------------------------
// -------------------------------------------------------------------------------

template<class T>
void Solver<T>::backSubstitution(Matrix<T>& upper_triangular_mat, Matrix<T>& output, Matrix<T>& vector_mat)
{
    // -- Description --------------------------------------------------------
    // Function that performs back substitution on an upper triangular matrix 
    // and modifies values at the output pointer. Solves the eqution U*o = v
    // where U is upper triangular, v is vector_mat and o is output. 
    // The solution is stored in the output matrix.
    // 
    // Input : Matrix, Matrix, Matrix
    // Output : void
    //-------------------------------------------------------------------------
    // populate the output with zero if the matrix is empty
    if (output.values == nullptr)
    {
        for (int i = 0; i < output.rows * output.cols; i++)
        {
            output.values[i] = 0;
        }
    }

    int nrows = upper_triangular_mat.rows;

    // iterate over the rows starting from the bottom and working way up
    for (int i = nrows - 1; i >= 0; i--) {
        double sub_element = 0;

        for (int j = i + 1; j < nrows; j++) {
            sub_element += upper_triangular_mat.values[i * nrows + j] * output.values[j];
        }
        output.values[i] = (vector_mat.values[i] - sub_element) / upper_triangular_mat.values[i * nrows + i];
    }
}

template<class T>
void Solver<T>::forwardSubstitution(Matrix<T>& lower_triangular_mat, Matrix<T>& output, Matrix<T>& vector_mat) {
    // -- Description --------------------------------------------------------
    // Function that performs forward substitution on a lower triangular matrix 
    // and modifies values at the output pointer. Solves the eqution L*o = v
    // where L is lower triangular, v is vector_mat and o is output.
    // The solution is stored in the output matrix.
    // 
    // Input : Matrix, Matrix, Matrix
    // Output : Void
    //-------------------------------------------------------------------------
    // 
    // Populate output with zero's if it's empty
    if (output.values == nullptr)
    {
        for (int i = 0; i < output.rows * output.cols; i++)
        {
            output.values[i] = 0;
        }
    }
    int nrows = lower_triangular_mat.rows;

    for (int i = 0; i < nrows; i++) {
        double sub_element = 0;
        for (int j = 0; j < i; j++) {
            sub_element += lower_triangular_mat.values[i * nrows + j] * output.values[j];
        }
        output.values[i] = (vector_mat.values[i] - sub_element) / lower_triangular_mat.values[i * nrows + i];
    }
}

template<class T>
bool Solver<T>::diagonal_dominant(Matrix<T>& A)
{
    // -- Description --------------------------------------------------------
    // Check if a matrix is diagonally dominant.
    // This is one of the conditions for application of Gauss-Siedel method.
    // The matrix must be diagonally dominant or symmetric definite positive
    // This condition checks for all rows:
    //                  abs(a[i][i])> abs(a[i][:] - a[i][i])
    // 
    // Input : Matrix
    // Output : Bool
    //-------------------------------------------------------------------------
   // Initiate diagonal
    T A_diag = 0;
    //iterate through each row
    for (int i = 0; i < A.rows; i++)
    {
        A_diag = abs(A.values[i * A.rows + i]);
        T sum = 0;
        for (int j = 0; j < A.cols; j++)
        {
            if (j != i)
            {
                // sum up over all values in the row except for diagonal
                sum += abs(A.values[i * A.rows + j]);
            }
        }
        // if  a single row is not diagonally dominant break and return false
        if (A_diag <= sum)
        {
            return false;
        }
    }
    return true;
}

template <class T>
double Solver<T>::conditionNumber(Matrix<T>& b_1, Matrix<T>& b_2, Matrix<T>& x_1, Matrix<T>& x_2)
{
    // ------------------------------------------------------------------------
     // 1.Condition number => || Ax -b || ------------------------------------
     // 
     // Condition number is an indication of how ill-posed a linear system is
     // ------------------------------------------------------------------------

     // relative error of the input
    double relativeError_b = errorRelativeCalc(b_1, b_2);

    // relative error of the output
    double relativeError_x = errorRelativeCalc(x_1, x_2);

    double conditionNumber = relativeError_x / relativeError_b;

    return conditionNumber;
}


// Error Calcuration -------------------------
template <class T>
double Solver<T>::errorAbsoluteCalc(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& vector_x)
{
    // ------------------------------------------------------------------------
     // 1.Error Caluclation => || Ax -b || ------------------------------------
     // 
     // Calculating the absolute error of derived solution 
     // ------------------------------------------------------------------------

    // create temp vector for residuals
    auto* output = new Matrix<T>(A.cols, 1, true);

    // calculate Ax -b
    A.matVecMult(vector_x, *output);
    *output -= vector_b;

    // calculate norm of error
    double error_maginude = output->L2norm();

    // delete the pointer and array
    delete[] output;

    return error_maginude;
}

template <class T>
double Solver<T>::errorRelativeCalc(Matrix<T>& X, Matrix<T>& Xprev)
{
    // ------------------------------------------------------------------------
     // 1.Error Caluclation => || X - Xprev||/||X||-----------------------------
     // 
     // Calculating the relative error of derived solution, where X is at the new
     // iteration and Xprev is previous iteration. 
     // This can also be used to calculate relative error to exact solution, 
     // where X is the exact solution and Xprev is the derived solution
     // 
     // Input : Matrix X, Matrix Xprev
     // Output: long double 
     // ------------------------------------------------------------------------

    // intiialising the variables
    double Xdenominator = 0;
    double Xresidual = 0;
    double RelativeError_X = 0;


    // Loop for squaring each values and adding
    for (int i = 0; i < X.rows; i++)
    {
        // squaring previous X and adding
        Xdenominator += std::pow(Xprev.values[i], 2);
        // calcurating difference of new X and previous X and squaring and adding them.
        Xresidual += std::pow((X.values[i] - Xprev.values[i]), 2);

        // check that the denominator is non- zero
        if (Xdenominator == 0)
        {
            cerr << "Error!! the Input X is zero vector" << endl;
            return -1;
        }
    }

    // Calcurating Relative error
    RelativeError_X = sqrt(Xresidual) / sqrt(Xdenominator);

    return RelativeError_X;
}

// --------------------------------------------------------------------------
// --------------------SOLVER FUNCTIONS -------------------------------------
// --------------------------------------------------------------------------

template<class T>
void Solver<T>::gaussElimination(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& output)
{
    // ------------------------------------------------------------------------
    // 1. Gauss Elimination --------------------------------------------------- 
    // ------------------------------------------------------------------------
    //  This method uses gaussian elimination with partial pivoting
    //  to solver the linear system A * x = b
    // 
    //  Input : Matrix
    //  Output : Void
    //-------------------------------------------------------------------------

    // ---- Partial Pivoting ----------------------
    // loop down the column to find the max column and 
    // swap the rows such that the greatest value is 
    // at the top. The loop runs for all rows below the
    // row under consideration for a swap

    for (int i = 0; i < A.rows - 1; i++)
    {
        // set index of max row as current row index
        int max_index = i;
        // Therefore. max row value is current value
        T max_value = A.values[i * A.cols + i];

        // Loop through the column below the current element
        for (int j = i + 1; j < A.rows; j++)
        {
            // loop to find row below the pivot with 
            // maximum absolute value
            if (fabs(A.values[j * A.cols + i]) > max_value) {
                max_value = A.values[j * A.cols + i];
                max_index = j;
            }
        }
        // Swapping rows using accessory function
        if (max_index != i) {
            A.swap_row(i, max_index);
            vector_b.swap_row(i, max_index);
        }
        // Updating values of matrices
        for (int m = i + 1; m < A.rows; m++) {
            T row_coeff = A.values[m * A.cols + i] / A.values[i * A.cols + i];
            for (int n = i; n < A.rows; n++) {
                A.values[m * A.cols + n] -= row_coeff * A.values[i * A.cols + n];
            }
            vector_b.values[m] = vector_b.values[m] - row_coeff * vector_b.values[i];
        }
        // do backward substitution to obtain the solution, X
        backSubstitution(A, output, vector_b);

    }
}

template<class T>
void Solver<T>::LUDecomposition_Solver(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T> output)
{
    // ------------------------------------------------------------------------
    // 1. Solver Using LU Decomposition of Matrix A --------------------------- 
    // ------------------------------------------------------------------------
    //  This method uses the decomposition of matrix A into upper and lower 
    //  triangular matrices to find x for the linear system A * x = b, by splitting 
    //  it into two related linear systems
    // 
    //  Ax=b => (PLU)x=b => L(Ux) =P_inverse*b
    //  => Ly=P_inverse*b where Ux = y;
    // 
    //  Input : Matrix
    //  Output : Void
    //-------------------------------------------------------------------------

    // Check the dimension
    if (vector_b.rows != A.rows || output.rows != vector_b.rows || output.rows != A.cols)
    {
        std::cerr << "Dimentions do not match" << std::endl;
        return;
    }
    if (A.rows != A.cols)
    {
        std::cerr << "The matrix is not square!" << std::endl;
        return;
    }

    // Decompose the matrix into Permutation, Upper and Lower matrices
    // use unique_ptr so it is deleted after the solver has finished running
    auto* P = new Matrix<T>(A.rows, A.cols, true);
    auto* L = new Matrix<T>(A.rows, A.cols, true);

    A.LUDecomposition_withPP(*P, *L);

    // Using forward substitution find y in the linear system
    // L*y = P_inverse * b
    // overwrite P with its inverse, which for permutation matrix is its transpose
    P->transpose();
    auto* P_inv_b = new Matrix<T>(A.rows, 1, true);
    auto* y = new Matrix<T>(A.rows, 1, true);

    P->matVecMult(vector_b, *P_inv_b);
    // we no longer need Permutation matrix, so we release the memory 
    // during runtime rather than at the end 
    delete[] P;

    forwardSubstitution(*L, *y, *P_inv_b);
    delete[] L;
    delete[] P_inv_b;
    // Using backward substitution find x in the linear system
    // U*x = y (A has been overwritten with U)
    backSubstitution(A, output, *y);

    delete[] y;
}


template <class T>
void Solver<T>::gaussSeidel(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& output, bool relativeError, double tolerance)
{
    // ------------------------------------------------------------------------
    // 1. Solver Using Gauss_Siedel Iteration Method --------------------------- 
    // ------------------------------------------------------------------------
    //  This method iterates  the decomposition of matrix A into upper and lower 
    //  triangular matrices to find x for the linear system A * x = b, by splitting 
    //  it into two related linear systems
    // Warning! The matrix should be a symmetric positive definite (meaning all 
    // eigenvalues are positive) or it must be diagonally dominant.
    // Please run function diagonal_dominant to verify if the property is unknown
    // 
    //  Input : Matrix A, 
    //          Matrix vector_b, 
    //          Matrix output (initial guess for X or empty)
    //          bool relativeError (false will use absolute error, true will use relative
    //              error to previous iteration. !!!!Note use of relative error will require
    //              extra storage of previous iteration of X). Default = false
    //          double tolerance (Default 0.000001)
    // 
    //  Output : Void
    //-------------------------------------------------------------------------


    // Populate output with zero's 
    for (int i = 0; i < output.rows; i++)
    {
        output.values[i] = 0;
    }

    Matrix<T>* Xprev(nullptr);
    // Create storage for previous iteration of solution and copy the initial guess
    if (relativeError)
    {
        Xprev = new Matrix<T>(A.rows, 1, true);
        for (int i = 0; i < A.rows; i++)
        {
            Xprev->values[i] = output.values[i];
        }
    }

    double omega = 1;
    int k_relax = 10;
    double original = 0;
    double residual = 1;
    double last_residuals = 1;
    double current_residuals = 1;

    // Perform maximum number of iterations as defined by the user 
    // in solver contructor or default value is used
    for (int k = 0; k < this->iterations; k++)
    {
        k_relax = k;
        // for each row in matrix A
        for (int i = 0; i < A.rows; i++)
        {
            double diag = A.values[i + i * A.cols];
            original = vector_b.values[i];
            //cerr << "vector_b" << vector_b.values[i] << endl;
            for (int j = 0; j < A.cols; j++)
            {
                if (j != i)
                {
                    original -= A.values[j + i * A.cols] * output.values[j] + (1 - omega) * output.values[j];

                }
            }
            output.values[i] = original * omega / diag;

        }
        // omega calculation for adjusted relaxation. No one constan omega fits all systems
        // hence we need to calculate it during runtime
        // the overelaxation allows for faster convergence
        if (k == k_relax)
        {
            last_residuals = residual;
        }
        if (k == k_relax + 1)
        {
            current_residuals = residual;
            omega = 2 / (1 + sqrt(1 - current_residuals / last_residuals));
        }
        // Calculate every 5th iteraton, this will help reduce the amount of work
        // during each iteration. This may mean we do up to 5 extra unneccessary iterations
        // but calculating error is costly. 

        if (k % 1 == 0)
        {
            if (relativeError)
            {
                // Calculate relative erro
                residual = errorRelativeCalc(output, *Xprev);
                if (k % 1000 == 0)
                {
                    std::cerr << "The iteration is: " << k << std::endl;
                    std::cerr << "The error magnitude is: " << residual << std::endl;
                }
            }
            else
            {
                residual = errorAbsoluteCalc(A, vector_b, output);

                if (k % 1000 == 0)
                {
                    std::cerr << "The iteration is: " << k << std::endl;
                    std::cerr << "The error magnitude is: " << residual << std::endl;
                }

            }

            //--- Check the convergence is less or more than tolerence --
            if (residual < tolerance)
            {
                std::cerr << "The solver has converged!" << std::endl;
                std::cerr << "The error magnitude is: " << residual << std::endl;
                this->numiterationsDone = k;
                delete[] Xprev;
                return;
            }
            else
            {
                if (relativeError)
                {
                    // copy the values only if required to do another iteration to save time
                    for (int i = 0; i < output.rows; i++)
                    {
                        Xprev->values[i] = output.values[i];

                    }
                }
            }
        }
    }
    cerr << "The solver did not converge, the maximum number of iterations is reached" << std::endl;
    cerr << "The error magnitude is:" << residual << std::endl;
    this->numiterationsDone = this->iterations;
    delete[] Xprev;
    return;
}




// Cholesky ------------------------------- 
template <class T>
void Solver<T>::choleskySolver(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& output)
{
    // ------------------------------------------------------------------------
    // 1. Solver Using Cholesky Decomposition of Matrix A --------------------- 
    // ------------------------------------------------------------------------
    //  This method uses Cholesky decomposition of matrix A into upper and lower 
    //  triangular matrices,to find x for the linear system A * x = b, 
    //  by splitting it into two related linear systems:
    // 
    //  Ax=b => (L*L_trans)x = b => L(L_trans*x) =b
    //  => Ly=b where L_trans*x = y;
    // 
    //  Input : Matrix
    //  Output : Void
    //-------------------------------------------------------------------------

    cerr << "Please note, this solver is for symmetric definite matrices only!!" << std::endl;

    // perform decomposition to find the lower triangular matrix L
    auto* L = new Matrix<T>(A.rows, A.cols, true);
    A.choleskyFactorisation(*L, true);
    // Using forward substitution find y in the linear system
    // L*y = b
    auto* y = new Matrix<T>(A.rows, 1, true);

    forwardSubstitution(*L, *y, vector_b);

    // Using backward substitution find x in the linear system
    // L_trans*x = y
    L->transpose();

    backSubstitution(*L, output, *y);

    delete[] L;
    delete[] y;

}

// 6. Conjugate Gradient --------------------- 
template <class T>
void Solver<T>::conjugateGradient(Matrix<T>& A, Matrix<T>& vector_b, Matrix<T>& X, double tolerance)
{
    // First guess at solution is 0
    for (int i = 0; i < X.rows; i++)
    {
        X.values[i] = 0;

    }
    //  residuals & residuals copy - storing residuals
    //  we will require both sets of errors from previous and current iterations
    auto* residuals_0 = new Matrix<T>(A.cols, 1, true);
    auto* residuals_1 = new Matrix<T>(A.cols, 1, true);

    // gradient_direction - current direction of gradient progression
    //  we will require both sets of gradient descend from previous and current iterations
    auto* gradient_direction_0 = new Matrix<T>(A.cols, 1, true);
    auto* gradient_direction_1 = new Matrix<T>(A.cols, 1, true);

    // Initialize the error term that will be used in stopping condition
    double residual_error = 1;

    // Initialize loop counters
    int iteration_counter = 0;

    // Calculate the first error terms from initial guess at X
    // Leave the below multiplication in, in case user passes in x vector that is other than zero
    A.matVecMult(X, *residuals_0);
    *(residuals_0) -= vector_b;

    // Calculated the first gradient 
    for (int i = 0; i < A.rows; i++)
    {
        // put the residuals - b in this loop
        gradient_direction_0->values[i] = -residuals_0->values[i];
    }

    // Set stopping conditions: max number of iterations or tolerance
    for (int iteration_counter = 0; iteration_counter < this->iterations; iteration_counter++)
    {
        // initiate the coefficients outside the loop
        long double alpha_0 = 0;
        long double alpha_1 = 0;
        long double beta = 0;

        // calculate A*gradient and alpha 
        for (int i = 0; i < A.rows; i++)
        {
            double Ap = 0;
            for (int j = 0; j < A.cols; j++)
            {
                Ap += A.values[i * A.cols + j] * gradient_direction_0->values[j];

            }

            alpha_0 += pow(residuals_0->values[i], 2);
            alpha_1 += Ap * gradient_direction_0->values[i];
        }
        alpha_0 = alpha_0 / alpha_1;

        // Calculate the next estimate of X 
        for (int j = 0; j < vector_b.rows; j++)
        {
            X.values[j] = X.values[j] + (alpha_0 * gradient_direction_0->values[j]);

        }

        // Calculate new residuals using alpha coefficient and A*gradient_descend
        for (int i = 0; i < A.rows; i++)
        {
            double Ap = 0;
            for (int j = 0; j < A.cols; j++)
            {
                Ap += A.values[i * A.cols + j] * gradient_direction_0->values[j];
            }

            residuals_1->values[i] = residuals_0->values[i] + (alpha_0 * Ap);

        }

        //Calculate the beta coefficient
        double r1 = 0;
        double r2 = 0;

        for (int i = 0; i < residuals_0->rows; i++)
        {
            r1 += pow(residuals_0->values[i], 2);
            r2 += pow(residuals_1->values[i], 2);

        }
        beta = r2 / r1;

        // gradient descend update based on new residuals and beta coefficient
        for (int i = 0; i < residuals_0->rows; i++)
        {
            gradient_direction_1->values[i] = -residuals_1->values[i] + beta * gradient_direction_0->values[i];
        }

        // Calculate the findal error magnitude from this iteration
        residual_error = residuals_1->L2norm() / sqrt(A.cols);

        //only update if another iteration is required
        // to save unnecessary calculations
        if (tolerance < residual_error)
        {
            for (int i = 0; i < A.rows; i++)
            {
                gradient_direction_0->values[i] = gradient_direction_1->values[i];
                residuals_0->values[i] = residuals_1->values[i];

            }

            if (iteration_counter % 1000 == 0)
            {
                std::cerr << "The iteration is: " << iteration_counter << std::endl;
                std::cerr << "The error magnitude is: " << residual_error << std::endl;

            }
        }
        else
        {
            std::cerr << "The solver has converged!" << std::endl;
            std::cerr << "The error magnitude is:" << residual_error << std::endl;
            this->numiterationsDone = iteration_counter;

            // delete the arrays no longer in use
            delete[] residuals_0;
            delete[] residuals_1;
            delete[] gradient_direction_1;
            delete[] gradient_direction_0;
            return;

        }
        iteration_counter++;

    }
    cerr << "The solver did not converge, the maximum number of iterations is reached" << std::endl;
    cerr << "The error magnitude is:" << residual_error << std::endl;
    this->numiterationsDone = this->iterations;

    // delete the arrays no longer in use
    delete[] residuals_0;
    delete[] residuals_1;
    delete[] gradient_direction_1;
    delete[] gradient_direction_0;
}




// Jacobi ---------------------------------
template <class T>
void Solver<T>::jacobiDense(Matrix<T>& A, Matrix<T>& B, Matrix<T>& X, bool relativeError, double tolerance)
{

    // set the output to zeros

    for (int i = 0; i < A.cols; i++)
    {
        X.values[i] = 0;
    }
    // Create storage for previous iteration of solution and copy the initial guess
    Matrix<T> Xprev = Matrix<T>(A.rows, 1, true);

    for (int i = 0; i < A.cols; i++)
    {
        Xprev.values[i] = X.values[i];
    }

    // Iteration counter
    int iteration_count = 0;
    double residuals = 10;

    // As long as convergence criteria is not met
    // run the loop
    while (iteration_count <= this->iterations && residuals > tolerance)
    {
        for (int i = 0; i < A.rows; i++)
        {
            double sum = 0;

            for (int j = 0; j < A.cols; j++)
            {
                if (i != j)
                {
                    // ASSUMING ROW MAJOR ORDERING HERE
                    sum += A.values[j + i * A.cols] * X.values[j];

                    //cerr << "sum" << sum << endl;
                }
            }
            // Copy/reassign the value to output
            X.values[i] = (B.values[i] - sum) / (A.values[i + i * A.cols]);

        }

        if (iteration_count % 5 == 0)
        {
            if (relativeError)
            {
                // Calculate relative erro
                residuals = errorRelativeCalc(X, Xprev);

            }
            else
            {
                residuals = errorAbsoluteCalc(A, B, X);
            }
            //--- Check the convergence is less or more than tolerence --
            if (residuals < tolerance)
            {
                std::cerr << "The solver has converged!" << std::endl;
                std::cerr << "The error magnitude is:" << residuals << std::endl;
                this->numiterationsDone = iteration_count;

                return;
            }
        }

        // Copy the values 
        for (int i = 0; i < A.rows; i++)
        {
            Xprev.values[i] = X.values[i];
        }

        iteration_count++;
    }
    // The loop did not break early, so we reached the maximum number of iterations :(
    cerr << "The solver did not converge, the maximum number of iterations is reached" << std::endl;
    cerr << "The error magnitude is:" << residuals << std::endl;
    this->numiterationsDone = iteration_count;

    return;
}


template <class T>
void Solver<T>::jacobiSparse(CSRMatrix<T>& A, Matrix<T>& B, Matrix<T>& X, double tolerance)
{

    //-----Description-----------------------------------------
    // Function that solves a linear system stored in sparse  
    // format using Jacobi iterative method.
    // Input: Matrix<T>, Matrix<T>
    // X : Void
    //--------------------------------------------------------

    // initialize vector to store the previous X values
    double* X_temp = new double[A.rows];

    for (int i = 0; i < A.rows; i++)
    {
        X_temp[i] = 0;
    }

    //set the X to zeros
    for (int i = 0; i < A.rows; i++)
    {
        X.values[i] = 0;
    }
    // initialize sum and counter to zero
    T sum = 0;
    int counter = 0;
    double residuals = 10;

    // An array to store the diagonal values of A and
    T* diag = new T[A.rows];
    // Fill in the array of diagonals
    for (int i = 0; i < A.rows; i++)
    {
        // Jacobi required diagonally dominant therefore the diagonal
        // value will always be present
        for (int j = A.row_position[i]; j < A.row_position[i + 1]; j++)
        {
            if (A.col_index[j] == i)
            {
                diag[i] = A.values[j];
            }
        }
    }
    // Run till convergence criteria is met OR
    //loop count expires (ie reaches 10000)

    while (residuals > tolerance)
    {
        counter++;
        // continue the loop only until the last row
        // since we cant reference to i+1 for the last row
        for (int i = 0; i < A.rows - 1; i++)
        {
            sum = 0;
            int nz_id = 0;
            // T diag_element = 0;
            int nz_start = A.row_position[i];
            int nz_end = A.row_position[i + 1];

            // Loop through the elements in the row
            for (nz_id = nz_start; nz_id < nz_end; nz_id++)
            {
                // For values other than diagonal values on the row
                if (i != A.col_index[nz_id])
                {
                    // for each element in row of A multiply by the corresponding X value
                    sum += A.values[nz_id] * X.values[A.col_index[nz_id]];

                }

            }
            // calculate residual values for next iteration
            // for X in row i
            X_temp[i] = (B.values[i] - sum) / diag[i];

            // store the current values of X for next iteration
            {
                X.values[i] = X_temp[i];

            }
        }
        // Repeat the process for the last row
        int nz_id = 0;
        // start of the last row (recall the last value in row position is the total number of non zeros)
        int nz_start = A.row_position[A.rows - 1];
        // will use in the loop with < size so we go up to and including the last index
        int nz_end = A.nnzs;
        // loop through the last row
        T sum_last = 0;
        for (int i = nz_start; i < nz_end;i++)
        {
            // for all values apart from on diagonal
            if (A.col_index[i] != A.rows - 1)
            {
                // add A*x for each element in the row
                sum_last += A.values[i] * X.values[A.col_index[i]];

            }
        }

        X_temp[A.rows - 1] = (B.values[A.rows - 1] - sum_last) / diag[A.rows - 1];
        // store the calculated X values for the next iteration
        X.values[A.rows - 1] = X_temp[A.rows - 1];

        // Set X values with X_temp values  for (int i = 0; i < A.rows; i++)
        residuals = residualsCSR(A, B, X);
        cerr << counter << endl;
    }
    //Check for Residuals
    //residuals = residualsCSR(A, B, X);
    delete[] X_temp;

    this->numiterationsDone = counter;
}

template <class T>
double residualsCSR(CSRMatrix<T>& A, Matrix<T>& b, Matrix<T>& x)
{
    T value = 0;
    double residual = 0;
    for (int i = 0; i < A.rows; i++)
    {
        // initiate the start and end of the elements in the row
        int nz_id = 0;
        int nz_start = A.row_position[i];
        int nz_end = A.row_position[i + 1];

        value = 0;
        // loop through the elements in the row
        for (int nz_id = nz_start; nz_id < nz_end; nz_id++)
        {
            value += A.values[nz_id] * x.values[A.col_index[nz_id]];
        }
        // calculate (A*x - b)^2 for the row and add it to the otherall resildual
        residual += pow((value - b.values[i]), 2);

    }
    residual = sqrt(residual);

    return residual;
}

// 6. Conjugate Gradient --------------------- 
template <class T>
void Solver<T>::conjugateGradientSparse(CSRMatrix<T>& A, Matrix<T>& vector_b, Matrix<T>& X, double tolerance)
{
    // First guess at solution is 0
    for (int i = 0; i < X.rows; i++)
    {
        X.values[i] = 0;

    }
    //  residuals & residuals copy - storing residuals
    //  we will require both sets of errors from previous and current iterations
    // Have to store as pointers, so that we dont use stack memory in case of large vectors 
    auto* residuals_0 = new Matrix<T>(A.cols, 1, true);
    auto* residuals_1 = new Matrix<T>(A.cols, 1, true);

    // gradient_direction - current direction of gradient progression
    //  we will require both sets of gradient descend from previous and current iterations
    // Have to store as pointers, so that we dont use stack memory in case of large vectors 
    auto* gradient_direction_0 = new Matrix<T>(A.cols, 1, true);
    auto* gradient_direction_1 = new Matrix<T>(A.cols, 1, true);

    // Initialize the error term that will be used in stopping condition
    double residual_error = 1;

    // Initialize loop counters
    int iteration_counter = 0;

    // Calculate the first error terms from initial guess at X
    // Leave the below multiplication in, in case user passes in x vector that is other than zero
    //This uses the overloaded function matVecmult for CSRMatrix due to input parameter being CSRMatrix 
    A.matVecMult(X, *residuals_0);
    *(residuals_0) -= vector_b;

    // Calculated the first gradient 
    for (int i = 0; i < A.rows; i++)
    {
        // put the residuals - b in this loop
        gradient_direction_0->values[i] = -residuals_0->values[i];
    }
    gradient_direction_0->printMatrix();
    // initiate the coefficients outside the loop


    // Set stopping conditions: max number of iterations or tolerance
    for (int iteration_counter = 0; iteration_counter < this->iterations; iteration_counter++)
    {
        long double alpha_0 = 0;
        long double alpha_1 = 0;
        double beta = 0;

        // calculate A*gradient and alpha 
        for (int i = 0; i < A.rows; i++)
        {
            // initiate the start and end of the elements in the row
            //int nz_id = 0;
            int nz_start = A.row_position[i];
            int nz_end = 0;
            //cerr << "i:   " << i << endl;
            // if we are at the last row then we dont go to the next value but to the end of array
            if (i != A.rows - 1)
            {
                nz_end = A.row_position[i + 1];
                // cerr << "nz_end  :  " << nz_end << endl;
            }
            else
            {
                nz_end = A.nnzs;
            }

            double Ap = 0;
            // loop through the elements in the row
            for (int nz_id = nz_start; nz_id < nz_end; nz_id++)
            {
                Ap += A.values[nz_id] * gradient_direction_0->values[A.col_index[nz_id]];
            }
            alpha_0 += pow(residuals_0->values[i], 2);
            alpha_1 += Ap * gradient_direction_0->values[i];
        }

        alpha_0 = alpha_0 / alpha_1;

        // Calculate the next estimate of X 
        // Vector operations no sparse matrix 
        for (int j = 0; j < vector_b.rows; j++)
        {
            X.values[j] = X.values[j] + (alpha_0 * gradient_direction_0->values[j]);

        }

        // Calculate new residuals using alpha coefficient and A*gradient_descend
        for (int i = 0; i < A.rows; i++)
        {
            double Ap = 0;
            // initiate the start and end of the elements in the row
            //int nz_id = 0;
            int nz_start = A.row_position[i];
            int nz_end = 0;
            // if we are at the last row then we dont go to the next value but to the end of array
            if (i != A.rows - 1)
            {
                nz_end = A.row_position[i + 1];
            }
            else
            {
                nz_end = A.nnzs;
            }
            // loop through the elements in the row
            for (int nz_id = nz_start; nz_id < nz_end; nz_id++)
            {
                Ap += A.values[nz_id] * gradient_direction_0->values[A.col_index[nz_id]];
            }

            residuals_1->values[i] = residuals_0->values[i] + (alpha_0 * Ap);

        }

        //Calculate the beta coefficient
        double r1 = 0;
        double r2 = 0;

        for (int i = 0; i < residuals_0->rows; i++)
        {
            r1 += pow(residuals_0->values[i], 2);
            r2 += pow(residuals_1->values[i], 2);

        }
        beta = r2 / r1;

        // gradient descend update based on new residuals and beta coefficient
        for (int i = 0; i < residuals_0->rows; i++)
        {
            gradient_direction_1->values[i] = -residuals_1->values[i] + beta * gradient_direction_0->values[i];
        }

        // Calculate the findal error magnitude from this iteration
        residual_error = residuals_1->L2norm() / sqrt(A.cols);

        cerr << "residual_error" << residual_error << "   " << endl;
        //only update if another iteration is required
        // to save unnecessary calculations
        if (tolerance < residual_error)
        {
            for (int i = 0; i < A.rows; i++)
            {
                gradient_direction_0->values[i] = gradient_direction_1->values[i];
                residuals_0->values[i] = residuals_1->values[i];

            }
        }
        else
        {
            std::cerr << "The solver has converged!" << std::endl;
            std::cerr << "The error magnitude is:" << residual_error << std::endl;
            this->numiterationsDone = iteration_counter;
            delete[] residuals_0;
            delete[] residuals_1;
            delete[] gradient_direction_1;
            delete[] gradient_direction_0;
            return;

        }
        iteration_counter++;

    }
    cerr << "The solver did not converge, the maximum number of iterations is reached" << std::endl;
    cerr << "The error magnitude is:" << residual_error << std::endl;
    this->numiterationsDone = this->iterations;

    // delete the arrays no longer in use
    delete[] residuals_0;
    delete[] residuals_1;
    delete[] gradient_direction_1;
    delete[] gradient_direction_0;
}

// for defining the solvers to help across different compliers
template class Solver<double>;
template class Solver<int>;