#include <iostream>
#include<memory>
#include "Matrix.h"

template<class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate) : rows(rows), cols(cols), size_of_values(rows* cols), preallocated(preallocate)
{
    //-----Description------------------------------------------
    // Default Constructor 1  of the class
    //----------------------------------------------------------

    // If memory is preallocated, initialize it to a new pointer
    if (this->preallocated)
    {
        this->values = std::shared_ptr<T[]>(new T[size_of_values]);
    }
}

template<class T>
Matrix<T>::Matrix(int rows, int cols, std::shared_ptr<T[]> values_ptr) : rows(rows), cols(cols), size_of_values(rows* cols), values(values_ptr)
{
    //-----Description------------------------------------------
    // Default Constructor 2  of the class
    //----------------------------------------------------------
}

template<class T>
Matrix<T>::~Matrix() {
    // ------Description-----------------------------
    // Destructor to prevent memory leaks
    //-----------------------------------------------
    // By using the shared pointers we don't need to delete the values
    // they will be deleted by sharedpointer once the last copy of the matrix 
    // is out of scope
}

template<class T>
void Matrix<T>::printValues()
{
    // ------Description-----------------------------
    // Print flattened 1-D array just as values
    //-----------------------------------------------

    std::cout << "Printing out values:" << std::endl;
    for (int i = 0; i < this->size_of_values; i++)
    {
        std::cout << this->values[i] << " ";
    }
    std::cout << std::endl;
}

template<class T>
void Matrix<T>::printMatrix()
{
    // ----Description-----------------------------------
    // Print array values as defined but
    // in the form of a matrix
    // ---------------------------------------------------

    std::cerr << "Printing out matrix:" << std::endl;
    // Row Major Ordering
    for (int j = 0; j < this->rows; j++)
    {
        for (int i = 0; i < this->cols; i++)
        {
            std::cerr << this->values[i + j * this->cols] << " ";
        }
        std::cerr << std::endl;
    }
}

template<class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_left, Matrix<T>& output) {
    // ----Description--------------------------------
    // Computes Matrix - Matrix Multiplication
    // output = mat_left * this
    //-----------------------------------------------

    // ---- Note -------------------------------------
    // Assumption is that output has correct format and
    // it has memory allocated
    // -----------------------------------------------

    // Check dimension match
    //-------------------------------------------------
    // if output(r3*c3) = mat_left(r1*c1) * this(r2*c2)
    // then c3 == c2
    // and r3 == r1 and r2 == c1
    //------------------------------------------------- 


    // For mat_left
    if (this->cols != output.cols || this->rows != mat_left.cols)
    {
        std::cerr << "Error!!Input dimensions dont match!!!" << std::endl;
        return;
    }
    if (mat_left.rows != output.rows) {
        std::cerr << "Error!!Input dimensions dont match!!!" << std::endl;
        return;
    }

    // Check that the output matrix has memory allocated
    if (output.values != nullptr)
    {
        // If no memory is allocated, allocate it explicity
        output.values = std::shared_ptr<T[]>(new T[mat_left.rows * this->cols]);
        output.preallocated = true;
    }

    // Set values to zero beforehand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    // Calculating Matrix-Matrix Multiplication
    // Looping through the rows of mat_left
    for (int i = 0; i < mat_left.rows; i++)
    {
        // For each element in the column of this matrix (ie. RHS matrix)
        for (int j = 0; j < this->cols; j++)
        {
            // For each element in the row of mat_left
            for (int k = 0; k < mat_left.cols; k++)
            {
                output.values[i * this->cols + j] += mat_left.values[i * mat_left.cols + k] * this->values[k * this->cols + j];
            }
        }
    }
}

template<class T>
void Matrix<T>::matVecMult(Matrix<T>& vector, Matrix<T>& output)
{
    // ----Description--------------------------------
    // Computes Matrix - Vector Multiplication
    // output = this matrix * vector
    //-----------------------------------------------

    // ---- Note -------------------------------------
    // Assumption is that output has correct format and
    // it has memory allocated. 
    // -----------------------------------------------

    // Check dimension match
    //-------------------------------------------------
    // if output(r3*c3) = this(r1*c1) * vector(r2*c2)
    // then c3 == 1 (output is vector)
    // and r3 == r1 and r2 == c1
    //------------------------------------------------- 

    if (this->cols != vector.rows || output.cols != 1)
    {
        std::cerr << "Error!!Input dimensions dont match!!!" << std::endl;
        return;
    }

    // Check that the output matrix has memory allocated
    if (output.values != nullptr)
    {
        // Check the dimensions match
        if (output.rows != this->rows)
        {
            std::cerr << "Error!!Input dimensions dont match!!!" << std::endl;
            return;
        }
    }
    // If no memory is allocated, allocate it explicity
    else
    {
        output.values = std::shared_ptr<T[]>(new T[this->rows]);
        output.preallocated = true;
    }

    // Set values to zero beforehand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    // Calculating Matrix-Vector Multiplication
    // Looping through the rows of this matrix (LHS matrix)
    for (int j = 0; j < this->rows; j++)
    {
        // For each element in vector
        for (int k = 0; k < this->cols; k++)
        {
            output.values[j] += this->values[k + j * this->cols] * vector.values[k];
        }
    }
}

template<class T>
std::weak_ptr<T[]> Matrix<T>::getValues_ptr()
{
    // ----Description--------------------------------
    // Return a weak pointer to the values directly
    //-----------------------------------------------
    // Note: weak pointer does not grant ownership.
    // If used then beware whether the matrix has not been
    // deleted during execution
    //------------------------------------------------- 

    std::weak_ptr<T[]> weak_ptr_to_values(this->values);
    return weak_ptr_to_values;
}


template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& mat_left)
{
    // ----Description--------------------------------
    // Calculate Matrix - Matrix Addition
    //-----------------------------------------------
    // Note: The mat_left is added to the current matrix
    // so it will overwrite it. If the original matrix is 
    // required later, the user should make a hard copy
    //------------------------------------------------- 

    // Check the dimensions
    if (mat_left.cols != this->cols || mat_left.rows != this->rows)
    {
        std::cout << "The matrix dimensions are not the same" << std::endl;
    }

    // loop over the elements and add individual values
    for (int i = 0; i < this->size_of_values; i++)
    {
        this->values[i] = this->values[i] + mat_left.values[i];

    }
    return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& mat_left)
{
    // ----Description--------------------------------
    // Calculate Matrix - Matrix Subtraction
    //-----------------------------------------------
    // Note: The mat_left is subtracted from the current matrix
    // so it will overwrite it. If the original matrix is 
    // required later, the user should make a hard copy
    //------------------------------------------------- 
    // 
    // check dimensions
    if (mat_left.cols != this->cols || mat_left.rows != this->rows)
    {
        std::cout << "The matrix dimensions are not the same" << std::endl;
    }
    // loop over the elements and add
    for (int i = 0; i < this->size_of_values; i++)
    {
        this->values[i] = this->values[i] - mat_left.values[i];

    }
    return *this;
}


// Matrix x Scalar multiplication
template<class T>
void Matrix<T>::matScalarMult(T number)
{
    // ----Description--------------------------------
    // Multiply each element of a matrix by a scalar
    // NOTE: 
    // -THIS WILL OVERWRITE THE MATRIX
    // -the scalar has to be of the same type as matrix
    //-----------------------------------------------

    for (int i = 0; i < this->cols * this->rows; i++)
    {
        this->values[i] *= number;
    }
}

template<class T>
void Matrix<T>::swap_row(int a, int b)
{
    // ----Description--------------------------------
    // Swap row a and row b in the matrix
    //-----------------------------------------------
    // Note: This operation will overwrite the current matrix.
    // If the original matrix is required later, 
    // the user should make a hard copy
    //------------------------------------------------- 
    for (int i = 0; i < this->cols; i++)
    {
        //take copy of the element
        T element = this->values[a * this->cols + i];
        this->values[a * this->cols + i] = this->values[b * this->cols + i];
        this->values[b * this->cols + i] = element;
    }
}

template<class T>
void Matrix<T>::LUDecomposition_withPP(Matrix<T>& permutation_mat, Matrix<T>& lower_triangular_mat)
{
    // ------------------------------------------------------------------------
    // 1. LU Decomposition with Partial Pivoting------------------------------- 
    // ------------------------------------------------------------------------
    //  This method splits the matrix A in a lower and upper triangular
    //  matrices and a permutation matrix as result of partial pivoting (A = PLU). 
    // 
    //  Input : Memory allocated for permutation matrix and lower triangular matrix
    //  Output : Void
    //
    //  NOTE: THIS WILL OVERWRITE THE MATRIX A with an upper triangular matrix
    //  Please make a copy of the matrix before running the algorithm if 
    //  matrix A is required after the solution for x is found
    //-------------------------------------------------------------------------
    // Make the permutation matrix into identity matrix
    for (int i = 0; i < this->cols; i++)
    {
        for (int j = 0; j < this->rows; j++)
        {
            if (i == j)
            {
                permutation_mat.values[i * this->cols + j] = 1;
                lower_triangular_mat.values[i * this->cols + j] = 0;
            }
            else
            {
                permutation_mat.values[i * this->cols + j] = 0;
                lower_triangular_mat.values[i * this->cols + j] = 0;
            }
        }
    }
    // For pivotting find the row with maximum value in column k on or below the diagonal
    // we don't need to check the last row
    for (int k = 0; k < this->rows - 1; k++)
    {
        // set the current row to max and iterate through the lower rows to find
        // a larger absolute number - if exists
        int kmax = k;
        double currentmax = pow(this->values[kmax * this->cols + k], 2);
        for (int i = k + 1; i < this->rows; i++)
        {
            double element_sq = pow(this->values[i * this->cols + k], 2);
            if (element_sq > currentmax)
            {
                kmax = i;
                currentmax = element_sq;
            }
        }
        // swap the kmax row in place of k in matrix A or keep as is
        if (kmax != k)
        {
            this->swap_row(kmax, k);
            permutation_mat.swap_row(kmax, k);
            lower_triangular_mat.swap_row(kmax, k);
        }
        // loop over subsequent rows in the same column to turn them to zero
        // start from off diagonals 
        for (int i = k + 1; i < this->rows; i++)
        {
            // define the multiplication factor
            // for row major
            double s = this->values[k + i * this->cols] / this->values[k * this->cols + k];

            // loop over the elements in rows requiring update
            // we start from k on diagonal, because we have already converted all elements before k to 0
            for (int j = k; j < this->rows; j++)
            {
                // this iterates through the columns that are stored next to each other 
                // so it satisfies the row major structure for speed. For column major structure calculating 
                // Upper triangulat would be more efficient
                this->values[i * this->cols + j] = this->values[i * this->cols + j] -
                    s * this->values[k * this->cols + j];
            }

            //set the multiplication factor in lower traingular matrix
            lower_triangular_mat.values[i * this->cols + k] = s;
        }
    }
    for (int i = 0; i < this->cols; i++)
    {
        lower_triangular_mat.values[i * this->cols + i] += 1;
    }
    permutation_mat.transpose();
}
template<class T>
void Matrix<T>::LUDecomposition(Matrix<T>& lower_triangular_mat)
{
    // ------------------------------------------------------------------------
    // 1. LU Decomposition without Partial Pivotting--------------------------- 
    // ------------------------------------------------------------------------
    //  This method splits the matrix A in a lower and upper triangular
    //  matrices (A = LU). (No partial pivoting is implemented)
    //  STRONG RECOMENDATION: USE LU DECOMPOSITION WITH PARTIAL PIVOTING AS IT 
    //  IS MORE STABLE.
    //  
    //  Input : Memory allocated for permutation matrix and lower triangular matrix
    //  Output : Void
    //
    //  NOTE: THIS WILL OVERWRITE THE MATRIX A with an upper triangular matrix
    //  Please make a copy of the matrix before running the algorithm if 
    //  matrix A is required after the solution for x is found
    //-------------------------------------------------------------------------

    for (int i = 0; i < this->cols; i++)
    {
        for (int j = 0; j < this->rows; j++)
        {
            lower_triangular_mat.values[i * this->cols + j] = 0;
        }
    }

    // we don't need to check the last row
    for (int k = 0; k < this->rows - 1; k++)
    {

        // loop over subsequent rows in the same column to turn them to zero
        // start from off diagonals 
        for (int i = k + 1; i < this->rows; i++)
        {
            // define the multiplication factor
            // for row major
            double s = this->values[k + i * this->cols] / this->values[k * this->cols + k];

            // loop over the elements in rows requiring update
            // we start from k on diagonal, because we have already converted all elements before k to 0
            for (int j = k; j < this->rows; j++)
            {
                // this iterates through the columns that are stored next to each other 
                // so it satisfies the row major structure for speed. For column major structure calculating 
                // Upper triangulat would be more efficient
                this->values[i * this->cols + j] = this->values[i * this->cols + j] -
                    s * this->values[k * this->cols + j];
            }

            //set the multiplication factor in lower traingular matrix
            lower_triangular_mat.values[i * this->cols + k] = s;
        }
    }
    for (int i = 0; i < this->cols; i++)
    {
        lower_triangular_mat.values[i * this->cols + i] += 1;
    }
}


template<class T>
void Matrix<T>::transpose()
{
    // ------------------------------------------------------------------------
   // 1. Matrix Transpose ---------------------------------------------------- 
   // ------------------------------------------------------------------------
   //  This method transposes the current matrix.
   // 
   //  Input : Void
   //  Output : Void
   //
   //  NOTE: This will overwrite the matrix with its transpose!!!!!!
   //  Please make a copy of the matrix before running the algorithm if 
   //  matrix A is required after the solution for x is found
   //-------------------------------------------------------------------------
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = i + 1; j < this->cols; j++)
        {
            if (i != j)
            {
                // assuming row major
                T copyelement_1 = this->values[j * this->cols + i];
                this->values[j * this->cols + i] = this->values[i * this->cols + j];
                this->values[i * this->cols + j] = copyelement_1;
            }
        }
    }
}

template <class T>
void Matrix<T>::choleskyFactorisation(Matrix<T>& output, bool symmetric)
{

    // ------------------------------------------------------------------------
    // 1. Solver using Cholesky factorisation for Symmetric Matrix ------------ 
    // ------------------------------------------------------------------------
    //  This method uses Cholesky decomposition of matrix A into upper and lower 
    //  triangular matrices, where the upper and lower matrices are transposes 
    //  of each other. For that reason only the lower triangular matrix is stored
    //  It is also required that the matrix A is symmetric and positive 
    //  definite, i.e. A = A_transpose and it has positive eigenvalues
    // 
    //  Input : Matrix
    //  Output : Void
    //-------------------------------------------------------------------------
    if (!symmetric)
    {
        std::cout << "Matrix is not symmetric definite!" << std::endl;
        std::cout << "It cannot be decomposed using Cholesky factorisation" << std::endl;
        return;
    }

    if (output.values == nullptr)
    {
        for (int i = 0; i < output.rows * output.cols; i++)
        {
            output.values[i] = 0;
        }
    }
    // ----------------------------------------------------
    // iterate through the rows

    for (int i = 0; i < this->rows; i++)
    {
        // for values before the diagonal 
        for (int k = 0; k < i; k++)
        {
            // store the A matrix value so we can suubtract from it 
            T original = this->values[i * this->cols + k];
            // for loop for row major
            for (int j = 0; j < k; j++)
            {
                original -= output.values[i * this->cols + j] * output.values[k * this->cols + j];
            }

            output.values[i * this->cols + k] = original / output.values[k * this->cols + k];
        }
        T original = this->values[i * this->cols + i];
        for (int j = 0; j < i; j++)
        {
            original -= output.values[i * this->cols + j] * output.values[i * this->cols + j];
        }
        output.values[i * this->cols + i] = sqrt(original);
    }
    // ----------------------------------------------------


}

template <class T>
double Matrix<T>::L2norm()
{
    // ------------------------------------------------------------------------
     // 1. Matrix Transpose ---------------------------------------------------- 
     // ------------------------------------------------------------------------
     //  This method calculates the L2 norm of the matrix/vector.
     //  L2 norm = sqrt(sum(A_ij * Aij))
     // 
     //  Input : Void
     //  Output : Void
     //-------------------------------------------------------------------------

    double sum = 0;
    double l2norm = 0;
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            sum += pow(this->values[i * this->cols + j], 2);
        }
    }
    return sqrt(sum);
}

template <class T>
int Matrix<T>::sparsity()
{
    // calculate number of non zero values in a Matrix
    int nnzv = 0;
    double l2norm = 0;
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            if (this->values[i * this->cols + j] != 0)
            {
                nnzv += 1;
            }
        }
    }
    return nnzv;
}

template class Matrix<double>;
template class Matrix<int>;