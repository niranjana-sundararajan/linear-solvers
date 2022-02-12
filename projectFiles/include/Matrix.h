#pragma once
#include<memory>

template <class T>
class Matrix
{
public:

	// constructor for allocating the memory
	Matrix(int rows, int cols, bool preallocate);

	// constructor when the memory has already been allocated
	Matrix(int rows, int cols, std::shared_ptr<T[]> values_ptr);

	// destructor
	virtual ~Matrix();

	// Print out the values in our matrix
	void printValues();
	virtual void printMatrix();

	// Matrix operations with values
	std::weak_ptr<T[]> getValues_ptr();

	// Swap rows a and b in the matrix
	void swap_row(int a, int b);

	// Overloading operator += and -= 
	// This will overwrite the matrix
	Matrix<T>& operator+=(const Matrix<T>& mat_left);
	Matrix<T>& operator-=(const Matrix<T>& mat_left);

	// Matrix x Matrix multiplication
	void matMatMult(Matrix<T>& mat_left, Matrix<T>& output);

	// Matrix x Vector multiplication
	void matVecMult(Matrix<T>& vector, Matrix<T>& output);

	// Matrix x Scalar multiplication
	void matScalarMult(T number);

	// Calculating matrix properties
	void transpose();
	double L2norm();
	int sparsity();

	// Decomposition of matrix into upper and(or) lower triangular matrices	
	void LUDecomposition_withPP(Matrix<T>& permutation_mat, Matrix<T>& lower_triangular_mat);
	void choleskyFactorisation(Matrix<T>& output, bool symmetric);

	// The LU Decomposition without Partial pivoting is not stable, therefore we implemented with partial pivotting
	// The decomposition below is for analysis only. Therefore it should be protected method, so only subclasses have access to it
	// However our test class (where we are testing the LU decomp is not a subclass at the moment and if we make it protected we will not be able to 
	// access it  
	void LUDecomposition(Matrix<T>& lower_triangular_mat);


	// Initiate class members
	std::shared_ptr<T[]> values{ nullptr };
	int rows = -1;
	int cols = -1;

protected:
	int size_of_values = -1;
	bool preallocated = false;

private:

};
