#include <iostream>
#include<memory>
#include "Matrix.h"
#include "CSRMatrix.h"

template<class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate) : Matrix<T>(rows, cols, false), nnzs(nnzs)
{
	//--------------------------------------
	//  Constructor 1
	//  with memory allocation
	//_______________________________________
	// If preallocate is true 
	// then the Matrix contructor will allocate
	// memory the size of rows*cols, which is not
	// what we want since  we have a sparse matrix
	// So we overwrite it to false in Matrix constructor 
	// and allocate memory int he CSRMatrix constructor

	this->preallocated = preallocate;
	this->size_of_values = nnzs;

	// allocate memory for sparse matrix if preallocate was set to true
	if (this->preallocated)
	{
		this->row_position = std::shared_ptr<int[]>(new int[this->rows + 1]);
		this->col_index = std::shared_ptr<int[]>(new int[this->nnzs]);
		this->values = std::shared_ptr<T[]>(new T[this->nnzs]);
	}
}

template<class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, std::shared_ptr<T[]> values_ptr, std::shared_ptr<int[]>  row_position, std::shared_ptr<int[]> col_index)
	: Matrix<T>(rows, cols, values_ptr), row_position(row_position), col_index(col_index)
	//CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index) : Matrix<T>(rows, cols, values_ptr), row_position(row_position), col_index(col_index)
{
	//--------------------------------------
	//  Constructor 2
	//  memory location passed by user
	//_______________________________________

	// need to overload the size of matrix
	this->size_of_values == nnzs;
}

template<class T>
CSRMatrix<T>::~CSRMatrix()
{
	//--------------------------------------
	//  Destructor
	// 
	//Since we are using shared pointer for Matrix
	// and CSRMatrix classes we dont need a specific 
	// destructor to account on whether we have ownership of
	// memory or not. Any object of these two classes will
	// be destroyed once the last copy goes out of scope
	//_______________________________________
}

template<class T>
void CSRMatrix<T>::printValues()
{
	//--------------------------------------
	//  Print Values
	// --------------------------------------
	std::cout << "Printing values:" << std::endl;
	for (int i = 0; i < nnzs; i++)
	{
		std::cout << this->values[i];
	}
	std::cout << std::endl;
}

template<class T>
void CSRMatrix<T>::printMatrix()
{	//--------------------------------------
	//  Print three arrays: values, col_index and 
	//  row_position
	// --------------------------------------
	std::cout << "Printing matrix" << std::endl;
	std::cout << "Values: ";
	for (int j = 0; j < this->nnzs; j++)
	{
		std::cout << this->values[j] << " ";
	}
	std::cout << std::endl;
	std::cout << "row_position: ";
	for (int j = 0; j < this->rows + 1; j++)
	{
		std::cout << this->row_position[j] << " ";
	}
	std::cout << std::endl;
	std::cout << "col_index: ";
	for (int j = 0; j < this->nnzs; j++)
	{
		std::cout << this->col_index[j] << " ";
	}
	std::cout << std::endl;
}

template<class T>
// The output of matrix by vector multiplication is a vector
// vector is small enough that we will consider it populated in dense form
void CSRMatrix<T>::matVecMult(Matrix<T>& vector, Matrix<T>& output)
// why did in class we use:  void CSRMatrix<T>::matVecMult(T* input, T* output)
{
	// Check if input and output are defined
	if (vector.values == nullptr || output.values == nullptr)
	{
		std::cout << "Vector or output are not created" << std::endl;
		return;
	}

	// initialise the output vector with zero
	for (int i = 0; i < this->rows; i++)
	{
		output.values[i] = 0;
		// why in class we use  output[i] = 0.0;?
	}

	// loop through each row
	for (int i = 0; i < this->rows; i++)
	{
		// for each element within that row (the indeces of first and last non-zero values in the row
		// are identified by row_position[i] and row_position[i+1]
		for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
		{
			output.values[i] += this->values[val_index] * vector.values[this->col_index[val_index]];
		}
	}

}

template<class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& output)
{
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

	// ?????????shouldnt we create the output in this function because we cant initiate the output object ahead of time

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
		output.values = std::shared_ptr<T[]>(new T[mat_left.rows * this->cols]); //new T[mat_left.rows * this->cols];
		output.preallocated = true;
	}

	for (int i = 0; i < output.nnzs; i++)
	{
		output.values[i] = 0;
	}
	std::cout << "start multi " << std::endl;

	int count_nnzs = 0;
	for (int i = 0; i < mat_left.rows; i++)
	{

		output.row_position[i] = 0 + count_nnzs;
		// loop through the non-zero elements in the row of mat_left
		for (int j = mat_left.row_position[i]; j < mat_left.row_position[i + 1]; j++)
		{
			// for each non zero column in A in this row find the corresponding row in "this" matrix
			int row_in_this = mat_left.col_index[j];
			for (int k = this->row_position[row_in_this]; k < this->row_position[row_in_this + 1]; k++)
			{
				output.values[count_nnzs] += mat_left.values[j] * this->values[k];
				output.col_index[count_nnzs] = mat_left.col_index[j];
				count_nnzs++;
			}

		}
	}
	output.row_position[mat_left.rows] = count_nnzs;
	std::cout << "there are nnzs: " << count_nnzs << std::endl;
}

template class CSRMatrix<double>;
template class CSRMatrix<int>;