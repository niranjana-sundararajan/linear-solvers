#pragma once
#include<iostream>
#include<memory>
#include "Matrix.h"

template<class T>
class CSRMatrix :public Matrix<T>
{
public:

	// define the class members
	int nnzs = nnzs;

	// make them integers because they are position values
	std::shared_ptr<int[]> row_position{ nullptr };
	std::shared_ptr<int[]> col_index{ nullptr };

	// Declare constructors
	CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
	//CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);
	CSRMatrix(int rows, int cols, int nnzs, std::shared_ptr<T[]> values_ptr, std::shared_ptr<int[]>  row_position, std::shared_ptr<int[]> col_index);
	// Declare destructor
	~CSRMatrix();

	// Declare member functions
	void printValues();
	void printMatrix();
	void matVecMult(Matrix<T>& vector, Matrix<T>& output);
	void matMatMult(CSRMatrix<T>& mat_left, CSRMatrix<T>& output);
};