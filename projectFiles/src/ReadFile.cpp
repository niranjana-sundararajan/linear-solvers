#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include "CSRMatrix.h"
#include "ReadFile.h"

using namespace std;

template<class T>
// Read in a densely populated matrix from a comma delimited file
void readMatrix(std::string filename, Matrix<T>& output)
{
	// ------Description----------------------------------
	// Read a file in TXT format of a matrix, where each line 
	// corresponds to a row in a matrix. 
	// Populate the values of output matrix with the elements 
	// from the file.
	// *** Make sure the number of rows and columns in output
	// *** matrix matches to the number of lines and elements 
	// *** in the file 
	//-----------------------------------------------

	// Open the file 
	fstream matrixData;
	matrixData.open(filename);

	// Check if file read successfully
	if (matrixData.fail())
	{
		cerr << "Failed to open file " << filename << std::endl;
		//system("pause>0");
	}

	string line;
	int matrix_index = 0;

	// Loop through each line in the file
	while (matrixData.good())
	{
		// parse the line as string and convert string
		// to stringstream to allow conversion to double
		getline(matrixData, line);
		stringstream streamLine(line);
		string value;

		// split the line into elements with comma delimeter
		while (getline(streamLine, value, ','))
		{
			// fill in the given matrix location with values
			output.values[matrix_index] = stod(value);
			matrix_index++;
		}
	}
	// check number of elements in data matches to the output size
	if (matrix_index != (output.rows * output.cols))
	{
		cout << "Warning! The matrix size doesn't match the number of elements in file!" << endl;
		cout << "There are" << matrix_index << endl;
	}

	// Close the file
	matrixData.close();

	return;
}

template<class T>
// Read in a densely populated matrix from a comma delimited file
void writeMatrix(std::string filename, Matrix<T>& output)
{
	// ------Description----------------------------------
	// Write a matrix out to a CSV file, where each line 
	// corresponds to a row in a matrix. 
	//-----------------------------------------------

	fstream outputFile;
	outputFile.open(filename, ios::out);
	for (int i = 0; i < output.rows; i++)
	{
		for (int j = 0; j < output.cols; j++)
		{
			outputFile << output.values[i * output.cols + j] << ",";

		}
		outputFile << endl;

	}
	// close the file
	outputFile.close();

}

template <class T>
void readSparse(Matrix<T>& dense, CSRMatrix<T>& sparse)
{
	// ------Description----------------------------------
	// Convert a matrix in dense format into a CSR matrix
	// (row major sparse matrix format) 
	//-----------------------------------------------
	
	// Set counter to count number of non-zero values
	int count_nnzv = 0;
	// Set variable for dense matrix elements
	int size_of_values = dense.rows * dense.cols;
	int rows = dense.rows;
	int cols = dense.cols;

	// std::cout << " Rows " << rows << " cols " << cols << std :: endl;
	 // Calculate number of non-zero values in dense matrix
	for (int i = 0; i < size_of_values; i++)
	{

		cout << "i : " << i << endl;
		if (dense.values[i] != 0)
		{
			count_nnzv++;
		}

		// Set variables for sparse format
		// values , col_index and row_position in templated form
		// with appropriate sizes according to format

	   // set first row position as zero
		sparse.row_position[0] = 0;
		int nnzv = 0;

		// loops through the rows
		for (int i = 0; i < dense.rows; i++)
		{
			// loop through the columns
			for (int j = 0; j < dense.cols; j++)
			{
				// element index in dense matrix
				int element_index = i * cols + j;
				if (dense.values[element_index] != 0)
				{
					sparse.col_index[nnzv] = j;
					sparse.values[nnzv] = dense.values[element_index];

					nnzv++;
				}
			}
			// add the starting point of the next row, if at last row it will
			// add the total number of non zeros
			sparse.row_position[i + 1] = nnzv;
		}
	}
	sparse.nnzs = count_nnzv;
}

template class Matrix<double>;
template class Matrix<int>;