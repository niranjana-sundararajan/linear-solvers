#pragma once
#include <iostream>
#include <fstream>
#include <string.h>
#include "Matrix.h"
#include "CSRMatrix.h"

template<class T>
// Read in a densely populated matrix from a comma delimited file
// MAKE SURE THE ROWS AND COLS of MATRIX CLASS MATCH THE NUMBER OF ELEMENTS IN A ROW AND NUMBER OF LINES
void readMatrix(std::string filename, Matrix<T>& output);


template<class T>
// Read in a densely populated matrix from a comma delimited file
void writeMatrix(std::string filename, Matrix<T>& output);


template<class T>
// Convert dense matrix into sparse format
void readSparse(Matrix<T>& dense, CSRMatrix<T>& sparse);