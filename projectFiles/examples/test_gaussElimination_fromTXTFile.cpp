#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>

// #include "../include/Matrix.h"
#include "../include/Solver.h"
#include "../src/Solver.cpp"
#include "../include/ReadFile.h"
#include "../src/ReadFile.cpp"

using namespace std;


int main()
{
	int size = 4;

	// Read values

	auto* A = new Matrix<double>(size, size, true);
	auto* X = new Matrix<double>(size, 1, true);
	auto* B = new Matrix<double>(size, 1, true);

	// Read Matrix
	readMatrix("../tests/mat4_A.txt", *A );
	readMatrix("../tests/mat4_B.txt", *B );

	// Print the read matrix to terminal
	A->printMatrix();
	B->printMatrix();



    // set values of matrices
	for (int i = 0; i < size; i++){
		X->values[i] = 0;
	}

    // Use solver
	Solver<double> gauss;
	gauss.gaussElimination(*A, *B, *X);
   
   
	writeMatrix("gauss_output_4x4.txt", *X);

    delete[] A;
    delete[] X;
    delete[] B;

}
