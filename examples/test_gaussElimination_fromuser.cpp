#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>

// #include "../include/Matrix.h"
#include "../include/Solver.h"
#include "../src/Solver.cpp"


using namespace std;


int main()
{
    // // initialize values
	// int size = 3;
	// double a[] = { 2, 1, 1 , -3, -1, 2, -2 , 1, 2};
	// double x[] = { 0, 0, 0 };
	// double b[] = { 8, -11, -3 };

	// Read values
	auto* A = new Matrix<double>(size, size, true);
	auto* X = new Matrix<double>(size, 1, true);
	auto* B = new Matrix<double>(size, 1, true);


	readMatrix("../tests/mat4_A.txt", *A );
	A->printMatrix();

    // create matrix pointers
	auto* A = new Matrix<double>(size, size, true);
	auto* X = new Matrix<double>(size, 1, true);
	auto* B = new Matrix<double>(size, 1, true);

    // set values of matrices
	for (int i = 0; i < size * size; i++){
		A->values[i] = a[i];
	}

	for (int i = 0; i < size; i++){
		B->values[i] = b[i];
		X->values[i] = x[i];
	}

    // Use solver

	Solver<double> gauss;
	gauss.gaussElimination(*A, *B, *X);
    X->printMatrix(); 

    delete[] A;
    delete[] X;
    delete[] B;

}
