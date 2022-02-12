#pragma once
#include "Matrix.h"
#include "Solver.h"

template <class T>
class Test : public Solver<T>
{
public:

	// These matrices can be access by user. They are populated with 
	//  random matrices when a certain constructor is called
	Matrix<T>* A = (nullptr);
	Matrix<T>* b = (nullptr);
	Matrix<T>* x = (nullptr);

	// Constructor for read in matrix (sparse/dense, size) --------
	Test(int matrixSize, bool dense);
	// Constructor to create a random matrix (symmetric/non-symmetric)
	Test();
	// Constructor for tolerance & iterations (from solver) -------
	Test(int iterations);

	// Constructor to generate random dense linear system
	Test(int size, bool dense, bool symmetric);


	//----- TESTS FOR SPEED AND PERFORMANCE -------------------------

	// Test of timing, absolute and relative errors and number of iterations for dence
	void TestDenseSolver();

	// Test of timing, absolute and relative errors and number of iterations for sparse
	void TestSparseSolver();

	// Test the solvers across different ill conditioned linear systems 
	void Test_IllCondition();

	// change some values in the vector b for ill conditioning
	void  changeB(Matrix<double>* b_name, Matrix<double>* originalB, int multiplier, int substract, int freq);

	// Test the solvers on one specific matrix to make sure there are no bugs
	void TestforAllSolvers();


	// Protected as the tests for timing and accuracy are using specific matrices
	// So we want to keep them protected from user and only allow for use within class
	// and allowing for potential subclasses in the future
protected:
	Matrix<double>* testAll_A = nullptr;
	Matrix<double>* testAll_b = nullptr;
	Matrix<double>* testAll_x = nullptr;

};