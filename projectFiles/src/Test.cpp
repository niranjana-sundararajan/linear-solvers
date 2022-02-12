#include <random>
#include <string>
#include <iostream>
#include <algorithm>
#include "Matrix.h"
#include "Test.h"
#include "ReadFile.h"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include "ReadFile.cpp"
#include <cstdlib>
#include <ctime>

using namespace std;

// ===========================================================================================
// -------- Coustructor of test class --------------------------------------------------------
// ===========================================================================================


template <class T>
// Constructor to generate random matrices--------------------------------
// user must confirm they want dense matrix
Test<T>::Test(int matrixSize, bool dense) : Solver<T>()
{

	// Create objects of Matrix class for each matrix in linear system
	this->A = new Matrix<T>(matrixSize, matrixSize, true);
	this->b = new Matrix<T>(matrixSize, 1, true);
	this->x = new Matrix<T>(matrixSize, 1, true);

	// Read the datafile

	string fileName_A = "mat4_A.txt";
	string fileName_b = "mat4_b.txt";
	string filename_x = "mat4_x.txt";

	// read the matrix and save to the location
	readMatrix(fileName_A, *A);
	readMatrix(fileName_b, *b);
	readMatrix(filename_x, *x);
}

template <class T>
Test<T>::Test() : Solver<T>()
{
	// Empty constructor to allow for access to both solver functions and test functions
}

template <class T>
Test<T>::Test(int iterations) : Solver<T>(iterations)
{
	// Constructor to allow for changes in number of max iterations
}


template <class T>
Test<T>::Test(int size, bool dense, bool symmetric)
{
	// --------------------------------------------------------------------------
	// Constructor to generate random matrices of different sizes
	// 
	// --------------------------------------------------------------------------

	if (symmetric)
	{
		this->A = new Matrix<T>(size, size, true);
		this->b = new Matrix<T>(size, 1, true);
		this->x = new Matrix<T>(size, 1, true);

		// Set a seed number for the random number generator
		srand(time(0));

		// Create a random dense symmetric matrix A
		// Allocate memory for transpose
		Matrix<T> A_transpose = Matrix<T>(size, size, true);
		Matrix<T> A_temp = Matrix<T>(size, size, true);
		for (int i = 0; i < size * size; i++)
		{
			std::random_device random_generator;
			std::mt19937 random_engine(random_generator());
			std::uniform_int_distribution<> distr(1, 12984);

			double const randomNumber = distr(random_engine);

			// generate random integer from 1 to 100 for each matrix entry
			A_temp.values[i] = randomNumber / 9984.;
			A_transpose.values[i] = A_temp.values[i];
		}
		// Transpose the matrix and add to self and half to get the symmetric matrix
		A_transpose.transpose();
		A_temp.matMatMult(A_transpose, *this->A);

		vector<double> vs = { 1, 2, 3, 4, 5, 60000, 70000, 80000, 90000, 100000 };
		//vector<double> vs = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		// Generate the random vector x
		for (int i = 0; i < size; i++)
		{
			// generate random integer from 1 to 100 for each matrix entry
			//this->x->values[i] = rand() % 40 + 1;
			this->x->values[i] = vs[i];
		}
		// Caluclate vector b in A*x = b
		(this->A)->matVecMult(*this->x, *this->b);

	}

}
// ===========================================================================================
// -------- Test for All Solvers Function ----------------------------------------------------
// ===========================================================================================
template <class T>
void Test<T>::TestforAllSolvers()
{
	// --------------------------------------------------------------------------
	// Function to test all solvers on a 10x10 linear system with ---------------
	// known exact solutioin. This linear system was created by using randon 
	// symmetric positive definite matrix generator and x vector. that are multiplied 
	// together to find the vector b. This vector b was then used within the solvers
	// to backengineer the vector x.
	// The required tolerance set here is 0.00001--------------------------------
	// --------------------------------------------------------------------------

	// Read the datafile

	string fileName_A = "output_A.txt";
	string fileName_b = "output_b.txt";
	string filename_x = "output_x.txt";

	this->testAll_A = new Matrix<double>(10, 10, true);
	this->testAll_b = new Matrix<double>(10, 1, true);
	this->testAll_x = new Matrix<double>(10, 1, true);

	// read the matrix and save to the location
	readMatrix(fileName_A, *(this->testAll_A));
	readMatrix(fileName_b, *(this->testAll_b));
	readMatrix(filename_x, *(this->testAll_x));

	// Set up clean system for solver
	Matrix<double> A_copy = Matrix<double>(testAll_A->rows, testAll_A->rows, true);
	Matrix<double> b_copy = Matrix<double>(testAll_b->rows, 1, true);
	Matrix<double> output = Matrix<double>(testAll_x->rows, 1, true);

	// loop through all the solvers testing on the same linear system
	for (int i = 1; i < 7; i++)
	{

		// Copy over the values and overwrite output with 0 as required by iterative solvers
		for (int i = 0; i < testAll_A->rows * testAll_A->cols; i++)
		{
			if (i < testAll_A->cols)
			{
				A_copy.values[i] = testAll_A->values[i];
				b_copy.values[i] = testAll_b->values[i];
				output.values[i] = 0;
			}
			else
			{
				A_copy.values[i] = testAll_A->values[i];
			}
		}
		// based on the iteration run the solvers in a sequence
		switch (i)
		{
		case 1:
			cerr << "// -------------------------------------------------------------" << endl;
			cerr << "// ----Test Gauss Elimination-----------------------------------" << endl;
			cerr << "// ---------------------------------------------------------------" << endl;

			//this->gaussElimination(A_copy, b_copy, output);
			break;

		case 2:
			cerr << "// -------------------------------------------------------------" << endl;
			cerr << "// ----Test LU Decomposition------------------------------------" << endl;
			cerr << "// ---------------------------------------------------------------" << endl;

			//this->LUDecomposition_Solver(A_copy, b_copy, output);
			break;

		case 3:
			cerr << "// -------------------------------------------------------------" << endl;
			cerr << "// ----Test gaussSeidel-----------------------------------------" << endl;
			cerr << "// ---------------------------------------------------------------" << endl;

			//this->gaussSeidel(A_copy, b_copy, output);
			break;

		case 4:
			cerr << "// -------------------------------------------------------------" << endl;
			cerr << "// ----Test choleskySolver--------------------------------------" << endl;
			cerr << "// ---------------------------------------------------------------" << endl;

			//this->choleskySolver(A_copy, b_copy, output);
			break;

		case 5:
			cerr << "// -------------------------------------------------------------" << endl;
			cerr << "// ----Test Jacobi----------------------------------------------" << endl;
			cerr << "// ---------------------------------------------------------------" << endl;

			//this->jacobiDense(A_copy, b_copy, output);
			break;

		case 6:
			cerr << "// -------------------------------------------------------------" << endl;
			cerr << "// ----Test Conjugate Gradient-----------------------------------" << endl;
			cerr << "// ---------------------------------------------------------------" << endl;

			//this->conjugateGradient(A_copy, b_copy, output);
			break;

		}

		// compare the results from solvers against the exact solution. 
		output -= *this->testAll_x;

		// if the error magnitude is below tolerance
		if (output.L2norm() < 0.00001)
		{
			cerr << "Yays! It works!!!" << endl;
			cerr << "The number of iterations is: " << this->numiterationsDone << endl;
			cerr << "error norm is:" << output.L2norm() << endl << endl;
			cerr << "// ---------------------------------------------------------------" << endl;
		}
		else
		{
			cerr << "Oh no! It failed has failed :((((((((" << endl;
			cerr << "the number of iterations is: " << this->numiterationsDone << endl;
			cerr << "error norm is:" << output.L2norm() << endl << endl;
			cerr << "// ---------------------------------------------------------------" << endl;
		}

	}
	return;
}


template <class T>
void Test<T>::changeB(Matrix<double>* b_name, Matrix<double>* originalB, int multiplier, int substract, int freq)
{
	//------------------------------------------------
	// Function to change certain values in a vector
	// by multiplying it and subtracting from each value
	// ------------------------------------------------

	for (int i = 0; i < originalB->rows; i += freq)
	{
		b_name->values[i] = originalB->values[i] * multiplier - substract;
		cout << b_name->values[i] << "   " << originalB->values[i] << endl;
	}
}

// ===========================================================================================
// -------- Test for Dense Solver Function ---------------------------------------------------
// ===========================================================================================

template <class T>
void Test<T>::TestDenseSolver()
{
	//------------------------------------------------
	// Test function that reads in different size matrices
	// and for each matrix it runs each solvers 5 times to
	// measure the time taken to calculate the solution
	// (we calculate 5 times so we can report on the average 
	// time for the analysis)
	// ------------------------------------------------
	// 
	// Read the datafile and run tests---------------------------------

	// Read the datafile and run tests ----------------------------------------------
	for (int sizeMat = 2; sizeMat < 3; sizeMat++)
	{
		int size = 0;
		string fileName_A;
		string fileName_b;
		string filename_x;
		// Read text files of size 10 random dense matrices
		if (sizeMat == 0)
		{
			// File name
			 fileName_A = "output_10_10_A.txt";
			 fileName_b = "output_10_10_b.txt";
			 filename_x = "output_10_10_x.txt";

			// Create a Matrix A (10x10), Matrix b (10x1) and initial guess of x (10x1)
			size = 50;

			// Tell a user when it finishes reading size 10 matrices
			cout << "Finished reading size 10 matrices" << endl;
		}
		// Read text files of size 50 random dense matrices
		else if (sizeMat == 1)
		{
			// File name
			 fileName_A = "output_50_50_A.txt";
			 fileName_b = "output_50_50_b.txt";
			 filename_x = "output_50_50_x.txt";

			size = 50;

			// Tell a user when it finishes reading size 50 matrices
			cout << "Finished reading size 50 matrices" << endl;
		}
		// Read text files of size 100 random dense matrices
		else if (sizeMat == 2)
		{
			// File name
			 fileName_A = "output_100_100_A.txt";
			 fileName_b = "output_100_100_b.txt";
			 filename_x = "output_100_100_x.txt";

			size = 100;

			// Tell a user when it finishes reading size 100 matrices
			cout << "Finished reading size 100 matrices" << endl;
		}
		// Read text files of size 200 random dense matrices
		else if (sizeMat == 3)
		{
			// File name
			 fileName_A = "output_200_200_A.txt";
			 fileName_b = "output_200_200_b.txt";
			 filename_x = "output_200_200_x.txt";

			size = 200;

			// Tell a user when it finishes reading size 200 matrices
			cout << "Finished reading size 200 matrices" << endl;
		}
		// Read text files of size 400 random dense matrices
		else if (sizeMat == 4)
		{
			// File name
			 fileName_A = "output_400_400_A.txt";
			 fileName_b = "output_400_400_b.txt";
			 filename_x = "output_400_400_x.txt";

			size = 400;

			// Tell a user when it finishes reading size 400 matrices
			cout << "Finished reading size 400 matrices" << endl;
		}
		// Read text files of size 600 random dense matrices
		else if (sizeMat == 5)
		{
			// File name
			 fileName_A = "output_600_600_A.txt";
			 fileName_b = "output_600_600_b.txt";
			 filename_x = "output_600_600_x.txt";

			size = 600;

			// Tell a user when it finishes reading size 600 matrices
			cout << "Finished reading size 600 matrices" << endl;
		}
		// Read text files of size 800 random dense matrices
		else if (sizeMat == 6)
		{
			// File name
			 fileName_A = "output_800_800_A.txt";
			 fileName_b = "output_800_800_b.txt";
			 filename_x = "output_800_800_x.txt";

			size = 800;

			// Tell a user when it finishes reading size 800 matrices
			cout << "Finished reading size 800 matrices" << endl;
		}

		// Read text files of size 150 random dense matrices


		// Create a Matrix A (800x800), Matrix b (800x1) and initial guess of x (800x1)
		this->testAll_A = new Matrix<double>(size, size, true);
		this->testAll_b = new Matrix<double>(size, 1, true);
		this->testAll_x = new Matrix<double>(size, 1, true);

		// read the matrix and save to the location
		readMatrix(fileName_A, *(this->testAll_A));
		readMatrix(fileName_b, *(this->testAll_b));
		readMatrix(filename_x, *(this->testAll_x));

		cout << endl;

		// calculate k times for same size of matrices
		for (int k = 0; k < 1; k++)
		{
			// --------------------------------------------------------------------
			// ----Test Gauss Elimination------------------------------------------
			// --------------------------------------------------------------------

			// Create copies that may be overwritten within a solver
			Matrix<double>* A_copyGE = new Matrix<double>(testAll_A->rows, testAll_A->rows, true);
			Matrix<double>* b_copyGE = new Matrix<double>(testAll_b->rows, 1, true);
			Matrix<double>* outputGE = new Matrix<double>(testAll_x->rows, 1, true);

			for (int i = 0; i < testAll_A->rows * testAll_A->cols; i++)
			{
				if (i < testAll_A->cols)
				{
					A_copyGE->values[i] = testAll_A->values[i];
					b_copyGE->values[i] = testAll_b->values[i];
				}
				else
				{
					A_copyGE->values[i] = testAll_A->values[i];
				}
			}

			// Time (sec) ------------------
			clock_t start1 = clock();
			//this->gaussElimination(*A_copyGE, *b_copyGE, *outputGE);
			clock_t end1 = clock();
			cout << "Gauss Elimination\t" << k << "  " << (double)(end1 - start1) / (double)(CLOCKS_PER_SEC) << endl;

			// Output ----------------------
			writeMatrix("output_400_400_output_GE.txt", *outputGE);
			outputGE->printMatrix();

			// Absolute error --------------
			//cout << "GE absolute error : " << this->errorAbsoluteCalc(*this->testAll_A, *this->testAll_b, *outputGE) << endl;

			// Rerative error --------------
			//cout << "GE relative error : " << this->errorRelativeCalc(*outputGE, *this->testAll_x) << endl;


			// Delete memories -------------
			delete[] A_copyGE;
			delete[] b_copyGE;
			delete[] outputGE;

			cout << endl;

			//// --------------------------------------------------------------------
			//// ----Test LU Decomposition solver -----------------------------------
			//// --------------------------------------------------------------------

			// Create copies that may be overwritten within a solver
			Matrix<double>* A_copyLU = new Matrix<double>(this->testAll_A->rows, this->testAll_A->rows, true);
			Matrix<double>* b_copyLU = new Matrix<double>(this->testAll_b->rows, 1, true);
			Matrix<double>* outputLU = new Matrix<double>(this->testAll_x->rows, 1, true);

			for (int i = 0; i < this->testAll_A->rows * this->testAll_A->cols; i++)
			{
				if (i < this->testAll_A->cols)
				{
					A_copyLU->values[i] = this->testAll_A->values[i];
					b_copyLU->values[i] = this->testAll_b->values[i];
				}
				else
				{
					A_copyLU->values[i] = this->testAll_A->values[i];
				}
			}

			// Time (sec) ------------------
			clock_t start2 = clock();
			//this->LUDecomposition_Solver(*A_copyLU, *b_copyLU, *outputLU);
			clock_t end2 = clock();
			cout << "LU Decomposition\t" << k << "  " << (double)(end2 - start2) / (double)(CLOCKS_PER_SEC) << endl;

			// Output ----------------------
			writeMatrix("output_400_400_output_LU.txt", *outputLU);
			outputLU->printMatrix();

			// Absolute error --------------
			//cout << "LU absolute error : " << this->errorAbsoluteCalc(*this->testAll_A, *this->testAll_b, *outputLU) << endl;

			// Lirative error --------------
			//cout << "LU relative error : " << this->errorRelativeCalc(*outputLU, *this->testAll_x) << endl;

			// Delete memories -------------
			delete[] A_copyLU;
			delete[] b_copyLU;
			delete[] outputLU;

			cout << endl;

			// --------------------------------------------------------------------
			// ----Test gaussSeidel------------------------------------------------
			// --------------------------------------------------------------------

			// Create copies that may be overwritten within a solver
			Matrix<double>* A_copyGS = new Matrix<double>(testAll_A->rows, testAll_A->rows, true);
			Matrix<double>* b_copyGS = new Matrix<double>(testAll_b->rows, 1, true);
			Matrix<double>* outputGS = new Matrix<double>(testAll_x->rows, 1, true);

			for (int i = 0; i < testAll_A->rows * testAll_A->cols; i++)
			{
				if (i < testAll_A->cols)
				{
					A_copyGS->values[i] = testAll_A->values[i];
					b_copyGS->values[i] = testAll_b->values[i];
				}
				else
				{
					A_copyGS->values[i] = testAll_A->values[i];
				}
			}

			// Time (sec) ----------------
			clock_t start3 = clock();
			//this->gaussSeidel(*A_copyGS, *b_copyGS, *outputGS);
			clock_t end3 = clock();
			cout << "Gauss Seidel\t" << k << "  " << (double)(end3 - start3) / (double)(CLOCKS_PER_SEC) << endl;

			// Output --------------------
			//writeMatrix("output_400_400_output_GS.txt", *outputGS);
			//*outputGS->printMatrix();

			// Absolute error ------------
			//cout << "GS absolute error : " << this->errorAbsoluteCalc(*this->testAll_A, *this->testAll_b, *outputGS) << endl;

			// Relative error ------------
			//cout << "GS relative error : " << this->errorRelativeCalc(*outputGS, *this->testAll_x) << endl;

			// Number of iterations ------
			cout << "Number of iteration done = " << this->numiterationsDone << endl;

			// Delete memories -----------
			delete[] A_copyGS;
			delete[] b_copyGS;
			delete[] outputGS;

			cout << endl;

			// --------------------------------------------------------------------
			// ----Test choleskySolver --------------------------------------------
			// --------------------------------------------------------------------

			// Create copies that may be overwritten within a solver
			Matrix<double>* A_copyCH = new Matrix<double>(testAll_A->rows, testAll_A->rows, true);
			Matrix<double>* b_copyCH = new Matrix<double>(testAll_b->rows, 1, true);
			Matrix<double>* outputCH = new Matrix<double>(testAll_x->rows, 1, true);

			for (int i = 0; i < testAll_A->rows * testAll_A->cols; i++)
			{
				if (i < testAll_A->cols)
				{
					A_copyCH->values[i] = testAll_A->values[i];
					b_copyCH->values[i] = testAll_b->values[i];
				}
				else
				{
					A_copyCH->values[i] = testAll_A->values[i];
				}
			}

			// Time ----------------------
			clock_t start4 = clock();
			//this->choleskySolver(*A_copyCH, *b_copyCH, *outputCH);
			clock_t end4 = clock();
			cout << "Cholesky Solver\t" << k << "  " << (double)(end4 - start4) / (double)(CLOCKS_PER_SEC) << endl;

			// Output --------------------
			writeMatrix("output_400_400_output_CH.txt", *outputCH);

			// Absolute error ------------
			//cout << "CH absolute error : " << this->errorAbsoluteCalc(*A_copyCH, *b_copyCH, *outputCH) << endl;

			// Relative error ------------
			//cout << "CH relative error : " << this->errorRelativeCalc(*outputCH, *this->testAll_x) << endl;

			// Delete memories -----------
			delete[] A_copyCH;
			delete[] b_copyCH;
			delete[] outputCH;

			cout << endl;

			// --------------------------------------------------------------------
			// ----Test Conjugate Gradient-----------------------------------------
			// --------------------------------------------------------------------

			// Create copies that may be overwritten within a solver
			Matrix<double>* A_copyCG = new Matrix<double>(testAll_A->rows, testAll_A->rows, true);
			Matrix<double>* b_copyCG = new Matrix<double>(testAll_b->rows, 1, true);
			Matrix<double>* outputCG = new Matrix<double>(testAll_x->rows, 1, true);

			for (int i = 0; i < testAll_A->rows * testAll_A->cols; i++)
			{
				if (i < testAll_A->cols)
				{
					A_copyCG->values[i] = testAll_A->values[i];
					b_copyCG->values[i] = testAll_b->values[i];
					testAll_x->values[i] = testAll_x->values[i];
				}
				else
				{
					A_copyCG->values[i] = testAll_A->values[i];
				}
			}

			// Time (sec) ---------------
			clock_t start5 = clock();
			//this->conjugateGradient(*A_copyCG, *b_copyCG, *outputCG);
			clock_t end5 = clock();
			cout << "conjugateGradient\t" << k << "  " << (double)(end5 - start5) / (double)(CLOCKS_PER_SEC) << endl;

			// Output -------------------
			writeMatrix("output_400_400_output_CG.txt", *outputCG);
			outputCG->printMatrix();

			// Absolute error ------------
			//cout << "CG absolute error : " << this->errorAbsoluteCalc(*A_copyCG, *b_copyCG, *outputCG) << endl;

			// Relative error ------------
			//cout << "CG relative error : " << this->errorRelativeCalc(*outputCG, *this->testAll_x) << endl;

			// Number of iterations ------
			cout << "Number of iteration done = " << this->numiterationsDone << endl;

			// Delete memories -----------
			delete[] A_copyCG;
			delete[] b_copyCG;
			delete[] outputCG;

			cout << endl;

			// ------------------------------------------------------------------------
			// ----Test jacobiDense ---------------------------------------------------
			// ------------------------------------------------------------------------

			// Create copies that may be overwritten within a solver
			Matrix<double>* A_copyJA = new Matrix<double>(testAll_A->rows, testAll_A->rows, true);
			Matrix<double>* b_copyJA = new Matrix<double>(testAll_b->rows, 1, true);
			Matrix<double>* outputJA = new Matrix<double>(testAll_x->rows, 1, true);

			for (int i = 0; i < testAll_A->rows * testAll_A->cols; i++)
			{
				if (i < testAll_A->cols)
				{
					A_copyJA->values[i] = testAll_A->values[i];
					b_copyJA->values[i] = testAll_b->values[i];
				}
				else
				{
					A_copyJA->values[i] = testAll_A->values[i];
				}
			}

			// Time (sec) ----------------
			clock_t start6 = clock();
			//this->jacobiDense(*A_copyJA, *b_copyJA, *outputJA);
			clock_t end6 = clock();
			cout << "jacobiDense\t" << (double)(end6 - start6) / (double)(CLOCKS_PER_SEC) << endl;

			// Output --------------------
			writeMatrix("output_800_800_output_JA.txt", *outputJA);

			// Absolute error ------------
			//cout << "JA absolute error : " << this->errorAbsoluteCalc(*A_copyJA, *b_copyJA, *outputJA);

			// Relative error ------------
			//cout << "JA relative error : " << this->errorRelativeCalc(*outputJA, *this->testAll_x) << endl;

			// Number of iterations ------
			cout << "Number of iteration done = " << this->numiterationsDone << endl;

			// Delete memories -----------
			delete[] A_copyJA;
			delete[] b_copyJA;
			delete[] outputJA;
		}
		cout << endl;
	}
}

// ===========================================================================================
// -------- Ill Conditioned Problems----------------------------------------------------------
// ===========================================================================================

template <class T>
void Test<T>::Test_IllCondition()
{

	for (int i = 0; i < 6; i++)
	{
		//------------------------------------------------------------------------------------------------------
		// test how solvers work on ill condditioned linear systems
		// We explored: Wilson, Hilbert, Wilkinson and Rutishauser. 
		// All matrices have different increasing condition number
		//------------------------------------------------------------------------------------------------------

		// Uncomment the matrices to run the one of interest
		// ------------------------Wilson---------------------------
		//std::shared_ptr<double[]> A_Array1(new double[16]{ 5.,7.,6.,5.,7.,10.,8.,7.,6.,8.,10.,9.,5.,7.,9.,10. }, [](double* a) {delete[] a;});
		//Matrix<double>* A = new Matrix<double>(4, 4, A_Array1);
		//cerr << "----------------Wilson------------------" << endl;

		//// --------------------------------------------------------
		////-----------------------Hilbert --------------------------

		//*double Arrays[16] = {1.,1. / 2.,1. / 3,1. / 4,1. / 2,1. / 3,1. / 4,1. / 5,1. / 3.,1. / 4.,1. / 5.,1. / 6,1. / 4.,1. / 5,1. / 6,1. / 7};
		//Matrix<double>* A = new Matrix<double>(4, 4, true);
		//cerr << "----------------Hilbert------------------" << endl;

		////-----------------------Wilkinsons matrix --------------------------
		std::shared_ptr<double[]> A_Array1(new double[16]{ 0.00009143, 0.,0.,0.,0.8762,0.00007156,0.,0.,0.7943,0.8143,0.00009504,0,0.8017,0.6123,0.7165,0.00007123 }, [](double* a) {delete[] a;});
		Matrix<double> A = Matrix<double>(4, 4, A_Array1);
		//// --------------------------------------------------------
		// 
		//std::shared_ptr<double[]> A_Array1(new double[16]{ 10.,1.,4.,0.,1.,10.,5.,-1.,4.,5.,10.,7.,0.,-1.,7.,9., }, [](double* a) {delete[] a;});
		//Matrix<double>* A = new Matrix<double>(4, 4, A_Array1);
		//cerr << "----------------Rutishauser ------------------" << endl;
		//// --------------------------------------------------------

		//matrix x for the solution
		std::shared_ptr<double[]> x_array(new double[4]{ 0.,0.,0.,0. }, [](double* a) {delete[] a;});
		Matrix<double> x = Matrix<double>(4, 1, x_array);

		// Allow for different b vectors
		//std::shared_ptr<double[]> b_Array(new double[16]{ 1.,1.,1.,1. }, [](double* a) {delete[] a;});
		std::shared_ptr<double[]> b_Array(new double[16]{ 0.00009143,0.87627156, 1.60869504,2.13057123 }, [](double* a) {delete[] a;});
		Matrix<double> b = Matrix<double>(4, 1, b_Array);

		A.printMatrix();
		b.printMatrix();

		switch (i)
		{
		case 0:
			cerr << "---------------------------" << endl;
			cout << "Cholesky" << endl;
			//this->choleskySolver(A, b, x);
			x.printMatrix();
			cerr << "---------------------------" << endl;
			break;
		case 1:
			cerr << "---------------------------" << endl;
			cout << "LU Decomposition" << endl;
			//this->LUDecomposition_Solver(A, b, x);
			//x.printMatrix();
			cerr << "---------------------------" << endl;
			break;
		case 2:
			cerr << "---------------------------" << endl;
			cout << "Conjugate Gradient" << endl;
			//this->conjugateGradient(A, b, x, 0.000001);
			x.printMatrix();
			cerr << "---------------------------" << endl;
			break;
		case 3:
			cerr << "---------------------------" << endl;
			cout << "Gauss Elimination" << endl;
			//this->gaussElimination(A, b, x);
			x.printMatrix();
			cerr << "---------------------------" << endl;
			break;
		case 4:
			cerr << "---------------------------" << endl;
			cout << "Gauss Siedel" << endl;
			//this->gaussSeidel(A, b, x, false, 0.000001);
			x.printMatrix();
			cerr << "---------------------------" << endl;
			break;
		case 5:
			cerr << "---------------------------" << endl;
			cout << "LU without PP" << endl;
			auto L = Matrix<double>(A.rows, A.cols, true);
			A.LUDecomposition(L);
			// Using forward substitution find y in the linear system
			Matrix<double> y = Matrix<double>(A.rows, 1, true);
			//this->forwardSubstitution(L, y, b);

			//this->backSubstitution(A, x, y);
			x.printMatrix();
			cerr << "---------------------------" << endl;
			break;

		}
	}


}

// ===========================================================================================
// -------- Test for Sparse Solver Function --------------------------------------------------
// ===========================================================================================

template <class T>
void Test<T>::TestSparseSolver()
{

	// Read the datafile and run tests----------------------------------------------
	for (int sizeMat = 0; sizeMat < 1; sizeMat++)
	{
		string fileName_A;
		string fileName_b;
		string filename_x;

		int size = 0;
		// Read a text file of size 1000 random sparse matrices
		if (sizeMat == 0)
		{
			// File name
			string fileName_A = "randomCRS_1000A.txt";
			string fileName_b = "randomCRS_1000b.txt";
			string filename_x = "randomCRS_1000x.txt";

			size = 1000;
		}
		// Read a text file of size 5000 random sparse matrices
		else if (sizeMat == 1)
		{
			// File name
			string fileName_A = "randomCRS_5000A.txt";
			string fileName_b = "randomCRS_5000b.txt";
			string filename_x = "randomCRS_5000x.txt";

			size = 5000;

		}
		// Read a text file of size 100 random matrices
		else if (sizeMat == 2)
		{
			// File name
			string fileName_A = "randomCRS_100A.txt";
			string fileName_b = "randomCRS_100b.txt";
			string filename_x = "randomCRS_100x.txt";

			size = 100;

		}
		// Read a text file of size 200 random matrices
		else if (sizeMat == 3)
		{
			// File name
			string fileName_A = "randomCRS_200A.txt";
			string fileName_b = "randomCRS_200b.txt";
			string filename_x = "randomCRS_200x.txt";

			size = 200;

		}
		// Read a text file of size 400 random matrices
		else if (sizeMat == 4)
		{
			// File name
			string fileName_A = "randomCRS_400A.txt";
			string fileName_b = "randomCRS_400b.txt";
			string filename_x = "randomCRS_400x.txt";

			size = 400;
		}
		// Read a text file of size 600 random matrices
		else if (sizeMat == 5)
		{
			// File name
			string fileName_A = "randomCRS_600A.txt";
			string fileName_b = "randomCRS_600b.txt";
			string filename_x = "randomCRS_600x.txt";

			size = 600;
			
		}
		// Read a text file of size 800 random matrices
		else if (sizeMat == 6)
		{
			// File name
			string fileName_A = "randomCRS_800A.txt";
			string fileName_b = "randomCRS_800b.txt";
			string filename_x = "randomCRS_800x.txt";

			size = 800;

		}

		// Create a Matrix A (10x10), Matrix b (10x1) and initial guess of x (10x1)
		this->testAll_A = new Matrix<double>(size, size, true);
		this->testAll_b = new Matrix<double>(size, 1, true);
		this->testAll_x = new Matrix<double>(size, 1, true);
		CSRMatrix<double> sparse = CSRMatrix<double>(size, size, this->testAll_A->sparsity(), true);

		// read the matrices and save to the location
		readMatrix(fileName_A, *(this->testAll_A));
		readMatrix(fileName_b, *(this->testAll_b));
		readMatrix(filename_x, *(this->testAll_x));

		readSparse(*testAll_A, sparse);
		cout << endl;

		// number of times to run
		for (int k = 0; k < 1; k++)
		{
			// ------------------------------------------------------------------------
			// ----Test jacobi Sparse -------------------------------------------------
			// ------------------------------------------------------------------------

			//Copy the A matrix again
			Matrix<double> outputJS = Matrix<double>(testAll_x->rows, 1, true);

			// Time (sec) ----------------
			clock_t start6 = clock();
			//this->jacobiSparse(sparse, *(this->testAll_b), outputJS);
			clock_t end6 = clock();
			cout << "jacobiSparse\t" << k << "  " << (double)(end6 - start6) / (double)(CLOCKS_PER_SEC) << endl;

			// Absolute error ------------
			//cout << "JS absolute error : " << this->errorAbsoluteCalc(sparse, *(this->testAll_b), *outputJS);

			// Number of iterations ------
			cout << "Number of iteration done = " << this->numiterationsDone << endl;

			cout << endl;

			// --------------------------------------------------------------------
			// ----Test Conjugate Gradient Sparse ---------------------------------
			// --------------------------------------------------------------------
			Matrix<double> outputCGS = Matrix<double>(testAll_x->rows, 1, true);

			// Time (sec) ---------------
			clock_t start5 = clock();
			//this->conjugateGradientSparse(sparse, *this->testAll_b, outputCGS);
			clock_t end5 = clock();
			cout << "conjugateGradient\t" << k << "  " << (double)(end5 - start5) / (double)(CLOCKS_PER_SEC) << endl;

			//// Output -------------------
			writeMatrix("randomCRS_1000output_CGS.txt", outputCGS);


			// Absolute error ------------
			//cout << "CGS absolute error : " << this->errorAbsoluteCalc(*this->testAll_A, *this->testAll_b, outputCGS) << endl;
			// Number of iterations ------
			cout << "Number of iteration done = " << this->numiterationsDone << endl;

			cout << endl;
		}

	}
}
template Test<double>;
template Test<int>;