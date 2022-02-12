#pragma once
#include <iostream>
#include <iomanip> 
#include <string.h>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "Matrix.h"
#include "Solver.h"
#include "Interface.h"
#include "ReadFile.h"
#include "ReadFile.cpp"

using namespace std;



// INTERFACE (terminal prompt) class
// Constructor
Interface::Interface() 
{

}
Interface::~Interface() {}
// Function to ask user to launch an set up example or input their own data
void Interface::interfaceSelectInputOptions() {

	// clear screen
	system("CLS");

	cout << endl;
	cout << "*=============================================*" << endl;
	cout << "||                                           ||" << endl;
	cout << "||       	 LINEAR SOLVERS FOR          ||" << endl;
	cout << "||        	  DENSE AND SPARSE 	     ||" << endl;
	cout << "||       	    MATRIX SYSTEMS           ||" << endl;
	cout << "||                                           ||" << endl;
	cout << "*=============================================*" << endl;
	cout << endl;

	// user input 
	cout << "-------------------------------------------------" << endl;
	cout << "Please pick an option to input the data :" << std::endl;
	cout << " 1: Read from our Test Files." << endl;
	cout << " 2: Input your own Matrices." << endl;
	cout << " 3: Exit." << endl;
	cout << "-------------------------------------------------" << endl;
	cout << endl;

	char answer;
	cin >> answer;

	switch (answer)
	{
	case '1': {
		// launch function to determine dense or sparse structure
		interfaceDenseOrSparse();
		break;
	}
	case '2': {
		// interfaceReadFile();
		break;
	}
	case '3': {
		exit(0);
		break;
	}
	default: {
		invalidSelection(answer);
		interfaceSelectInputOptions();
	}
	}
}
// template<class T>
void Interface::invalidSelection(char c) {
	system("CLS");
	cout << endl << " Hi! Your selection : " << c << " is invalid." << endl;
	cout << endl;
	cout << " Press any key to restart." << endl;
	cin.get();
}

// Function to select the input is sparse matrix or dense matrix
// This will impact how the data is stored in memory
// ** use template so can return the object matrix (either Matrix or CSRMatrix -> this contains the size and location)
// template<class T>
void Interface::interfaceDenseOrSparse()
{

	// user input 
	cout << "-------------------------------------------------" << endl;
	std::cout << "Please pick one of the two types of Matrix solver" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << " 1: Dense Solvers." << endl;
	cout << " 2: Sparse Solvers." << endl;
	cout << "-------------------------------------------------" << endl;
	char answer;
	cin >> answer;

	// matrix set up
	switch (answer) {
	case '1': {
		interfaceSelectMatrixSizeDense();

	}
	case '2': {
		cout << "--------------------------------------------------------------------" << endl;
		cout << "Sparse matrix on interface is not available on the free version" << endl;
		cout << "Please upgrade to premium service to access this option on interface" << endl;
		cout << "OR" << endl;
		cout << "Please use as library through the examples found in examples folder." << endl;
		cout << "------------------------------------------------------------------" << endl;
		system("pause");
		interfaceDenseOrSparse();
	}
	default: {
		invalidSelection(answer);
		interfaceDenseOrSparse();
	}
	}
}
// template<class T>
void Interface::interfaceSelectMatrixSizeDense() {

	// Selection palatte for Matrix Sizes
	cout << "-------------------------------------------------" << endl;
	cout << "Please choose one of the following matrix sizes :" << endl;
	cout << "-------------------------------------------------" << endl;
	cout << " 1: 10x10(4x4)" << endl;
	cout << " 2: 50x50" << endl;
	cout << " 3: 100x100" << endl;
	cout << " 4: 200x200" << endl;
	cout << " 5: 400x400" << endl;
	cout << " 6: 800x800." << endl;
	cout << "-------------------------------------------------" << endl;
	char matrix_size;
	cin >> matrix_size;
	switch (matrix_size) {
		// matrix size = 10x10
	case '1':
	{
		this->size = 4;
		this->A = new Matrix<double>(size, size, true);
		readMatrix("../tests/mat4_A.txt", *(this->A));
		cout << "-------------------------------------------------" << endl;
		cout << " The matrix A you have chosen is as follows" << endl;
		this->A->printMatrix();
		cout << "-------------------------------------------------" << endl;
		this->B = new Matrix<double>(size, 1, true);
		readMatrix("../tests/mat4_b.txt", *(this->B));
		cout << " The matrix B you have chosen is as follows" << endl;
		this->B->printMatrix();
		cout << "-------------------------------------------------" << endl;
		system("pause");
		interfaceDenseSolver();
		break;
	}
	// matrix size = 50x50
	case '2': { // 
		this->size = 50;
		this->A = new Matrix<double>(size, size, true);
		cout << "-------------------------------------------------" << endl;
		readMatrix("../tests/output_100_100_A.txt", *(this->A));
		cout << " The matrix A has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;
		this->B = new Matrix<double>(size, 1, true);

		readMatrix("../tests/output_100_100_A.txt", *(this->B));
		cout << " The matrix B has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;

		system("pause");
		interfaceDenseSolver();
		break;
	}

	case '3': {
		this->size = 100;
		this->A = new Matrix<double>(size, size, true);
		cout << "-------------------------------------------------" << endl;
		readMatrix("../tests/output_100_100_A.txt", *(this->A));
		cout << " The matrix A has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;
		this->B = new Matrix<double>(size, 1, true);
		readMatrix("../tests/output_100_100_A.txt", *(this->B));
		cout << " The matrix B has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;

		system("pause");
		interfaceDenseSolver();
		break;
	}
			// matrix size = 200x200
	case '4': {
		this->size = 200;
		this->A = new Matrix<double>(size, size, true);
		cout << "-------------------------------------------------" << endl;
		readMatrix("../tests/output_200_200_A.txt", *(this->A));
		cout << " The matrix A has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;
		this->B = new Matrix<double>(size, 1, true);
		readMatrix("../tests/output_200_200_A.txt", *(this->B));
		cout << " The matrix B has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;

		system("pause");
		interfaceDenseSolver();
		break;
	}
			// matrix size = 400x400
	case '5': {
		this->size = 400;
		this->A = new Matrix<double>(size, size, true);
		cout << "-------------------------------------------------" << endl;
		readMatrix("../tests/output_400_400_A.txt", *(this->A));
		cout << " The matrix A has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;
		this->B = new Matrix<double>(size, 1, true);
		readMatrix("../tests/output_400_400_A.txt", *(this->B));
		cout << " The matrix B has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;

		system("pause");
		interfaceDenseSolver();
		break;
	}
			// matrix size = 800x800
	case '6': {
		this->size = 800;
		this->A = new Matrix<double>(size, size, true);
		cout << "-------------------------------------------------" << endl;
		readMatrix("../tests/output_800_800_A.txt", *(this->A));
		cout << " The matrix A has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;
		this->B = new Matrix<double>(size, 1, true);
		readMatrix("../tests/output_800_800_A.txt", *(this->B));
		cout << " The matrix B has been loaded" << endl;
		cout << "-------------------------------------------------" << endl;

		system("pause");
		interfaceDenseSolver();
		break;
	}

	default: {
		invalidSelection(matrix_size);
		// interfaceDenseOrSparse();
	}
	}
}


void Interface::interfaceDenseSolver() {

	// Print the available solvers for this matrix
	// Ask for user input
	char solverType;
	this->X = new Matrix<double>(size, 1, true);
	cout << "Choose of the following solvers" << endl;
	cout << " 1: Gauss Elimination" << endl;
	cout << " 2: LU Decomposition" << endl;
	cout << " 3: Cholesky" << endl;
	cout << " 4: Conjugate Gradient" << endl;
	cout << " 5: Return" << endl;

	cin >> solverType;
	switch (solverType) {
	case '1': {
		cout << ".............................................." << endl;
		cout << " ......Running Gauss Elimination............." << endl;
		cout << "............................................." << endl;
		// Set values of the X matrix to zero
		for (int i = 0; i < this->size; i++) {
			this->X->values[i] = 0;
		}
		// Run gaussian elimination on the matrices stored in class variables 
		// using the Solver class 
		this->denseSolvers->gaussElimination(*(this->A), *(this->B), *(this->X));
		cout << ".............................................." << endl;
		cout << " ..........THE LINEAR SYSTEM IS SOLVED......." << endl;
		cout << "............................................." << endl;
		cout << " Your output for X using Gaussian Elimination is :" << endl;
		this->X->printMatrix();
		cout << endl;
		break;
	}
	case '2': {
		cout << ".............................................." << endl;
		cout << " ......Running LU Decomposition............." << endl;
		cout << "............................................." << endl;
		// Set values of the X matrix to zero
		for (int i = 0; i < this->size; i++) {
			this->X->values[i] = 0;
		}
		// Run gaussian elimination on the matrices stored in class variables 
		// using the Solver class 
		this->denseSolvers->LUDecomposition_Solver(*(this->A), *(this->B), *(this->X));
		cout << ".............................................." << endl;
		cout << " ..........THE LINEAR SYSTEM IS SOLVED......." << endl;
		cout << "............................................." << endl;
		cout << " Your output for X using LU Decomposition is :" << endl;
		this->X->printMatrix();
		cout << endl;
		break;
	}
	case '3': {
		cout << ".............................................." << endl;
		cout << " ......Running Cholesky Decomposition............" << endl;
		cout << "............................................." << endl;
		// Set values of the X matrix to zero
		for (int i = 0; i < this->size; i++) {
			this->X->values[i] = 0;
		}
		// Run gaussian elimination on the matrices stored in class variables 
		// using the Solver class 
		this->denseSolvers->choleskySolver(*(this->A), *(this->B), *(this->X));
		cout << ".............................................." << endl;
		cout << " ..........THE LINEAR SYSTEM IS SOLVED......." << endl;
		cout << "............................................." << endl;
		cout << " Your output for X using Cholesky Decomposition is :" << endl;
		this->X->printMatrix();
		cout << endl;
		break;
	}
	case '4': {
		cout << ".............................................." << endl;
		cout << " ......Running Conjugat Gradient............" << endl;
		cout << "............................................." << endl;
		// Set values of the X matrix to zero
		for (int i = 0; i < this->size; i++) {
			this->X->values[i] = 0;
		}
		// Run gaussian elimination on the matrices stored in class variables 
		// using the Solver class 
		this->denseSolvers->conjugateGradient(*(this->A), *(this->B), *(this->X));
		cout << ".............................................." << endl;
		cout << " ..........THE LINEAR SYSTEM IS SOLVED......." << endl;
		cout << "............................................." << endl;
		cout << " Your output for X using Conjugate Gradient is :" << endl;
		this->X->printMatrix();
		cout << endl;
		break;
	}
	default: {
		invalidSelection(solverType);
		interfaceDenseSolver();
	}
	}
}
