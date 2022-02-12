#pragma once
#include <cstdio>
#include <iostream>
#include "Matrix.h"
#include "ReadFile.h"
#include "Solver.h"

// Add other header files as required
using namespace std;

// template <class T>
class Interface
{
private:
	// declare variables here!!



public:

	Matrix<double>* A = nullptr;
	Matrix<double>* B = nullptr;
	Matrix<double>* X = nullptr;
	Solver<double>* denseSolvers;
	Solver<double>* sparseSolvers;
	int size = -1;
	// -------------- CONSTRUCTOR AND DESTRUCTOR ----------------------
	Interface();
	~Interface();
	//-----------------------------------------------------------------

	// ---------------INTERFACE FUNCTIONS --------------------------------
	void interfaceSelectMatrixSizeDense(); // Function to select matrix size

	void interfaceSelectMatrixSizeSparse(); // Function to select matrix size

	void interfaceSelectInputOptions(); // Function to either upload own file or select demo matrix from matrix_data 

	void interfaceDenseOrSparse();

	void invalidSelection(char message); // Invaild input selection/file message

	void interfaceDenseSolver(); // Function to specify which dense solver to use
	//--------------------------------------------------------------------

	// ------------COMPUTING SOLUTION USING SELECTED LINEAR SOLVERS --------
	void findSolution();
	//---------------------------------------------------------------------

	//-------------- GENEREATE REPORT FROM THE SOLVER PERFORMANCE----------
	void interfaceGenerateReport();
	//----------------------------------------------------------------------

};