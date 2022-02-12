#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>

#include "../include/Matrix.h"
#include "../src/Matrix.cpp"
#include "../include/CSRMatrix.h"
#include "../src/CSRMatrix.cpp"

using namespace std;



int main(){
    // Set values for dense matrix
    int rows = 5;
    int cols = 5;
    auto *dense_mat = new Matrix<double>(rows, cols,true);

    // Set values for dense matrix
    vector<double> vs = {6, 0, 1, 1, 0,1, 0, 0, 0, 1, 1, 0, 5, 1, 0, 1, 0, 0, 7, 7, 0, 0, 0, 7, 1};
    
    // store values for dense matrices
    for (int i = 0; i < rows*cols; i++){
        dense_mat->values[i] = vs[i];
    }
    // print values for dense matrices
    // dense_mat->printMatrix();
    
    // set vales for sparse matrices
    int nnzs = 12;
    int srows = 5;
    int scols = 5;
    
    // create a sparse matrix
    auto* A = new CSRMatrix<double>(srows, scols, nnzs, true);

    // convert dense matrix to sparse matrix
    A->readSparse(*dense_mat, A);

    // print sparse matrix
    A->printMatrix(); 

    // initialize values of B and X
    auto *B = new Matrix<double>(rows, 1 , true);
    auto *X = new Matrix<double>(rows, 1, true);

    // set values for B and X
    for(int i = 0; i < rows; i++){
        B->values[i] = i+1;
        X->values[i] = 0;
    }

    // print your linear system
    cout<< " Matrix A" << endl;
    A->printValues();
    cout << " Matrix B" << endl;
    B->printMatrix();
    cout << " Matrix X" << endl;
    X->printMatrix();
    
    A->jacobiSparse(*A, *X, *B);
    delete dense_mat;
    delete A;

}
