#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>

#include "../include/Matrix.h"
#include "../include/Solver.h"

using namespace std;
//template <class T>
int main()
{
	int rows = 4;
	int cols = 4;

	double a[] = { 4.4376643, 1.6318403, -0.0708927, -1.8417670, -1.1942025, 0.2047227, -1.6580309, 0.7796127, -2.8642194, -1.7816485, 17.3547046, 0.9597557, -4.2666389, -0.2932443, -2.9780207, -14.4609881 };
	double x[] = { 0, 0, 0, 0 };
	double b[] = { -4.1492027, -4.5481267, -3.0348628, -2.4409358 };

	auto* A = new Matrix<double>(rows, cols, true);
	auto* X = new Matrix<double>(rows, cols, true);
	auto* B = new Matrix<double>(rows, cols, true);

	for (int i = 0; i < rows * cols; i++)
	{
		A->values[i] = a[i];
	}
	for (int i = 0; i < rows; i++)
	{
		B->values[i] = b[i];
		X->values[i] = x[i];
	}

	Solver<double> Sol;
	Sol.gaussSeidel(*A, *B, *X);
}
