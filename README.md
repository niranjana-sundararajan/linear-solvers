# Linear Solvers for Dense and Sparse Matrix Systems
![alt text](https://www.brighton.ac.uk/images/School-of-Computing-Engineering-and-Mathematics/Artificial-intelligence-banner-Cropped-990x259.jpg?f=webp&q=80&w=1300) \
This software library can be used to solve a linear system AX = B, a descretized version of many PDEs. Where A is a posititive definite matrix of a defined size, and X and B are vectors of different sizes. The library accounts for dense and sparse matrices as inputs to the linear system. Contained witthin it are different solvers that can be used based on type of input matrices and considerations for memory and time optimizations. Read our report linked [here](report/report.pdf) for further details.


## Solvers in this Library
**Dense Direct**
- Gaussian Elimination
- LU Decomposition
- Cholesky Decomposition 

**Dense Iterative**
- Gauss Siedel
- Jacobi
- Conjugate gradient 

**Sparse** 
- Jacobi
- Conjugate Gradient

## Setup Instructions

### Windows Recommendations
C++ : Version 11/17 \
IDE : Visual Studio Community/ Visual Studio Code \
Compiler : g++ (v 13.0)

### Installation Instructions
Download/clone the github repository into your local IDE of choice. \
Open the project/folder and follow the example below- replace the solver given with the solver of your choice 

### Using the Library through the Terminal Interface

Step 1 : Select 1 to read Matrix Data from our Demo Text Files or input your own Matrix in a .txt format delimited by commas \
<img width="214" alt="interfaceterminal" src="https://user-images.githubusercontent.com/88569855/151702774-72e46d88-b897-4fe0-83fa-751c93f55151.png"> \
-- Assuming you have chosen to read in from our Demo files ---\
Step 2:  Select the type of Solver you want to test - Dense or Sparse\
<img width="221" alt="interfaceterminal2" src="https://user-images.githubusercontent.com/88569855/151702876-bae9138d-12ab-408b-9dc5-4c46d1d4c2e5.png">\
Step 3 : Select a Matrix Size to test on \
<img width="204" alt="interfaceterminal3" src="https://user-images.githubusercontent.com/88569855/151702913-5bd57a09-3e6a-48a2-96c9-0d335bcff9c2.png">\
Step 4 :  Verify whether matrix has been read-in from the text file. For the purposes of this demo we have chosen a 4x4 matrix to diplay the read-ins to terminal \
<img width="206" alt="interfaceterminal4" src="https://user-images.githubusercontent.com/88569855/151702932-150f8184-30a2-4c82-a1a9-719e7c7849ce.png">\
Step 5 : Choose a Solver and generate the output\
<img width="205" alt="interfaceterminal5" src="https://user-images.githubusercontent.com/88569855/151702993-adf296db-1ac2-4d95-aded-e6a511132760.png">


### Compiling The Files
(To run the terminal tests )In the runInterface folder please execute :\
```
g++ " Matrix.h Matrix.cpp CSRMatrix.h CSRMAtrix.cpp Solver.h Test.h main.cpp "
```

##  Examples of Use as a Library
For examples of implementation of individual solvers, check the examples folder in the github repository) \

```
#include "../include/Solver.h"
#include "../src/Solver.cpp"
#include "../include/ReadFile.h"
#include "../src/ReadFile.cpp"

using namespace std;


int main()
{
	int size = 4;

	// Initialize martices with pointers

	auto* A = new Matrix<double>(size, size, true);
	auto* X = new Matrix<double>(size, 1, true);
	auto* B = new Matrix<double>(size, 1, true);

	// Read Matrix from text file
	readMatrix("../tests/mat4_A.txt", *A );
	readMatrix("../tests/mat4_B.txt", *B );

	// Print the read matrix to terminal (if needed)
	A->printMatrix();
	B->printMatrix();


    // set values of output matrix to zero
	for (int i = 0; i < size; i++){
		X->values[i] = 0;
	}

    // Use solver of choice by creating a solver object
	Solver<double> gauss;
	gauss.gaussElimination(*A, *B, *X);
   
   
	writeMatrix("gauss_output_4x4.txt", *X);

    delete[] A;
    delete[] X;
    delete[] B;

}

```


### Results of Some Testing on a Solvers and their performances


<img width="301" alt="Times for Solvers" src="https://user-images.githubusercontent.com/88569855/151703062-d600de80-4937-4b34-b573-c251cea67023.png"> 

For detailed review on analyis and code design and structure, please refer to the linked [report](report/report.pdf)

## Version History

[Version Files](version.md)
## License

This project is licensed under the MIT License  - see the [LICENSE.md](LICENSE) file for details
