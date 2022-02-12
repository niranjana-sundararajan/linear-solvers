#include <stdio.h>
 #include <string.h>
 #include <iostream>

 #include "Interface.h"
 #include "Matrix.h"
 using namespace std;

 int main(int argc, char** argv) {

     //To access the programme through interface please run this main
     // Alternatively you can access the solvers and the functionalities directly
     // For more examples please refer to README
     
 	 Interface* inter = new Interface();
 	 inter->interfaceSelectInputOptions();
 	 delete inter;


 	return 0;
 }