#include <stdio.h>
#include <vector>
#include "Matrix.hpp"
#include "LUDecomposition.hpp"

int main(void){

    std::vector<std::vector<double>> array_A = {{1, 2, 3, 4, 5, 6}, {6, 1, 2, 3, 4, 5}, {5, 6, 1, 2, 3, 4}, {4, 5, 6, 1, 2, 3}, {3, 4, 5, 6, 1, 2}, {2, 3, 4, 5, 6, 1}};
    std::vector<double> b = {0, 1, 2, 3, 4, 5};
    
    //Create matrix from array
    Matrix A = Matrix(array_A);
    
    //Create processor of the LUD
    LUDecomposition processor = LUDecomposition(A, b);

    //Apply steps
    processor.decomposite(true);
    processor.apply_FSP(true);
    processor.apply_BSP(true);

    //Get solution
    std::vector<double> x = processor.get_solution();

    return 0;
}