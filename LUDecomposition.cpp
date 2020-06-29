#include "LUDecomposition.hpp"

#include "Matrix.hpp"

#include <stdio.h>
#include <vector>

LUDecomposition::LUDecomposition(){}

LUDecomposition::LUDecomposition(Matrix matrix, std::vector<double> vector){
    if(!matrix.is_square() || matrix.dim_r() != vector.size()){
        printf("Unable to initializide LUDecomposition with non square matrix or unmatching vector.\n");
    }
    else{
        matrix_o = matrix;
        matrix_l = matrix.zeros(matrix.dim_r());
        matrix_u = matrix.zeros(matrix.dim_r());
        vector_b = vector;
        vector_y = std::vector<double>(matrix.dim_r());
        vector_s = std::vector<double>(matrix.dim_r());
        step = LUDStep::INITIALIZED;
    }
}
void LUDecomposition::decomposite(bool is_debug){
    if(step != LUDStep::INITIALIZED){
        printf("Already decomposited.\n");
    }
    int n = matrix_o.dim_r();
    for(int k = 1; k <= n - 1; k++){
        for(int i = k + 1; i <= n; i++){
            matrix_o.set(i, k, (matrix_o.at(i, k) / matrix_o.at(k, k)));
            for(int j = k + 1; j <= n; j++){
                matrix_o.set(i, j, (matrix_o.at(i, j) - matrix_o.at(i, k) * matrix_o.at(k, j)));
            }
        }
    }
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            if(i == j){
                matrix_l.set(i, j, 1);
            }
            else if(i > j){
                matrix_l.set(i, j, matrix_o.at(i, j));
            }
            else{
                matrix_l.set(i, j, 0);
            }
        }
    }
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            if(i <= j){
                matrix_u.set(i, j, matrix_o.at(i, j));
            }
            else{
                matrix_u.set(i, j, 0);
            }
        }
    }
    step = LUDStep::DECOMPOSITED;
    if(is_debug){
        matrix_l.show("Matrix L");
        matrix_u.show("Matrix U");
        matrix_l.multiply(matrix_u).show("Matrix LU");
        matrix_o.show("Matrix O");
    }
}
void LUDecomposition::apply_FSP(bool is_debug){
    if(step != LUDStep::DECOMPOSITED){
        printf("Unable to apply FSP before applying steps.\n");
    }
    else{
        int n = matrix_l.dim_r();
        vector_y.at(1 - 1) = vector_b.at(1 - 1);
        for(int k = 2; k <= n; k++){
            double sum = 0;
            for(int i = 1; i <= k - 1; i++){
                sum += matrix_l.at(k, i) * vector_y.at(i - 1);
            }
            vector_y.at(k - 1) = vector_b.at(k - 1) - sum;
        }
        step = LUDStep::FSP_APPLIED;
        if(is_debug){
            printf("================ Vector y ================\n");
            for(int i = 1; i  <= vector_y.size(); i++){
                printf("%1.2e\n", vector_y.at(i - 1));
            }
            printf("================ Vector y ================\n");
        }
    }
}
void LUDecomposition::apply_BSP(bool is_debug){
    if(step != LUDStep::FSP_APPLIED){
        printf("Unable to apply BSP before applying steps.\n");
    }
    else{
        int n = matrix_u.dim_r();
        for(int i = 1; i <= n; i++){
            vector_s.at(i - 1) = 0;
        }
        vector_s.at(n - 1) = (1 / matrix_u.at(n, n)) * vector_y.at(n - 1);
        for(int k = n - 1; k >= 1; k--){
            double sum = 0;
            for(int j = k + 1; j <= n; j++){
                sum += matrix_u.at(k, j) * vector_s.at(j - 1);
            }
            vector_s.at(k - 1) = (1 / matrix_u.at(k, k)) * (vector_y.at(k - 1) - sum);
        }
        step = LUDStep::BSP_APPLIED;
        if(is_debug){
            printf("================ Vector s ================\n");
            for(int i = 1; i <= vector_s.size(); i++){
                printf("%1.2e\n", vector_s.at(i - 1));
            }
            printf("================ Vector s ================\n");
        }
    }
}
std::vector<double> LUDecomposition::get_solution(){
    if(step == LUDStep::BSP_APPLIED){
        return vector_s;
    }
    else{
        printf("Unable to get solution before applying steps.\n");
        return {};
    }
}
LUDStep LUDecomposition::get_step(){
    return step;
}