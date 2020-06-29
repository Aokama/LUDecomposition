#include "Matrix.hpp"

#include <vector>
#include <stdio.h>

Matrix::Matrix(){}

Matrix::Matrix(int i, int j){
    matrix = std::vector<std::vector<double>>(i, std::vector<double>(j));
}
Matrix::Matrix(std::vector<std::vector<double>> arrays){
    matrix = arrays;
}
double Matrix::at(int i, int j){
    if(i <= dim_r() && j <= dim_c()){
        return matrix.at(i - 1).at(j - 1);
    }
    else{
        printf("(%d, %d) is index out of bound. returned (0.0).\n", i, j);
        return 0;
    }
}
void Matrix::set(int i, int j, double value){
    if(i <= dim_r() && j <= dim_c()){
        matrix[i - 1][j - 1] = value;
    }
    else{
        printf("(%d, %d) is index out of bound. set (0.0).\n", i, j);
        matrix[i - 1][j - 1] = 0;
    }
}
int Matrix::dim_r(){
    return matrix.size();
}
int Matrix::dim_c(){
    return matrix.at(0).size();
}
Matrix Matrix::zeros(int n){
    Matrix mat = Matrix(n, n);
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            mat.set(i, j, 0);
        }
    }
    return mat;
}
bool Matrix::is_vector(){
    return dim_c() == 1 ? true : false;
}
bool Matrix::is_square(){
    return (dim_r() == dim_c() && dim_r() != 0) ? true : false;
}
Matrix Matrix::multiply(Matrix mat){
    if(dim_c() != mat.dim_r()){
        printf("Unable to multiply in this dimensions.\n");
        return Matrix();
    }
    else{
        Matrix res = Matrix(dim_r(), mat.dim_c());
        for(int i = 1; i <= dim_r(); i++){
            for(int j = 1; j <= mat.dim_c(); j++){
                double sum = 0;
                for(int k = 1; k <= dim_c(); k++){
                    sum += at(i, k) * mat.at(k, j);
                }
                res.set(i, j, sum);
            }
        }
        return res;
    }
}
void Matrix::show(const char name[]){
    printf("================ %s ================\n", name);
    for(int i = 1; i <= dim_r(); i++){
        for(int j = 1; j <= dim_c(); j++){
            printf("%1.2e\t", matrix.at(i - 1).at(j - 1));
        }
        printf("\n");
    }
    printf("================ %s ================\n", name);
}