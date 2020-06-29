#ifndef _DEF_12487291957384283719231231_
#define _DEF_12487291957384283719231231_

#include <vector>

class Matrix{
    private:
        std::vector<std::vector<double>> matrix;
    public:
        Matrix();
        Matrix(int i, int j);
        Matrix(std::vector<std::vector<double>> arrays);
        double at(int i, int j);
        void set(int i, int j, double value);
        int dim_r();
        int dim_c();
        bool is_vector();
        bool is_square();
        Matrix zeros(int n);
        Matrix multiply(Matrix mat);
        void show(const char name[]);
};

#endif