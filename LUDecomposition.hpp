#ifndef _DEF_1248729412419283719231231_
#define _DEF_1248729412419283719231231_

#include "Matrix.hpp"
#include <vector>

enum LUDStep{
    INITIALIZED = 1,
    DECOMPOSITED = 2,
    FSP_APPLIED = 3,
    BSP_APPLIED = 4,
};

class LUDecomposition{
    private:
        Matrix matrix_o;
        Matrix matrix_l;
        Matrix matrix_u;
        LUDStep step;
        std::vector<double> vector_b;
        std::vector<double> vector_y;
        std::vector<double> vector_s;
    public:
        LUDecomposition();
        LUDecomposition(Matrix matrix, std::vector<double> vector);
        void decomposite(bool is_debug);
        void apply_FSP(bool is_debug);
        void apply_BSP(bool is_debug);
        std::vector<double> get_solution();
        LUDStep get_step();

};

#endif