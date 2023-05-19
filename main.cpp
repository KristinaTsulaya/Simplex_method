#include <iostream>
#include "solver.h"
int main() {

    std::vector<std::vector<double>> table = {{1, 5, 6, 3, 0},
                                              {2, -1, 2, -8, 1},
                                              {-3, 6, -3, 4, -7}};

    std::vector<double> L = {-4, 3, -9, 1, 5};
    std::vector<double> b = {15, -4, -3};

    Solver solver(table, L, b);
    solver.b_column_to_positive();
    solver.create_L1_row_and_full();
    solver.artificial_basis_2();
    solver.Simplex_method();

    return 0;
}
