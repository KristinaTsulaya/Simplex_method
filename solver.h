//
// Created by tsulayakris on 11.05.23.
//

#ifndef HW_3_SOLVER_H
#define HW_3_SOLVER_H
#include <iostream>
#include <vector>
#include <iomanip>

class Solver {
private:
    std::vector<std::vector<double>> simplex_table;
    size_t rows;
    size_t cols;
    std::vector<double> L, b;

public:
    struct Data{
        double min_koef = 1000000; // специально задаю любой неотрицательный коэф, который чисто теор может быть в этой таблицу минимальным
        double prev_min_koef = 1000000;
        double cur_koef = 0;
        double min_elem_L1 = 0; // по сути максимальный из возможных, так как на него огранич (<=0)
        size_t chose_col = 0;
        size_t chose_row = 0;
    };
    Solver(std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&);
    Data data;
    void artificial_basis_2(size_t, size_t);
    void find_chose_col();
    void find_chose_row(size_t);
    void b_column_to_positive();
    void create_L1_row_and_full();
};

#endif //HW_3_SOLVER_H
