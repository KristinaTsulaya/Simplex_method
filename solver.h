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
        size_t nL_rows = 2; // const так как начинаем с двух L строк
    };
    Solver(std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&);
    Data data;
    void artificial_basis_2();
    void chose_col();
    void jordan();
    void chose_row(size_t number_of_l_rows);
    void b_column_to_positive();
    void create_L1_row_and_full();
    void print();
    void matrix_resize(size_t);
    void Simplex_method();
};
#endif //HW_3_SOLVER_H
