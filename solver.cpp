#include "solver.h"

Solver::Solver(std::vector<std::vector<double>> &table_, std::vector<double> &L_, std::vector<double> &b_) : L(L_),
                                                                                                             b(b_) {

    rows = table_.size() + 1;
    cols = table_[0].size() + 1;

    simplex_table.resize(rows, std::vector<double>());
    for (size_t i = 0; i < simplex_table.capacity(); ++i) {
        simplex_table[i].resize(cols);
    }

    for (size_t i = 0; i < rows - 1; ++i) {
        simplex_table[i][cols - 1] = b[i];
        for (size_t j = 0; j < cols - 1; ++j) {
            simplex_table[i][j] = table_[i][j];
        }
    }

    simplex_table[rows - 1][cols - 1] = 0.0;
    for (size_t j = 0; j < cols - 1; ++j) {
        simplex_table[rows - 1][j] = L[j] * (-1);
    }
}

void Solver::b_column_to_positive() { // 1 step

    for (size_t i = 0; i < rows - 1; ++i) {
        if (simplex_table[i][cols - 1] < 0) {
            for (size_t j = 0; j < cols; ++j) {
                simplex_table[i][j] *= (-1);
            }
        }
    }

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << simplex_table[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Solver::create_L1_row_and_full() { // 2 step

    std::vector<std::vector<double>> simplex_table_copy;
    simplex_table_copy.resize(rows + 1, std::vector<double>());
    for (size_t i = 0; i < rows + 1; ++i) {
        simplex_table_copy[i].resize(cols);
    }

    ++rows;

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (i != rows - 1) {
                simplex_table_copy[i][j] = simplex_table[i][j];
            }
        }
    }

    for (size_t j = 0; j < cols; ++j) { // идем по столбцам
        for (size_t i = 0; i < rows - 2; ++i) { // строки до L
            simplex_table_copy[rows - 1][j] += simplex_table[i][j];
        }
        simplex_table_copy[rows - 1][j] *= (-1);
    }
    simplex_table = simplex_table_copy;
    for (const auto& i: simplex_table) {
        for (auto j: i) {
            std::cout << j << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::vector<std::vector<double>>().swap(simplex_table_copy); // free memory
}

void Solver::find_chose_col(){
    for (size_t j = 0; j < cols - 1; ++j) {
        if (simplex_table[rows - 1][j] <= data.min_elem_L1) {
            data.chose_col = j;
        }
        data.min_elem_L1 = std::min(data.min_elem_L1, simplex_table[rows - 1][j]);
    }
}
void Solver::find_chose_row(size_t number_of_l_rows){

    for (size_t k = 0; k < rows - number_of_l_rows; ++k) { // так как не нужно делить на элементы строки L L1 // todo 2
        double cur_b_elem = simplex_table[k][cols - 1];
        for (size_t i = 0; i < rows - number_of_l_rows; ++i) {
            for (size_t j = 0; j < cols - 1; ++j) {
                if (simplex_table[i][j] != 0) {
                    data.cur_koef = cur_b_elem / simplex_table[i][j];
                    if (data.cur_koef >= 0) { // обязательное условие для них
                        data.min_koef = std::min(data.min_koef, data.cur_koef);
                        if (data.prev_min_koef >
                            data.min_koef) { // если на данном шаге предыдущее значение коэфа больше чем на этом, то нам предыдузий тогда нафиг не нужен
                            data.chose_row = i; // в моем кейсе вторая строка или 1 по системе отсчета с нуля //
                        }
                        data.prev_min_koef = data.min_koef;
                    }
                }
            }
        }
    }
}

void Solver::artificial_basis_2(size_t chose_row_, size_t chose_col_) {

    bool do_basis_change = true;

    size_t counter_of_artificial_vars = 0;
    size_t nL_rows = 2; // const так как начинаем с двух L строк

    while (do_basis_change) {

        ++counter_of_artificial_vars;

        find_chose_col();


//        if(do_basis_change) {
//            find_chose_row(nL_rows);
//        } else if(do_optimum_decision) {
//            find_chose_row(nL_rows);
//        }




        for (size_t k = 0; k < rows - nL_rows; ++k) { // так как не нужно делить на элементы строки L
            double cur_b_elem = simplex_table[k][cols - 1];
            for (size_t i = 0; i < rows - nL_rows; ++i) {
                for (size_t j = 0; j < cols - 1; ++j) {
                    if (simplex_table[i][j] != 0) {
                        data.cur_koef = cur_b_elem / simplex_table[i][j];
                        if (data.cur_koef >= 0) { // обязательное условие для них
                            data.min_koef = std::min(data.min_koef, data.cur_koef);
                            if (data.prev_min_koef > data.min_koef) { // если на данном шаге предыдущее значение коэфа больше чем на этом, то нам предыдузий тогда нафиг не нужен
                                data.chose_row = i; // в моем кейсе вторая строка или 1 по системе отсчета с нуля //
                            }
                            data.prev_min_koef = data.min_koef;
                        }
                    }
                }
            }
        }

        double bas_val = simplex_table[data.chose_row][data.chose_col];
        std::cout << bas_val << std::endl;
        simplex_table[data.chose_row][data.chose_col] = 1;

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                if (i != data.chose_row && j != data.chose_col) {
                    simplex_table[i][j] =
                            simplex_table[i][j] * bas_val -
                            simplex_table[data.chose_row][j] * simplex_table[i][data.chose_col];
                }
            }
        }

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                simplex_table[i][j] /= bas_val;
            }
        }

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << std::setprecision(3) << simplex_table[i][j] << std::fixed << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        data.min_elem_L1 = 0;
        data.min_koef = 100000;

        if (simplex_table[rows - 1][cols - 1] == 0) {
            do_basis_change = false;
            std::cout << "Artificial variables are derived (^_^)" << std::endl;
        } else if (simplex_table[rows - 1][cols - 1] < 0) {
            do_basis_change = true;
        }
    }


    ////////////////////////////////////////
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
                std::cout << std::setprecision(3) << simplex_table[i][j] << std::fixed << " ";
            }
        std::cout << std::endl;
    } std::cout << std::endl;


    std::vector<std::vector<double>> simplex_table_copy;
    simplex_table_copy.resize(rows - 1, std::vector<double>());
    for (size_t i = 0; i < rows; ++i) {
        simplex_table_copy[i].resize(cols - counter_of_artificial_vars);
    }
    for (size_t i = 0; i < simplex_table.size() - 1; ++i) {
        for (size_t j = 0; j < simplex_table[i].size(); ++j) {
            if (j < 2) {
                simplex_table_copy[i][j] = simplex_table[i][j];
            } else if (j == cols - 1) {
                simplex_table_copy[i][cols - counter_of_artificial_vars - 1] = simplex_table[i][j];
            }
        }
    }

    rows = simplex_table_copy.size();
    cols = simplex_table_copy[0].size();

    simplex_table = simplex_table_copy;
    std::vector<std::vector<double>>().swap(simplex_table_copy);
////////////////////////////////////

    --nL_rows;
    bool do_optimum_decision = true;
    while (do_optimum_decision) {

    find_chose_col();

        for (size_t k = 0; k < rows - nL_rows; ++k) { // так как не нужно делить на элементы строки L L1
            double cur_b_elem = simplex_table[k][cols - 1];
            for (size_t i = 0; i < rows - nL_rows; ++i) {
                for (size_t j = 0; j < cols - 1; ++j) {
                    if (simplex_table[i][j] != 0) {
                        data.cur_koef = cur_b_elem / simplex_table[i][j];
                        if (data.cur_koef >= 0) { // обязательное условие для них
                            data.min_koef = std::min(data.min_koef, data.cur_koef);
                            if (data.prev_min_koef > data.min_koef) { // если на данном шаге предыдущее значение коэфа больше чем на этом, то нам предыдузий тогда нафиг не нужен
                                data.chose_row = i; // в моем кейсе вторая строка или 1 по системе отсчета с нуля //
                            }
                            data.prev_min_koef = data.min_koef;
                        }
                    }
                }
            }
        }

        double bas_val = simplex_table[data.chose_row][data.chose_col];
        std::cout << bas_val << std::endl;
        simplex_table[data.chose_row][data.chose_col] = 1;

        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                if (i != data.chose_row && j != data.chose_col) {
                    simplex_table[i][j] =
                            simplex_table[i][j] * bas_val -
                            simplex_table[data.chose_row][j] * simplex_table[i][data.chose_col];
                }
            }
        }

        for (auto &i: simplex_table) {
            for (double &j: i) {
                j /= bas_val;
            }
        }

        for (size_t i = 0; i < rows; ++i) {
            if (i != data.chose_row) {
                simplex_table[i][data.chose_col] *= (-1);
            }
        }

        for (auto &i: simplex_table) {
            for (double j: i) {
                std::cout << std::setprecision(3) << j << std::fixed << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        data.min_elem_L1 = 0;
        data.min_koef = 100000;


        for (size_t j = 0; j < cols; ++j) {
            if (simplex_table[rows - 1][j] >= 0) {
                do_optimum_decision = false;
            } else {
                do_optimum_decision = true;
                break;
            }
        }
    }
}