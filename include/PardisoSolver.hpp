#pragma once

#include "MKL_Const.hpp"
#include "SparseMatrix.hpp"

using namespace std;

class PardisoSolver {

    llong size;
    llong n_right_hand_sides;
    llong max_factors;
    llong matrix_type;
    llong matrix_number;
    llong error;
    llong message_level;
    llong phase;

    SparseMatrix& m;
    lvec control;
    lvec internal;
    dvec workspace;


public:
    PardisoSolver(SparseMatrix& A, llong matrix_type, llong message_level, llong n_rhs = 1);
    ~PardisoSolver();
    void analyze();
    dvec solve(dvec& b, bool transpose = false);
    dvec get_diagonal();
};
