#pragma once

#include "common.hpp"

class SparseMatrix
{
protected:
    llong current_row;
    bool full;

public:
    llong non_zeros;
    llong n_row;
    llong n_col;

    double* data;
    llong* columns;
    llong* row_index;


    // SparseMatrix(llong current_row, )
    SparseMatrix(llong n_row, llong n_col);

    ~SparseMatrix();

    //       r0 = 1;   r1 = 4;
    // r: a  b  c  d  e  f  g  h 
    //      /        /
    //     /        /
    // m: b  c  d  e
    //    m0 = 0;  m1 = 3;  
    // insert row from r0 to r1 (row index) into matrix from m0 to m1 (matrix index)
    void append_data(dvec& row, llong m0, llong m1, llong r0, llong r1, bool new_row = true, llong slack = 0);
    void finalize_row();

    void debug_print();
    dmat dense();
    dvec multiply(dvec& x, bool transpose = false);
    void multiply_inplace(dvec& x, bool transpose = false);

    // friend ostream& operator<< (ostream& os, const SparseMatrix& M);
    void save_market(std::string path);
    void subtract_identity();
    double operator() (llong i, llong j);
};

std::ostream& operator<<(std::ostream& os, const SparseMatrix& M);

typedef SparseMatrix smat;
