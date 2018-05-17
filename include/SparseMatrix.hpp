#pragma once

#include "common.hpp"
#include "NumericVector.hpp"
#include "NumericMatrix.hpp"

class SparseMatrix
{
protected:
    mutable bool valid = true;
    llong current_row;
    bool full;

    SparseMatrix(llong current_row, bool full, llong non_zeros, llong n_row, llong n_col, double* values, llong* columns, llong* row_index):
        current_row(current_row), full(full), non_zeros(non_zeros), n_row(n_row), n_col(n_col), values(values), columns(columns), row_index(row_index) {}

public:



    llong non_zeros;
    llong n_row;
    llong n_col;

    double* values;
    llong* columns;
    llong* row_index;


    // SparseMatrix(llong current_row, )
    SparseMatrix(llong n_row, llong n_col);

    ~SparseMatrix();

    SparseMatrix(const SparseMatrix& r): SparseMatrix(r.current_row, r.full, r.non_zeros, r.n_row, r.n_col, r.values, r.columns, r.row_index) {
        r.valid = false;
    }

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
