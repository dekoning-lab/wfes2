#pragma once

#include "common.hpp"
#include "util.hpp"

class SparseMatrix
{
protected:
    llong current_row;
    bool full;
    llong row_index_start;

public:
    llong non_zeros;
    llong n_row;
    llong n_col;

    double* data;
    llong* cols;
    llong* row_index;
    sparse_matrix_t handle;

    // Empty constructor
    SparseMatrix(llong n_row, llong n_col);
    // Copy constructor
    SparseMatrix(dmat& dense);
    // Identity matrix
    static SparseMatrix LeftPaddedDiagonal(llong dim, double x = 1, llong pad_left = 0);

    ~SparseMatrix();

    //       r0 = 1;   r1 = r0 + size;
    // r: a  b  c  d  e  f  g  h 
    //      /        /
    //     /        /
    // m: b  c  d  e
    //    m0 = 0;  m1 = m0 + size;  
    // insert row from r0 to r1 (row index) into matrix from m0 to m1 (matrix index)

    void append_row(dvec& row, llong col_start, llong size);
    void append_chunk(dvec& row, llong m0, llong r0, llong size);
    void next_row();
    void append_value(double value, llong i, llong j);

    void debug_print();
    dmat dense();
    dvec multiply(dvec& x, bool transpose = false);
    void multiply_inplace_rep(dvec& x, llong times, bool transpose = false);
    SparseMatrix multiply(SparseMatrix& B, bool transpose = false);

    // friend ostream& operator<< (ostream& os, const SparseMatrix& M);
    void save_market(std::string path);
    dvec get_diag_copy();
    void subtract_identity(); // I - Q
    double operator() (llong i, llong j);
    dvec col_copy(llong j);
    dvec row(llong i);

    bool approx_eq(const SparseMatrix& rhs, double tol = 1e-10, bool verbose = false);
private:
};

std::ostream& operator<<(std::ostream& os, const SparseMatrix& M);

typedef SparseMatrix smat;
