#pragma once

#include "common.hpp"
#include "Vector.hpp"
#include "DenseMatrix.hpp"

using namespace std;

class SparseMatrix
{
public:

    double* _values;
    llong* _columns;
    llong* _row_index;

    llong _current_row;
    bool _full;

    llong non_zeros;
    llong n_row;
    llong n_col;


    SparseMatrix(llong n_row, llong n_col);

    ~SparseMatrix();

    //       r0 = 1;   r1 = 4;
    // r: a  b  c  d  e  f  g  h 
    //      /        /
    //     /        /
    // m: b  c  d  e
    //    m0 = 0;  m1 = 3;  
    // insert row from r0 to r1 (row index) into matrix from m0 to m1 (matrix index)
    void append_data(Vector<double>& row, llong m0, llong m1, llong r0, llong r1, bool new_row = true);

    void debug_print();
    DenseMatrix<double> dense();
    Vector<double> multiply(Vector<double> x);

    friend ostream& operator<< (ostream& os, const SparseMatrix& M);
    void save_market(string path);
    void subtract_identity();
    double operator() (llong i, llong j);
};

typedef SparseMatrix smat;
