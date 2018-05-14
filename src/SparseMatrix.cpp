#include "SparseMatrix.hpp"

// SparseMatrix::SparseMatrix(
//         llong n_row, llong n_col,
//         vector<double> values, vector<llong> columns, vector<llong> row_index):
//     _current_row(n_row), _full(true),
//     non_zeros((llong)values.size()), n_row(n_row), n_col(n_col)
// {
//     if(n_col + 1 != (llong)row_index.size()) throw runtime_error("SparseMatrix::SparseMatrix(): Size of row_index vector should be n_rows+1");
//     _values = alloc_vector_copy(values);
//     _columns = alloc_vector_copy(columns);
//     _row_index = alloc_vector_copy(row_index);
// }

SparseMatrix::SparseMatrix(llong n_row, llong n_col): _current_row(0), _full(false), non_zeros(0), n_row(n_row), n_col(n_col) {

    _values = (double*)malloc(sizeof(double));
    _columns = (llong*)malloc(sizeof(llong));
    _row_index = (llong*)malloc((n_row + 1) * sizeof(llong));
    _row_index[0] = 0;
}

SparseMatrix::~SparseMatrix() {
    free(_values);
    free(_columns);
    free(_row_index);
}

// closed range [start, stop]
Vector<llong> closed_range(llong start, llong stop) {
    Vector<llong> r(stop - start + 1);
    for(llong i = start; i <= stop; i++) r(i - start) = i;
    return r;
}

void SparseMatrix::append_data(Vector<double>& row, llong m0, llong m1, llong r0, llong r1, bool new_row) {
    
    assert((m1 - m0) == (r1 - r0));

    llong size = r1 - r0 + 1;
    llong nnz = non_zeros + size;

    // row index
    if(new_row) _row_index[_current_row + 1] = nnz;

    // columns
    llong* _columns_new = (llong*)realloc(_columns, nnz * sizeof(llong));
    if (_columns_new != NULL) _columns = _columns_new;
    else throw runtime_error("SparseMatrix::append_data(): Reallocation failed - columns");

    Vector<llong> range = closed_range(m0, m1);
    cout << "range: "; print_buffer(range._values, range.size);
    cout << m0 << ":" << m1 << "; " << r0 << ":" << r1 << endl;
    memcpy(&_columns[non_zeros], range._values, size * sizeof(llong));
    print_buffer(_columns, nnz);

    // values
    double* _values_new = (double*)realloc(_values, nnz * sizeof(double));
    if (_values_new != NULL) _values = _values_new;
    else throw runtime_error("SparseMatrix::append_data(): Reallocation failed - values");

    memcpy(&(_values[non_zeros]), &(row._values[r0]), size * sizeof(double));

    non_zeros += size;

    if(new_row) {
        _current_row ++;

        // last row
        if(_current_row == n_row) {
            _full = true;
            _row_index[n_row] = non_zeros;
        }
    }
}

void SparseMatrix::debug_print()
{
    cout << "values:    " << endl;
    print_buffer(_values, (size_t)non_zeros);
    cout << "columns:   " << endl;
    print_buffer(_columns, (size_t)non_zeros);
    cout << "row_index:  " << endl;
    print_buffer(_row_index, (size_t)(n_row + 1));
}

DenseMatrix<double> SparseMatrix::dense() {
    DenseMatrix<double> dns(n_row, n_col);

    llong info = 0;
    llong* j = (llong*)malloc(6 * sizeof(llong));
    j[0] = 1; j[1] = 0; j[2] = 0;
    j[3] = 2; j[4] = non_zeros; j[5] = 1;

    mkl_ddnscsr(j, &n_row, &n_col, dns._values, &n_col, _values, _columns, _row_index, &info);

    free(j);

    if(info != 0) throw runtime_error("SparseMatrix::dense(): Error processing line " + to_string(info));

    return dns;
}

Vector<double> SparseMatrix::multiply(Vector<double> x) {
    Vector<double> y(x.size);
    mkl_cspblas_dcsrgemv("N", &n_row, _values, _row_index, _columns, x._values, y._values);
    return y;
}

ostream& operator<< (ostream& os, const SparseMatrix& M) {
    os << "%%MatrixMarket matrix coordinate real general" << endl;
    os << M.n_row << "\t" << M.n_col << "\t" << M.non_zeros << endl;

    cout.precision(numeric_limits<double>::max_digits10 + 2);
    for (llong i = 0; i < M.n_row; ++i) {
        for (llong j = M._row_index[i]; j < M._row_index[i + 1]; ++j) {
            os << i + 1 << "\t" << M._columns[j] + 1 << "\t" << scientific << M._values[j] << endl;
        }
    }
    return os;
}

void SparseMatrix::save_market(const string path) {
    FILE* out = fopen(path.c_str(), "w");
    fprintf(out, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(out, LPF "\t" LPF "\t" LPF "\n", n_row, n_col, non_zeros);

    for (llong i = 0; i < n_row; ++i) {
        for (llong j = _row_index[i]; j < _row_index[i + 1]; ++j) {
            fprintf(out, LPF "\t" LPF "\t" DPF "\n", i + 1, _columns[j] + 1, _values[j]);
        }
    }
    fclose(out);
}

double SparseMatrix::operator() (llong row, llong col) {
    if(row >= _current_row) return NAN;
    for(llong j = _row_index[row]; j < _row_index[row + 1]; j++) {
        if (_columns[j] == col) {
            return _values[j];
        }
    }
    return 0; // was not found
}

void SparseMatrix::subtract_identity() {
    for (llong i = 0; i < n_row; ++i) {
        for (llong j = _row_index[i]; j < _row_index[i + 1]; ++j) {
            if (i == _columns[j]) _values[j] = 1.0 - _values[j];
            else _values[j] = -_values[j];
        }
    }
}
