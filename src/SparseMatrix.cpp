#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix(llong n_row, llong n_col): current_row(0), full(false), non_zeros(0), n_row(n_row), n_col(n_col) {

    values = (double*)malloc(sizeof(double));
    columns = (llong*)malloc(sizeof(llong));
    row_index = (llong*)malloc((n_row + 1) * sizeof(llong));
    row_index[0] = 0;
}

SparseMatrix::~SparseMatrix() {
    if (this->valid) {
        free(values);
        free(columns);
        free(row_index);
    }
}

void SparseMatrix::append_data(dvec& row, llong m0, llong m1, llong r0, llong r1, bool new_row) {
    
    assert((m1 - m0) == (r1 - r0));

    llong size = r1 - r0 + 1;
    llong nnz = non_zeros + size;

    // row index
    if(new_row) row_index[current_row + 1] = nnz;

    // columns
    llong* columns_new = (llong*)realloc(columns, nnz * sizeof(llong));
    if (columns_new != NULL) columns = columns_new;
    else throw runtime_error("SparseMatrix::append_data(): Reallocation failed - columns");

    lvec range = NumericVector<llong>::closed_range(m0, m1);
    memcpy(&columns[non_zeros], range.values, size * sizeof(llong));

    // values
    double* values_new = (double*)realloc(values, nnz * sizeof(double));
    if (values_new != NULL) values = values_new;
    else throw runtime_error("SparseMatrix::append_data(): Reallocation failed - values");

    memcpy(&(values[non_zeros]), &(row.values[r0]), size * sizeof(double));

    non_zeros += size;

    if(new_row) {
        current_row ++;

        // last row
        if(current_row == n_row) {
            full = true;
            row_index[n_row] = non_zeros;
        }
    }
}

void SparseMatrix::debug_print()
{
    cout << "values:    " << endl;
    print_buffer(values, (size_t)non_zeros);
    cout << "columns:   " << endl;
    print_buffer(columns, (size_t)non_zeros);
    cout << "row_index:  " << endl;
    print_buffer(row_index, (size_t)(n_row + 1));
}

DenseMatrix<double> SparseMatrix::dense() {
    DenseMatrix<double> dns(n_row, n_col);

    llong info = 0;
    llong* j = (llong*)malloc(6 * sizeof(llong));
    j[0] = 1; j[1] = 0; j[2] = 0;
    j[3] = 2; j[4] = non_zeros; j[5] = 1;

    mkl_ddnscsr(j, &n_row, &n_col, dns.values, &n_col, values, columns, row_index, &info);

    free(j);

    if(info != 0) throw runtime_error("SparseMatrix::dense(): Error processing row " + to_string(info));

    return dns;
}

dvec SparseMatrix::multiply(dvec& x, bool transpose) {
    dvec y(x.size);
    const char* tr = transpose ? "T" : "N";
    mkl_cspblas_dcsrgemv(tr, &n_row, values, row_index, columns, x.values, y.values);
    return y;
}

void SparseMatrix::multiply_inplace(dvec& x, bool transpose) {
    const char* tr = transpose ? "T" : "N";
    mkl_cspblas_dcsrgemv(tr, &n_row, values, row_index, columns, x.values, x.values);
}

ostream& operator<< (ostream& os, const SparseMatrix& M) {
    os << "%%MatrixMarket matrix coordinate real general" << endl;
    os << M.n_row << "\t" << M.n_col << "\t" << M.non_zeros << endl;

    os.precision(numeric_limits<double>::max_digits10 + 2);
    for (llong i = 0; i < M.n_row; ++i) {
        for (llong j = M.row_index[i]; j < M.row_index[i + 1]; ++j) {
            os << i + 1 << "\t" << M.columns[j] + 1 << "\t" << scientific << M.values[j] << endl;
        }
    }
    return os;
}

void SparseMatrix::save_market(const string path) {
    FILE* out = fopen(path.c_str(), "w");
    fprintf(out, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(out, LPF "\t" LPF "\t" LPF "\n", n_row, n_col, non_zeros);

    for (llong i = 0; i < n_row; ++i) {
        for (llong j = row_index[i]; j < row_index[i + 1]; ++j) {
            fprintf(out, LPF "\t" LPF "\t" DPF "\n", i + 1, columns[j] + 1, values[j]);
        }
    }
    fclose(out);
}

double SparseMatrix::operator() (llong row, llong col) {
    if(row >= current_row) return NAN;
    for(llong j = row_index[row]; j < row_index[row + 1]; j++) {
        if (columns[j] == col) {
            return values[j];
        }
    }
    return 0; // was not found
}

void SparseMatrix::subtract_identity() {
    for (llong i = 0; i < n_row; ++i) {
        for (llong j = row_index[i]; j < row_index[i + 1]; ++j) {
            if (i == columns[j]) values[j] = 1.0 - values[j];
            else values[j] = -values[j];
        }
    }
}
