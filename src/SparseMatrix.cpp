#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix(llong n_row, llong n_col): current_row(0), full(false), non_zeros(0), n_row(n_row), n_col(n_col) {

    data = (double*)malloc(sizeof(double));
    columns = (llong*)malloc(sizeof(llong));
    row_index = (llong*)malloc((n_row + 1) * sizeof(llong));
    row_index[0] = 0;
}

SparseMatrix::~SparseMatrix() {
    free(data);
    free(columns);
    free(row_index);
}

lvec closed_range(llong start, llong stop) {
    lvec r(stop - start + 1);
    for(llong i = start; i <= stop; i++) r(i - start) = i;
    return r;
}

void SparseMatrix::append_data(dvec& row, llong m0, llong m1, llong r0, llong r1, bool new_row, llong slack) {
    
    assert((m1 - m0) == (r1 - r0));

    llong size = r1 - r0 + 1;
    llong nnz = non_zeros + size + slack;

    // row index
    // if(new_row) row_index[current_row + 1] = nnz;

    // columns
    llong* columns_new = (llong*)realloc(columns, nnz * sizeof(llong));
    if (columns_new != NULL) columns = columns_new;
    else throw std::runtime_error("SparseMatrix::append_data(): Reallocation failed - columns");

    lvec range = closed_range(m0, m1);

    // is accessing data() like this ok?
    memcpy(&columns[non_zeros], range.data(), size * sizeof(llong));

    // data
    double* data_new = (double*)realloc(data, nnz * sizeof(double));
    if (data_new != NULL) data = data_new;
    else throw std::runtime_error("SparseMatrix::append_data(): Reallocation failed - data");

    memcpy(&(data[non_zeros]), &(row.data()[r0]), size * sizeof(double));

    non_zeros += (size + slack);

    if (new_row) finalize_row();
}

void SparseMatrix::finalize_row() {
    row_index[current_row + 1] = non_zeros;
    current_row ++;

    // last row
    if (current_row == n_row) {
        full = true;
        row_index[n_row] = non_zeros;
    }
}

void SparseMatrix::debug_print()
{
    std::cout << "data:    " << std::endl;
    print_buffer(data, (size_t)non_zeros);
    std::cout << "columns:   " << std::endl;
    print_buffer(columns, (size_t)non_zeros);
    std::cout << "row_index:  " << std::endl;
    print_buffer(row_index, (size_t)(n_row + 1));
}

dmat SparseMatrix::dense() {
    dmat dns(n_row, n_col);

    llong info = 0;
    llong* j = (llong*)malloc(6 * sizeof(llong));
    j[0] = 1; j[1] = 0; j[2] = 0;
    j[3] = 2; j[4] = non_zeros; j[5] = 1;

    mkl_ddnscsr(j, &n_row, &n_col, dns.data(), &n_col, data, columns, row_index, &info);

    free(j);

    if(info != 0) throw std::runtime_error("SparseMatrix::dense(): Error processing row " + std::to_string(info));

    return dns;
}

dvec SparseMatrix::multiply(dvec& x, bool transpose) {
    llong v_size = transpose ? n_col : n_row;
    transpose ? assert(x.size() == n_row) : assert(x.size() == n_col);
    dvec y(v_size);

    sparse_matrix_t A;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, n_row, n_col, row_index, row_index + 1, columns, data);
    struct matrix_descr A_descr = {.type = SPARSE_MATRIX_TYPE_GENERAL};
    sparse_operation_t op = transpose ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE;

    mkl_sparse_d_mv(op, 1, A, A_descr, x.data(), 0, y.data());
    
    return y;
}

void SparseMatrix::multiply_inplace(dvec& x, bool transpose) {
    assert(x.size() == n_row); 
    assert(x.size() == n_col);
    transpose ? assert(x.size() == n_row) : assert(x.size() == n_col);
    dvec workspace(x.size());

    sparse_matrix_t A;
    mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, n_row, n_col, row_index, row_index + 1, columns, data);
    struct matrix_descr A_descr = {.type = SPARSE_MATRIX_TYPE_GENERAL};
    sparse_operation_t op = transpose ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE;

    mkl_sparse_d_mv(op, 1, A, A_descr, x.data(), 0, workspace.data());

    for(llong i = 0; i < x.size(); i++) x(i) = workspace(i);
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix& M) {
    os << "%%MatrixMarket matrix coordinate real general" << std::endl;
    os << M.n_row << "\t" << M.n_col << "\t" << M.non_zeros << std::endl;

    os.precision(std::numeric_limits<double>::max_digits10 + 2);
    for (llong i = 0; i < M.n_row; ++i) {
        for (llong j = M.row_index[i]; j < M.row_index[i + 1]; ++j) {
            os << i + 1 << "\t" << M.columns[j] + 1 << "\t" << std::scientific << M.data[j] << std::endl;
        }
    }
    return os;
}

void SparseMatrix::save_market(const std::string path) {
    FILE* out = fopen(path.c_str(), "w");
    fprintf(out, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(out, LPF "\t" LPF "\t" LPF "\n", n_row, n_col, non_zeros);

    for (llong i = 0; i < n_row; ++i) {
        for (llong j = row_index[i]; j < row_index[i + 1]; ++j) {
            fprintf(out, LPF "\t" LPF "\t" DPF "\n", i + 1, columns[j] + 1, data[j]);
        }
    }
    fclose(out);
}

double SparseMatrix::operator() (llong row, llong col) {
    if(row >= current_row) return NAN;
    for(llong j = row_index[row]; j < row_index[row + 1]; j++) {
        if (columns[j] == col) {
            return data[j];
        }
    }
    return 0; // was not found
}

void SparseMatrix::subtract_identity() {
    for (llong i = 0; i < n_row; ++i) {
        for (llong j = row_index[i]; j < row_index[i + 1]; ++j) {
            if (i == columns[j]) data[j] = 1.0 - data[j];
            else data[j] = -data[j];
        }
    }
}
