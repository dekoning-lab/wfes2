#pragma once

#include "NumericVector.hpp"

template <typename T>
class DenseMatrix: public MoveOnly {
protected:
		DenseMatrix(llong n_row, llong n_col, T* values): 
			n_row(n_row), n_col(n_col), values(values) {}
public:

	llong n_row;
	llong n_col;

	T* values;

	DenseMatrix(llong n_row, llong n_col, T init = 0): 
		n_row(n_row), n_col(n_col), values(nullptr) {
		
		values = (T*)calloc(n_row * n_col, sizeof(T));
		assert(values != NULL);
		if (init != 0) {
			for(llong i = 0; i < (n_row * n_col); i ++ ) {
				values[i] = init;
			}
		}
	}

	DenseMatrix(DenseMatrix&& r): 
		DenseMatrix<T>(r.n_row, r.n_col, r.values) {
		r.valid = false;
	}

	T& operator()(llong i, llong j) { return values[(i * n_col) + j]; }
	const T& operator()(llong i, llong j) const { return values[(i * n_col) + j]; }

	NumericVectorView<T> row(llong i) {
		return NumericVectorView<T>(n_col, 1, values);
	}

	NumericVectorView<T> col(llong j) {
		return NumericVectorView<T>(n_row, n_col, values + j);
	}

	void print_debug() {
		for(llong i = 0; i < n_row * n_col; i++) {
			std::cout << values[i] << "\t";
		}
		std::cout << std::endl;
	}

	~DenseMatrix() {
		if (this->valid) free(values);	
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const DenseMatrix<T>& x) {
    // os.precision(std::numeric_limits<double>::max_digits10 + 2);
    os.precision(3);
	for(llong j = 0; j < x.n_col; j++) {
		os << x(0, j) << "\t";
	}

	for(llong i = 1; i < x.n_row; i++) {
		os << std::endl;
		for(llong j = 0; j < x.n_col; j++) {
			os << x(i, j) << "\t";
		}
	}
	return os;
}

typedef DenseMatrix<double> dmat;
typedef DenseMatrix<llong> lmat;
