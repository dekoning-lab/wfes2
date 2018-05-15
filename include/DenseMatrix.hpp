#pragma once

#include "NumericVector.hpp"

template <typename T>
class DenseMatrix: public MoveOnly {
protected:
		DenseMatrix(llong rows, llong cols, T* values): 
			rows(rows), cols(cols), values(values) {}
public:

	llong rows;
	llong cols;

	T* values;

	DenseMatrix(llong rows, llong cols, T init = 0): 
		rows(rows), cols(cols), values(nullptr) {
		
		values = (T*)calloc(rows * cols, sizeof(T));
		assert(values != NULL);
		if (init != 0) {
			for(llong i = 0; i < (rows * cols); i ++ ) {
				values[i] = init;
			}
		}
	}

	DenseMatrix(DenseMatrix&& r): 
		DenseMatrix<T>(r.rows, r.cols, r.values) {
		r.valid = false;
	}

	T& operator()(llong i, llong j) { return values[(i * rows) + j]; }
	const T& operator()(llong i, llong j) const { return values[(i * rows) + j]; }

	NumericVectorView<T> row(llong i) {
		return NumericVectorView<T>(cols, 1, values);
	}

	NumericVectorView<T> col(llong j) {
		return NumericVectorView<T>(rows, cols, values + j);
	}

	void print_debug() {
		for(llong i = 0; i < rows * cols; i++) {
			std::cout << values[i] << "\t";
		}
		std::cout << std::endl;
	}

	~DenseMatrix() {
		if (this->valid) {
			free(values);	
		}
		
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const DenseMatrix<T>& x) {
	for(llong j = 0; j < x.cols; j++) {
		os << x(0, j) << "\t";
	}

	for(llong i = 1; i < x.rows; i++) {
		os << std::endl;
		for(llong j = 0; j < x.cols; j++) {
			os << x(i, j) << "\t";
		}
	}
	return os;
}

typedef DenseMatrix<double> dmat;
typedef DenseMatrix<llong> lmat;
