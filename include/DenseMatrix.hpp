#pragma once

#include "Vector.hpp"

template <typename T>
class DenseMatrix {
public:
	T* _values;
	
	llong rows;
	llong cols;

	explicit DenseMatrix(llong rows, llong cols, T init = 0): 
		_values((T*)calloc(rows * cols, sizeof(T))),
		rows(rows), cols(cols) {

		if (init != 0) {
			for(llong i = 0; i < (rows * cols); i ++ ) {
				_values[i] = init;
			}
		}
	}

	T& operator()(llong i, llong j) { return _values[(i * rows) + j]; }
	const T& operator()(llong i, llong j) const { return _values[(i * rows) + j]; }

	VectorView<T> row(llong i) {
		return VectorView<T>(cols, 1, _values);
	}

	VectorView<T> col(llong j) {
		return VectorView<T>(rows, cols, _values + j);
	}

	void print_debug() {
		for(llong i = 0; i < rows * cols; i++) {
			std::cout << _values[i] << "\t";
		}
		std::cout << std::endl;
	}

	~DenseMatrix() {
		free(_values);
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
