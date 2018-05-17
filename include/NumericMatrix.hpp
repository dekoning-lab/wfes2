#pragma once

#include "NumericVector.hpp"

template <typename T>
class NumericMatrix : public NumericVectorBase<T> {
public:
	llong rows;
	llong cols;

	NumericMatrix(llong rows, llong cols, T init = 0): NumericVectorBase<T>(rows * cols, init, 1), rows(rows), cols(cols) {}

	T& operator()(llong i) = delete;
	const T& operator()(llong i) const = delete;

	T& operator()(llong i, llong j) { return this->values[(i * rows) + j]; }
	const T& operator()(llong i, llong j) const { return this->values[(i * rows) + j]; }

	NumericVectorView<T> row(llong i) {
		return NumericVectorView<T>(cols, 1, this->values);
	}

	NumericVectorView<T> col(llong j) {
		return NumericVectorView<T>(rows, cols, this->values + j);
	}

	void print_debug() {
		for(llong i = 0; i < rows * cols; i++) {
			std::cout << this->values[i] << "\t";
		}
		std::cout << std::endl;
	}

	~NumericMatrix() {
		if(this->valid) {
			free(this->values);
		}
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const NumericMatrix<T>& x) {
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

typedef NumericMatrix<double> dmat;