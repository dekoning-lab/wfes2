#pragma once

#include "NumericVectorBase.hpp"

// "Owning" vector, owns its own allocation, cleans up after itself
template <typename T>
class NumericVector : public NumericVectorBase<T> {
public:

	NumericVector(llong size, T init = 0): 
		NumericVectorBase<T>(size, init, 1) {}

	NumericVector(NumericVector&& rhs): 
		NumericVectorBase<T>(rhs.size, rhs.stride, rhs.values) {
			
		rhs.valid = false;
	}

	NumericVector(const std::vector<T>& r): NumericVector<T>(r.size()) {
		for(llong i = 0; i < this->size; i++) {
			(*this)(i) = r[i];
		}
	}

	~NumericVector() {
		if (this->valid) free(this->values);
	}

	static NumericVector zeros(llong size) {
		return NumericVector(size, 0);
	}

	static NumericVector ones(llong size) {
		return NumericVector(size, 1);
	}

	static NumericVector identity(llong size, llong row = 0) {
		NumericVector<T> id(size);
		id(row) = 1;
		return id;
	}

	static NumericVector closed_range(llong start, llong stop) {
		NumericVector<T> r(stop - start + 1);
		for(llong i = start; i <= stop; i++) {
			r(i - start) = (T)i;
		}
		return r;
	}

	void exp_inplace() {
		vdExp(this->size, this->values, this->values);
	}

};

// "Non-owning" vector, points to memory of existing vector, does not clean
template <typename T>
class NumericVectorView : public NumericVectorBase<T>
{
private:

public:

	NumericVectorView(llong size, llong stride, T* values): NumericVectorBase<T>(size, stride, values) {}
	NumericVectorView(NumericVector<T>& d): NumericVectorView(d.size, d.stride, d.values) {}

	~NumericVectorView() {}
};

typedef NumericVector<double> dvec;
typedef NumericVector<llong> lvec;
