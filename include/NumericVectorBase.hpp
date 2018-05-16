#pragma once

#include "common.hpp"

template <typename T>
class NumericVectorBase : public MoveOnly {
protected:	
	NumericVectorBase(llong size, llong stride, T* values): 
		stride(stride), size(size), values(values) {}

public:
	llong stride;
	llong size;
	T* values;

	NumericVectorBase(llong size, T init = 0, llong stride = 1): 
		stride(stride), size(size), values(nullptr) {

		values = (T*)calloc(size * stride, sizeof(T));
		assert(values != NULL);
		if (init != 0) {
			for(llong i = 0; i < (size * stride); i += stride) {
				values[i] = init;
			}
		}
	}

	T& operator()(llong i) { 
		return values[i * stride]; 
	}

	const T& operator()(llong i) const { 
		return values[i * stride]; 
	}

	T sum() const {
		T s = 0;
		for(llong i = 0; i < size; i++) {
			s += (*this)(i);
		}
		return s;
	}

	void operator*=(T v) {
		for(llong i = 0; i < size; i++) {
			(*this)(i) *= v;
		}
	}

	T normalize_inplace() {
		T s = sum();	
		for(llong i = 0; i < size; i++) {
			(*this)(i) /= s;
		}
		return s;
	}

	void negate_inplace() {
		for(llong i = 0; i < size; i++) {
			T x = (*this)(i);
			(*this)(i) = -x;
		}
	}

	void abs_inplace() {
		for(llong i = 0; i < size; i++) {
			T x = (*this)(i);
			(*this)(i) = fabs(x);
		}
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const NumericVectorBase<T>& x) {
	os << x(0);
	for(llong i = 1; i < x.size; i++) {
		os << "\t" << x(i);
	}
	return os;
}

template <typename T>
inline bool operator==(const NumericVectorBase<T>& lhs, const NumericVectorBase<T>& rhs) { 
	bool all = false;
	if(lhs.size != rhs.size) { return false; }
	for(llong i = 0; i < lhs.size; i++) {
		all = (lhs(i) == rhs(i));
	}
	return all;
}

template <typename T>
inline bool operator!=(const NumericVectorBase<T>& lhs, const NumericVectorBase<T>& rhs) { 
	return !(lhs == rhs); 
}
