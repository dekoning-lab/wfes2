#pragma once

#include "common.hpp"

template <typename T>
class VectorBase {
public:
	llong _stride;
	T* _values;

	llong size;

	explicit VectorBase(): _stride(0), _values(nullptr) {}
	explicit VectorBase(llong size, llong stride, T* values): _stride(stride), _values(values), size(size) {}
	
	explicit VectorBase(llong size, T init = 0, llong stride = 1): _stride(stride), _values(nullptr), size(size) {

		_values = (T*)calloc(size * _stride, sizeof(T));
		assert(_values != NULL);
		if (init != 0) {
			for(llong i = 0; i < (size * _stride); i += _stride) {
				_values[i] = init;
			}
		}
	}

	VectorBase& operator=(const VectorBase& other) {
		if(this == &other) { return *this; }
		assert(size == other.size);

		for(llong i = 0; i < size; i++) {
			_values[i * _stride] = other(i);
		}

		return *this;
	}

	T& operator()(llong i) { return _values[i * _stride]; }
	const T& operator()(llong i) const { return _values[i * _stride]; }

	T sum() const {
		T s = 0;
		for(llong i = 0; i < size; i++) {
			s += (*this)(i);
		}
		return s;
	}

	T normalize() {
		T s = sum();	
		for(long i = 0; i < size; i++) {
			(*this)(i) /= s;
		}
		return s;
	}

	void exp() {
		vdExp(size, _values, _values);
	}

	virtual ~VectorBase() = 0;

};

template <typename T>
VectorBase<T>::~VectorBase() {
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const VectorBase<T>& x) {
	os << x(0);
	for(llong i = 1; i < x.size; i++) {
		os << "\t" << x(i);
	}
	return os;
}

template <typename T>
inline bool operator==(const VectorBase<T>& lhs, const VectorBase<T>& rhs) { 
	bool all = true;
	if(lhs.size != rhs.size) { return false; }
	for(llong i = 0; i < lhs.size; i++) {
		all = (lhs(i) == rhs(i));
	}
	return all;
}

template <typename T>
inline bool operator!=(const VectorBase<T>& lhs, const VectorBase<T>& rhs) { 
	return !(lhs == rhs); 
}
