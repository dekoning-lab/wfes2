#pragma once

#include "VectorBase.hpp"

// "Owning" vector, owns its own allocation, cleans up after itself
template <typename T>
class Vector : public VectorBase<T>
{
public:
	explicit Vector(llong size, T init = 0, llong stride = 1): VectorBase<T>(size, init, stride) {}

	~Vector() {
		free(this->_values);
	}

	static Vector zeros(llong size) {
		return Vector(size, 0);
	}

	static Vector ones(llong size) {
		return Vector(size, 1);
	}

	static Vector identity(llong size, llong row = 0) {
		Vector<T> id(size);
		id(row) = 1;
		return id;
	}

};

// "Non-owning" vector, points to memory of existing vector, does not clean
template <typename T>
class VectorView : public VectorBase<T>
{
private:

public:

	explicit VectorView(llong size, llong stride, T* values): VectorBase<T>(size, stride, values) {}
	explicit VectorView(Vector<T>& d): VectorView(d.size, d.stride, d.values) {}

	~VectorView() {}
};

typedef Vector<double> dvec;