#pragma once

#include "common.hpp"
#include "rdist.hpp"
#include "Vector.hpp"
#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"

namespace WrightFisher {

	enum absorption_type {
		NEITHER = 0,
		EXTINCTION,
		FIXATION,
		BOTH
	};

	struct Row {
		dvec Q;
		llong start;
		llong end;
		llong size;
		double weight;

		Row(llong start, llong end, llong size): Q(size), start(start), end(end), size(size), weight(1) {}
	};

    struct Matrix {
        llong size;
        absorption_type mt;
        llong n_abs;
        smat Q;
        dmat R;
        Matrix(llong size, llong n_abs, absorption_type mt): size(size), mt(mt), n_abs(n_abs), Q(size, size), R(size, n_abs) {}
    };

	double psi_diploid(llong i, llong N, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9);
	Row binom_row(llong i, llong Nx, llong Ny, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20);

	Matrix Single(llong N, absorption_type mt, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20, llong block_size = 100);
};