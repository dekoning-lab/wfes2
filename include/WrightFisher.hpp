#pragma once

#include "common.hpp"
#include "rdist.hpp"
#include "NumericVector.hpp"
#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"

namespace WrightFisher {

	enum absorption_type {
		NEITHER = 0,
		EXTINCTION,
		FIXATION,
		BOTH
	};

	enum model_type {
		SINGLE = 0,
		SWITCHING
	};

	class Row {
	public:
		dvec Q;
		llong start;
		llong end;
		llong size;
		double weight;

		explicit Row(llong start, llong end): Q(end - start + 1), start(start), end(end), size(end - start + 1), weight(1) {}
		~Row() {
			std::cout << "Destroying WF row" << std::endl;
		}

		Row(Row&& r): Q(std::move(r.Q)), start(r.start), end(r.end), size(r.size), weight(r.weight) {}

		// static binom_row(llong i, llong Nx, llong Ny, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20);

	};

    struct Matrix {
        llong size;
        absorption_type at;
        llong n_abs;
        smat Q;
        dmat R;

        explicit Matrix(llong size, llong n_abs, absorption_type at): size(size), at(at), n_abs(n_abs), Q(size, size), R(size, n_abs) {
        	std::cout << "Creating WF matrix " << this << std::endl;
        }
        ~Matrix() {
        	std::cout << "Destroying WF matrix " << this << std::endl;
        }

        Matrix(Matrix&& m): size(m.size), at(m.at), n_abs(m.n_abs), Q(std::move(m.Q)), R(std::move(m.R)) {}

    };

	double psi_diploid(llong i, llong N, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9);
	Row binom_row(llong i, llong Nx, llong Ny, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20);
	Matrix Single(llong N, absorption_type mt, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20, llong block_size = 100);

};