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

	enum model_type {
		SINGLE = 0,
		SWITCHING
	};

	struct Row {
		dvec Q;
		llong start;
		llong end;
		llong size;
		double weight;

		explicit Row(): Q(), start(0), end(0), size(0), weigth(0) {}
		explicit Row(llong start, llong end, llong size): Q(size), start(start), end(end), size(size), weight(1) {}
		explicit Row(llong i, llong Nx, llong Ny, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20);
		~Row() {
			std::cout << "Destroying WF row" << std::endl;
		}
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

    };

	double psi_diploid(llong i, llong N, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9);
	// Row binom_row(llong i, llong Nx, llong Ny, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20);
	static Matrix Single(llong N, absorption_type mt, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20, llong block_size = 100);

};