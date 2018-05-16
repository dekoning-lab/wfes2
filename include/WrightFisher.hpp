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

    static llong n_abs(absorption_type a_t) {
        switch(a_t) {
            case NEITHER:
                return 0;
            case EXTINCTION:
                return 1;
            case FIXATION:
                return 1;
            case BOTH:
                return 2;
        }
    }

    enum model_type {
        SINGLE = 0,
        SWITCHING
    };

    class Row : public MoveOnly {
    public:
        
        llong start;
        llong end;
        llong size;
        double weight;
        dvec Q;

        Row(llong start, llong end): 
            start(start), end(end), size(end - start + 1), weight(1), Q(end - start + 1) {}

        ~Row() {}

        Row(Row&& r): 
            start(r.start), end(r.end), size(r.size), weight(r.weight), Q(std::move(r.Q)) {}
    };

    struct Matrix : public MoveOnly {
        llong n_row;
        llong n_col;
        absorption_type a_t;
        smat Q;
        dmat R;

        Matrix(llong n_row, llong n_col, absorption_type a_t): 
            n_row(n_row), n_col(n_col), a_t(a_t), Q(n_row, n_col), R(n_row, n_abs(a_t)) {}

        ~Matrix() {}

        Matrix(Matrix&& m): 
            n_row(m.n_row), n_col(m.n_col), a_t(m.a_t), Q(std::move(m.Q)), R(std::move(m.R)) {}

    };

    double psi_diploid(llong i, llong N, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9);
    Row binom_row(llong i, llong Nx, llong Ny, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20);
    Matrix Equilibrium(llong N, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20, llong block_size = 100);
    Matrix Single(llong Nx, llong Ny, absorption_type mt, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9, double alpha = 1e-20, llong block_size = 100);
    // Matrix Switching(lvec Nx, lvec Ny, std::vector<absorption_type>& mt, dvec s, dvec h, dvec u, dvec v, double alpha = 1e-20, llong block_size = 100);

};
