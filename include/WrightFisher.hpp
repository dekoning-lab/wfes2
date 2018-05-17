#pragma once

#include "common.hpp"
#include "rdist.hpp"
#include "SparseMatrix.hpp"

namespace WrightFisher {

    enum absorption_type {
        NON_ABSORBING = 0,
        EXTINCTION_ONLY,
        FIXATION_ONLY,
        BOTH_ABSORBING
    };

    inline llong n_absorbing(absorption_type a_t) {
        switch(a_t) {
            case NON_ABSORBING:
                return 0;
            case EXTINCTION_ONLY:
                return 1;
            case FIXATION_ONLY:
                return 1;
            case BOTH_ABSORBING:
                return 2;
            default:
                throw std::runtime_error("Unknown absorption type");
        }
    }

    class Row {
    protected:
        mutable bool valid = true;
    public:
        
        llong start;
        llong end;
        llong size;
        double weight;
        dvec Q;

        Row(llong start, llong end): start(start), end(end), size(end - start + 1), weight(1), Q(end - start + 1) {}

        Row(const Row& r): start(r.start), end(r.end), size(r.size), weight(r.weight), Q(r.Q) {
            r.valid = false;
        }
    };

    struct Matrix {
    protected:
        mutable bool valid = true;
    public:
        llong n_row;
        llong n_col;
        llong n_abs;
        smat Q;
        dmat R;

        Matrix(llong n_row, llong n_col, llong n_abs): n_row(n_row), n_col(n_col), n_abs(n_abs), Q(n_row, n_col), R(n_row, n_abs) {}

        Matrix(const Matrix& m): n_row(m.n_row), n_col(m.n_col), n_abs(m.n_abs), Q(m.Q), R(m.R) {
            m.valid = false;
        }

    };

    double psi_diploid(const llong i, const llong N, const double s = 0, const double h = 0.5, const double u = 1e-9, const double v = 1e-9);
    Row binom_row(const llong i, const llong Nx, const llong Ny, const double s = 0, const double h = 0.5, const double u = 1e-9, const double v = 1e-9, const double alpha = 1e-20);

    // Harrod matrix to solve for equilibrium distribution
    Matrix Equilibrium(const llong N, const double s = 0, const double h = 0.5, const double u = 1e-9, const double v = 1e-9, const double alpha = 1e-20, const llong block_size = 100);

    // Single - one matrix of a given absorption type
    Matrix Single(const llong Nx, const llong Ny, const absorption_type a_t, const double s = 0, const double h = 0.5, const double u = 1e-9, const double v = 1e-9, const double alpha = 1e-20, const llong block_size = 100);
    Matrix SingleAlt(const llong Nx, const llong Ny, const absorption_type a_t, const double s = 0, const double h = 0.5, const double u = 1e-9, const double v = 1e-9, const double alpha = 1e-20, const llong block_size = 100);

    // Switching - each sub-model is of the same absorbing type
    Matrix Switching(const lvec& N, const absorption_type a_t, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching, const double alpha = 1e-20, const llong block_size = 100);

    // Two-model Switching: A is of type `NON_ABSORBING`, B is of type `FIXATION_ONLY`
    Matrix NonAbsorbingToFixationOnly(const lvec& N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching, const double alpha = 1e-20, const llong block_size = 100);

    // Two-model Switching: A is of type `NON_ABSORBING`, B is of type `BOTH_ABSORBING`
    Matrix NonAbsorbingToBothAbsorbing(const lvec& N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching, const double alpha = 1e-20, const llong block_size = 100);

};
