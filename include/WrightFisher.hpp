#pragma once

#include "PardisoSolver.hpp"
#include "SparseMatrix.hpp"
#include "common.hpp"
#include "rdist.hpp"

namespace WrightFisher {

enum absorption_type { NON_ABSORBING = 0, EXTINCTION_ONLY, FIXATION_ONLY, BOTH_ABSORBING };

inline llong n_absorbing(absorption_type a_t) {
    switch (a_t) {
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

inline std::string absorption_type_desc(absorption_type a_t) {
    switch (a_t) {
    case NON_ABSORBING:
        return "No absorbing boundaries";
    case EXTINCTION_ONLY:
        return "Only extinction boundary is absorbing";
    case FIXATION_ONLY:
        return "Only fixation boundary is absorbing";
    case BOTH_ABSORBING:
        return "Both extinction and fixation boundaries are absorbing";
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

    Row(llong start, llong end)
        : start(start), end(end), size(end - start + 1), weight(1), Q(end - start + 1) {}
    Row() : start(0), end(0), size(0), weight(0), Q((llong)0) {}

    Row(const Row &r) : start(r.start), end(r.end), size(r.size), weight(r.weight), Q(r.Q) {
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

    Matrix(llong n_row, llong n_col, llong n_abs)
        : n_row(n_row), n_col(n_col), n_abs(n_abs), Q(n_row, n_col), R(n_row, n_abs) {
        R.setZero();
    }

    Matrix(const Matrix &m) : n_row(m.n_row), n_col(m.n_col), n_abs(m.n_abs), Q(m.Q), R(m.R) {
        m.valid = false;
    }
};

double psi_diploid(const llong i, const llong N, const double s = 0, const double h = 0.5,
                   const double u = 1e-9, const double v = 1e-9);
Row binom_row(const llong size, const double p, const double alpha = 1e-20);

// Harrod matrix to solve for equilibrium distribution
Matrix EquilibriumSolvingMatrix(const llong N, const double s = 0, const double h = 0.5,
                                const double u = 1e-9, const double v = 1e-9,
                                const double alpha = 1e-20, const bool verbose = false,
                                const llong block_size = 100);
dmat Equilibrium(llong N, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9,
                 double alpha = 1e-20, bool verbose = false);

// Single - one matrix of a given absorption type
Matrix Single(const llong Nx, const llong Ny, const absorption_type a_t, const double s = 0,
              const double h = 0.5, const double u = 1e-9, const double v = 1e-9,
              const bool recurrent_mutation = true, const double alpha = 1e-20,
              const bool verbose = false, const llong block_size = 100);

// Bounce - mutation return to a non-zero count
Matrix Bounce(const llong Nx, const llong Ny, const double s = 0, const double h = 0.5,
              const double u = 1e-9, const double v = 1e-9, const bool recurrent_mutation = true,
              const double alpha = 1e-20, const bool verbose = false, const llong block_size = 100);

// Dual mutation - 0 copies only absorbing after the first mutation
Matrix DualMutation(const llong Nx, const llong Ny, const double s = 0, const double h = 0.5,
                    const double u = 1e-9, const double v = 1e-9,
                    const bool recurrent_mutation = true, const double alpha = 1e-20,
                    const bool verbose = false, const llong block_size = 100);

// Single but with entries larger than `t` summed into the fixation state
Matrix Truncated(const llong Nx, const llong Ny, const llong t, const double s, const double h,
                 const double u, const double v, bool recurrent_mutation = true,
                 const double alpha = 1e-20, const bool verbose = false,
                 const llong block_size = 100);

// Switching - each sub-model is of the same absorbing type
Matrix Switching(const lvec &N, const absorption_type a_t, const dvec &s, const dvec &h,
                 const dvec &u, const dvec &v, const dmat &switching, const double alpha = 1e-20,
                 const bool verbose = false, const llong block_size = 100);

// Two-model Switching: A is of type `NON_ABSORBING`, B is of type `FIXATION_ONLY`
Matrix NonAbsorbingToFixationOnly(const llong N, const dvec &s, const dvec &h, const dvec &u,
                                  const dvec &v, const dmat &switching, const double alpha = 1e-20,
                                  const bool verbose = false, const llong block_size = 100);

// Two-model Switching: A is of type `NON_ABSORBING`, B is of type `BOTH_ABSORBING`
Matrix NonAbsorbingToBothAbsorbing(const llong N, const dvec &s, const dvec &h, const dvec &u,
                                   const dvec &v, const dmat &switching, const double alpha = 1e-20,
                                   const bool verbose = false, const llong block_size = 100);

}; // namespace WrightFisher
