#pragma once

#include "common.hpp"
#include "rdist.hpp"
#include "SparseMatrix.hpp"
#include "WrightFisher.hpp"


// This implements simple non-optimized versions of WFES functions. These should only be used in testing.
// Still using WF containers

namespace WrightFisherElementary {
	double psi_diploid(long i, long N, double s = 0, double h = 0.5, double u = 1e-9, double v = 1e-9);
    dvec binom_row(const llong size, const double p);
    std::pair<dmat, dmat> Single(const llong Nx, const llong Ny, const WrightFisher::absorption_type a_t, const double s = 0, const double h = 0.5, const double u = 1e-9, const double v = 1e-9);
};