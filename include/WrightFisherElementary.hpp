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
    std::pair<dmat, dmat> EquilibriumSolvingMatrix(const llong N, const double s = 0, const double h = 0.5, const double u = 1e-9, const double v = 1e-9);
    std::pair<dmat, dmat> NonAbsorbingToFixationOnly(const llong N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching);
    std::pair<dmat, dmat> NonAbsorbingToBothAbsorbing(const llong N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching);
    std::pair<dmat, dmat> SwitchingTwoByTwo(const lvec& N, const WrightFisher::absorption_type a_t, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching);
};
