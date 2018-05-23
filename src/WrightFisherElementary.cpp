#include "WrightFisherElementary.hpp"

namespace WF = WrightFisher;

double WrightFisherElementary::psi_diploid(long i, long N, double s, double h, double u, double v) { 
    long j = (2 * N) - i;
    double w_11 = 1 + s;
    double w_12 = 1 + (s * h);
    double w_22 = 1;
    double a = w_11 * i * i;
    double b = w_12 * i * j;
    double c = w_22 * j * j;
    double w_bar = a + (2 * b) + c;
    return (((a + b) * (1 - u)) + ((b + c) * v)) / w_bar;
}

dvec WrightFisherElementary::binom_row(const llong size, const double p) {
	dvec br(size + 1);
	for(llong i = 0; i <= size; i++) {
		br(i) = d_binom(i, size, p);
	}
	return br;
}

std::pair<dmat, dmat> WrightFisherElementary::Single(const llong Nx, const llong Ny, const WF::absorption_type a_t, const double s, const double h, const double u, const double v) {
    dmat Q;
    dmat R;
    
    switch(a_t) {
    case WF::NON_ABSORBING:
        Q = dmat((2 * Nx) + 1, (2 * Ny) + 1);
        R = dmat(0,0);
        for (long i = 0; i <= 2 * Nx; i++) {
            double p = WF::psi_diploid(i, Nx, s, h, u, v);
            for (long j = 0; j <= 2 * Ny; j++) {
            	Q(i, j) = d_binom(j, 2 * Ny, p);
            }
        }
        break;
    case WF::EXTINCTION_ONLY:
        Q = dmat(2 * Nx, 2 * Ny);
        R = dmat(2 * Nx, 1);
        for (long i = 1; i <= 2 * Nx; i++) {
            double p = WF::psi_diploid(i, Nx, s, h, u, v);
            R(i - 1, 0) = d_binom(0, 2 * Ny, p);
            for (long j = 1; j <= 2 * Ny; j++) {
                Q(i - 1, j - 1) = d_binom(j, 2 * Ny, p);
            }
        }
        break;
    case WF::FIXATION_ONLY:
        Q = dmat(2 * Nx, 2 * Ny);
        R = dmat(2 * Nx, 1);
        for (long i = 0; i <= (2 * Nx) - 1; i++) {
            double p = WF::psi_diploid(i, Nx, s, h, u, v);
            for (long j = 0; j <= (2 * Ny) - 1; j++) {
                Q(i, j) = d_binom(j, 2 * Ny, p);
            }
            R(i, 0) = d_binom(2 * Ny, 2 * Ny, p);
        }
        break;
    case WF::BOTH_ABSORBING:
        Q = dmat((2 * Nx) - 1, (2 * Ny) - 1);
        R = dmat((2 * Nx) - 1, 2);
        for (long i = 1; i <= (2 * Nx) - 1; i++) {
            double p = WF::psi_diploid(i, Nx, s, h, u, v);
            R(i - 1, 0) = d_binom(0, 2 * Ny, p);
            for (long j = 1; j <= (2 * Ny) - 1; j++) {
                Q(i - 1, j - 1) = d_binom(j, 2 * Ny, p);
            }
            R(i - 1, 1) = d_binom(2 * Ny, 2 * Ny, p);
        }
        break;
    }
    
   	return std::make_pair(Q, R); 
}