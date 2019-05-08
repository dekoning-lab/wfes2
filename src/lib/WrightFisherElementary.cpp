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

std::pair<dmat, dmat> WrightFisherElementary::EquilibriumSolvingMatrix(const llong N, const double s, const double h, const double u, const double v) {
	std::pair<dmat, dmat> W = WrightFisherElementary::Single(N, N, WF::NON_ABSORBING, s, h, u, v);
	llong N2 = 2 * N;
	W.first = dmat::Identity(N2 + 1, N2 + 1) - W.first;
	W.first.col(2 * N) = dvec::Ones(N2 + 1);
	return W;
}

std::pair<dmat, dmat> WrightFisherElementary::NonAbsorbingToFixationOnly(const llong N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching) {


    lvec sizes(2); sizes << (2 * N) + 1, 2 * N;
    llong size = sizes.sum();

    dmat Q = dmat::Constant(size, size, -1);
    dmat R(size, 1);

    std::pair<dmat, dmat> W00 = WrightFisherElementary::Single(N, N, WF::NON_ABSORBING, s(0), h(0), u(0), v(0));
    W00.first *= switching(0, 0); W00.second *= switching(0, 0);

    std::pair<dmat, dmat> W01 = WrightFisherElementary::Single(N, N, WF::FIXATION_ONLY, s(1), h(1), u(1), v(1));
    W01.first *= switching(0, 1); W01.second *= switching(0, 1);

    dmat last_row(1, (2*N)+1);
    last_row.row(0) = WrightFisherElementary::binom_row(2 * N, psi_diploid(2 * N, N, s(1), h(1), u(1), v(1))) * switching(0, 1);

    std::pair<dmat, dmat> W10 = WrightFisherElementary::Single(N, N, WF::NON_ABSORBING, s(1), h(1), u(1), v(1));
    W10.first *= switching(1, 0); W10.second *= switching(1, 0);

    std::pair<dmat, dmat> W11 = WrightFisherElementary::Single(N, N, WF::FIXATION_ONLY, s(1), h(1), u(1), v(1));
    W11.first *= switching(1, 1); W10.second *= switching(1, 1);

    // build matrix
    Q.block(0,       0,       (2*N)+1, (2*N)+1) = W00.first;
    Q.block(0,       (2*N)+1, (2*N),   (2*N))   = W01.first;
    Q.block((2*N),   (2*N)+1, 1,       (2*N))   = last_row.row(0).head(2*N);
    Q.block((2*N)+1, 0,       (2*N),   (2*N)+1) = W10.first.block(0,0,2*N,(2*N)+1);
    Q.block((2*N)+1, (2*N)+1, (2*N),   (2*N))   = W11.first;

    // R << W00.second, W01.second, W10.second, W11.second;

    return std::make_pair(Q, R);
}

std::pair<dmat, dmat> WrightFisherElementary::NonAbsorbingToBothAbsorbing(const llong N, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching) {


    lvec sizes(2); sizes << (2 * N) + 1, (2 * N) - 1;
    llong size = sizes.sum();

    dmat Q = dmat::Constant(size, size, -1);
    dmat R(size, 2);

    std::pair<dmat, dmat> W00 = WrightFisherElementary::Single(N, N, WF::NON_ABSORBING, s(0), h(0), u(0), v(0));
    W00.first *= switching(0, 0); W00.second *= switching(0, 0);
    
    std::pair<dmat, dmat> W01 = WrightFisherElementary::Single(N, N, WF::BOTH_ABSORBING, s(1), h(1), u(1), v(1));
    W01.first *= switching(0, 1); W01.second *= switching(0, 1);
    
    dmat first_row(1, (2*N)+1);
    first_row.row(0) = WrightFisherElementary::binom_row(2 * N, psi_diploid(0, N, s(1), h(1), u(1), v(1))) * switching(0, 1);
    
    dmat last_row(1, (2*N)+1);
    last_row.row(0) = WrightFisherElementary::binom_row(2 * N, psi_diploid(2 * N, N, s(1), h(1), u(1), v(1))) * switching(0, 1);


    std::pair<dmat, dmat> W10 = WrightFisherElementary::Single(N, N, WF::NON_ABSORBING, s(1), h(1), u(1), v(1));
    W10.first *= switching(1, 0); W10.second *= switching(1, 0);
    
    std::pair<dmat, dmat> W11 = WrightFisherElementary::Single(N, N, WF::BOTH_ABSORBING, s(1), h(1), u(1), v(1));
    W11.first *= switching(1, 1); W10.second *= switching(1, 1);

    // build matrix
    Q.block(0,       0,       (2*N)+1, (2*N)+1) = W00.first;
    Q.block(0,       (2*N)+1, 1,       (2*N)-1) = first_row.row(0).segment(1, (2*N)-1);
    Q.block(1,       (2*N)+1, (2*N)-1, (2*N)-1) = W01.first;
    Q.block((2*N),   (2*N)+1, 1,       (2*N)-1) = last_row.row(0).segment(1, (2*N)-1);
    Q.block((2*N)+1, 0,       (2*N)-1, (2*N)+1) = W10.first.block(0,0,(2*N)-1,(2*N)+1);
    Q.block((2*N)+1, (2*N)+1, (2*N)-1, (2*N)-1) = W11.first;

    return std::make_pair(Q, R);
}

std::pair<dmat, dmat> WrightFisherElementary::SwitchingTwoByTwo(const lvec& N, const WF::absorption_type a_t, const dvec& s, const dvec& h, const dvec& u, const dvec& v, const dmat& switching) {
	llong k = N.size();
    llong n_abs_total = WF::n_absorbing(a_t) * k;
    lvec sizes(k);
    for(llong i = 0; i < k; i++) sizes(i) = 2 * N(i) + 1;
    llong size = sizes.sum();
	dmat Q(size - n_abs_total, size - n_abs_total);
	dmat R(size - n_abs_total, n_abs_total);

	std::deque<std::pair<dmat,dmat>> W;
	for(llong i = 0; i < k; i++) {
		for(llong j = 0; j < k; j++) {
			std::pair<dmat, dmat> w = WrightFisherElementary::Single(N(i), N(j), a_t, s(j), h(j), u(j), v(j));
			w.first *= switching(i, j);
			w.second *= switching(i, j);
			W.push_back(w);
		}
	}

	// My, my, is this ugly
	Q << W[0].first, W[1].first, W[2].first, W[3].first;
	R << W[0].second, W[1].second, W[2].second, W[3].second;
	return std::make_pair(Q, R);
}
