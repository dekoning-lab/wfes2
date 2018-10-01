#include "WSMPSolver.hpp"

WSMPSolver::WSMPSolver(SparseMatrix& A, llong n_rhs, int n_threads):
	size(A.n_row),
	n_right_hand_sides(n_rhs),
	m(A),
	iparm(lvec::Zero(64)),
	dparm(dvec::Zero(64)),
	workspace(dvec::Zero(size * n_rhs))
{
	wsetmaxthrds_(&n_threads);
	iparm[0] = 0;  // default iparm
	iparm[4] = 0;  // C-style indexing
}

void WSMPSolver::analyze()
{
	double* dummy = NULL;
	// --- Initialization {{{
	iparm[1] = 0;
	iparm[2] = 0;
	//wgsmp_(llong*, llong*,      llong*, double*, double*, double*, int*,              double*, llong*, llong*)
	wgsmp_(&m.n_row, m.row_index, m.cols, m.data, dummy, &m.n_col, &n_right_hand_sides, dummy, iparm.data(), dparm.data());

	if (iparm[63] != 0) {
		printf("The following ERROR was detected in S0: %lld\n", iparm[63]);
		exit(1);
	}
	// }}}

	// --- Analysis {{{
	iparm[1] = 1;
	iparm[2] = 1;
	wgsmp_(&m.n_row, m.row_index, m.cols, m.data, dummy, &m.n_col, &n_right_hand_sides, dummy, iparm.data(), dparm.data());

	if (iparm[63] != 0) {
		printf("The following ERROR was detected in S1: %lld\n", iparm[63]);
		exit(1);
	}
	// printf("NNZ 1000x%lld; FLOPS: %f\n", iparm[23], dparm[23]);
	// }}}

	// --- Factorization {{{
	iparm[1] = 2;
	iparm[2] = 2;
	wgsmp_(&m.n_row, m.row_index, m.cols, m.data, dummy, &m.n_col, &n_right_hand_sides, dummy, iparm.data(), dparm.data());

	if (iparm[63] != 0) {
		printf("The following ERROR was detected in S2: %lld\n", iparm[63]);
		exit(1);
	}
	// }}}

}

dvec WSMPSolver::solve(const dvec& x, bool transpose)
{
	double* dummy = NULL;
	dvec sol(x);
	if (transpose) {
		iparm[29] = 4;
	} else {
		iparm[29] = 0;
	}

	// --- Back substitution {{{
	iparm[1] = 3;
	iparm[2] = 3;
	wgsmp_(&m.n_row, m.row_index, m.cols, m.data, sol.data(), &m.n_col, &n_right_hand_sides, dummy, iparm.data(), dparm.data());

	if (iparm[63] != 0) {
		printf("The following ERROR was detected in S3: %lld\n", iparm[63]);
		exit(1);
	}
	// }}}

	// --- Iterative refinement {{{
	iparm[1] = 4;
	iparm[2] = 4;
	wgsmp_(&m.n_row, m.row_index, m.cols, m.data, sol.data(), &m.n_col, &n_right_hand_sides, dummy, iparm.data(), dparm.data());

	if (iparm[63] != 0) {
		printf("The following ERROR was detected in S4: %lld\n", iparm[63]);
		exit(1);
	}	
	return sol;

}

/*
WSMPSolver::~WSMPSolver()
{
	// Should there be cleanup?
}
*/
