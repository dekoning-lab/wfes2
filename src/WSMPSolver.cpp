#include "WSMPSolver.hpp"

WSMPSolver::WSMPSolver(SparseMatrix& A, llong n_rhs):
	size(A.n_row),
	n_right_hand_sides(n_rhs),
	m(A),
	iparm(lvec::Zero(64)),
	dparm(dvec::Zero(64)),
	workspace(dvec::Zero(size * n_rhs))
{
	iparm[0] = 0;  // default iparm
	iparm[4] = 0;  // C-style indexing
}

