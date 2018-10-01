#pragma once

#include "MKL_Const.hpp"
#include "SparseMatrix.hpp"

extern "C" {
	void wgsmp_(llong *n, llong *ia, llong *ja, double *avals,  double *b, llong *ldb, llong *nrhs, double *rmisc, llong *iparm, double *dparm);
	void wsetmaxthrds_(int* nthrds);
}

class WSMPSolver {

	llong size;
	llong n_right_hand_sides;
	llong phase;

	SparseMatrix& m;
	lvec iparm;
	dvec dparm;
	dvec workspace;

	public:
	WSMPSolver(SparseMatrix& A, llong nrhs = 1, int n_threads = 16);
	//~WSMPSolver();
	void analyze();
	dvec solve(const dvec& b, bool transpose = false);
};

