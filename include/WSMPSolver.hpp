#pragma once

#include "MKL_Const.hpp"
#include "SparseMatrix.hpp"

extern void wgsmp(llong *n, llong *ia, llong *ja, double *avals,  double *b, double *ldb, int *nrhs, double *rmisc, llong *iparm, llong *dparm);
extern double wsmprtc();

class WSMPSolver {

	llong size;
	llong n_right_hand_sides;
	llong phase;

	SparseMatrix& m;
	lvec iparm;
	dvec dparm;
	dvec workspace;

	public:
	WSMPSolver(SparseMatrix& A, llong nrhs = 1);
	~WSMPSolver();
	void analyze();
	dvec solve(dvec& b, bool transpose = false);
};

