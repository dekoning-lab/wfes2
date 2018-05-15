#include "common.hpp"
#include "WrightFisher.hpp"
#include "PardisoSolver.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[])
{
    int Nx = atoi(argv[1]);
    int Ny = atoi(argv[2]);

    // WF::Matrix wf_ref = WF::Single(Nx, Ny, WF::NEITHER, 0, 0.5, 1e-9, 1e-9, 0);
    // cout << wf_ref.Q.n_row << "x" << wf_ref.Q.n_col << endl;
    // wf_ref.Q.save_market("ref.mtx");
    // cout << wf_ref.Q.dense() << endl << endl;

    // dvec c = dvec::identity(wf_ref.Q.n_row);
    // dvec x = wf_ref.Q.multiply(c, true);

    // cout << x << endl;
    // cout << x.sum() << endl;

    WF::Matrix wf_eq = WF::Equilibrium(Nx, 0, 0.5, 1e-3, 1e-3, 1e-20);

    // llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;
    llong msg_level = MKL_PARDISO_MSG_VERBOSE;

    PardisoSolver solver(wf_eq.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
    solver.analyze();

    dvec id = dvec::identity(wf_eq.Q.n_row, wf_eq.Q.n_row - 1);

    dvec eq = solver.solve(id, true);

    cout << eq << endl;


    // WF::Matrix wf_ext = WF::Single(Nx, Ny, WF::EXTINCTION, 0, 0.5, 1e-9, 1e-9, 0);
    // cout << wf_ext.Q.n_row << "x" << wf_ext.Q.n_col << endl;
    // wf_ext.Q.save_market("ext.mtx");
    // cout << wf_ext.Q.dense() << endl << endl;
    // cout << wf_ext.R << endl << endl;

    // WF::Matrix wf_fix = WF::Single(Nx, Ny, WF::FIXATION, 0, 0.5, 1e-9, 1e-9, 0);
    // cout << wf_fix.Q.n_row << "x" << wf_fix.Q.n_col << endl;
    // wf_fix.Q.save_market("fix.mtx");
    // cout << wf_fix.Q.dense() << endl << endl;
    // cout << wf_fix.R << endl << endl;

    // WF::Matrix wf_abs = WF::Single(Nx, Ny, WF::BOTH, 0, 0.5, 1e-9, 1e-9, 0);
    // cout << wf_abs.Q.n_row << "x" << wf_abs.Q.n_col << endl;
    // wf_abs.Q.save_market("abs.mtx");
    // cout << wf_abs.Q.dense() << endl << endl;
    // cout << wf_abs.R << endl << endl;

    return 0;
}