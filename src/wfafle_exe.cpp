#include "common.hpp"
#include "WrightFisher.hpp"
#include "PardisoSolver.hpp"
#include "parsing.hpp"
#include "args.hpp"

namespace WF = WrightFisher;
using namespace std;

dvec equilibrium(llong N, double s, double h, double u, double v, double alpha, bool verbose = false) {
    WF::Matrix wf_eq = WF::Equilibrium(N, s, h, u, v, alpha);

    llong msg_level = verbose ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    PardisoSolver solver(wf_eq.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
    solver.analyze();

    dvec id = dvec::Zero(wf_eq.Q.n_row);
    id(wf_eq.Q.n_row - 1) = 1;

    dvec eq = solver.solve(id, true);
    // eq = eq.abs();

    return eq;
}

void iterate_generations(dvec& x, llong N, llong t, double s, double h, double u, double v, double alpha, bool verbose = false) {
    WF::Matrix wf = WF::Single(N, N, WF::NON_ABSORBING, s, h, u, v, alpha);

    for(llong i = 0; i < t; i ++) {
        wf.Q.multiply_inplace(x, true);
    }
}

dvec switch_population_size(dvec& x, llong Nx, llong Ny, double s, double h, double u, double v, double alpha, bool verbose = false) {
    WF::Matrix wf = WF::Single(Nx, Ny, WF::NON_ABSORBING, s, h, u, v, alpha);
    dvec next = wf.Q.multiply(x, true);
    return next;
}

int main(int argc, char const *argv[])
{
    args::ArgumentParser parser("WFAFLE");
    args::ValueFlag<Eigen::Matrix<llong, Eigen::Dynamic, 1>, NumericVectorReader<llong>> population_sizes_f(parser, "k", "Population sizes", {'N', "pop-sizes"}, args::Options::Required);
    args::ValueFlag<Eigen::Matrix<llong, Eigen::Dynamic, 1>, NumericVectorReader<llong>> epoch_lengths_f(parser, "k", "Epoch lengths", {'E', "epochs"}, args::Options::Required);

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        cerr << parser;
        return EXIT_FAILURE;
    } catch (args::Error& e) {
        cerr << e.what() << endl;
        cerr << parser;
        return EXIT_FAILURE;
    }

    lvec pop_sizes(args::get(population_sizes_f));
    lvec epoch_gens(args::get(epoch_lengths_f));

    if(pop_sizes.size() != epoch_gens.size()) throw runtime_error("Population size vector should be the same length as epoch lengths");
    if(pop_sizes.size() < 2) throw runtime_error("Population size vector should be longer than 2");

    double s = 0;
    double h = 0.5;
    double mu = 1e-5;
    double alpha = 1e-20;

    llong k = pop_sizes.size();

    deque<dvec> d;
    d.push_back(equilibrium(pop_sizes(0), s, h, mu, mu, alpha));

    for(llong i = 0; i < k - 1; i++) {
        iterate_generations(d[i], pop_sizes(i), epoch_gens(i), s, h, mu, mu, alpha);
        d.push_back(switch_population_size(d[i], pop_sizes(i), pop_sizes(i + 1), s, h, mu, mu, alpha));
    }

    iterate_generations(d[k - 1], pop_sizes(k - 1), epoch_gens(k - 1), s, h, mu, mu, alpha);

    cout << d[k - 1] << endl;

    return EXIT_SUCCESS;
}
