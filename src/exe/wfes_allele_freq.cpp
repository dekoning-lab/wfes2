#include "common.hpp"
#include "WrightFisher.hpp"
#include "PardisoSolver.hpp"
#include "parsing.hpp"
#include "args.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
using namespace std;

dvec equilibrium(llong N, double s, double h, double u, double v, double alpha, bool verbose = false) {
    WF::Matrix wf_eq = WF::Equilibrium(N, s, h, u, v, alpha, verbose);

    llong msg_level = verbose ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    PardisoSolver solver(wf_eq.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
    solver.analyze();

    dvec id = dvec::Zero(wf_eq.Q.n_row);
    id(wf_eq.Q.n_row - 1) = 1;

    dvec eq = solver.solve(id, true);
    eq = eq.array().abs();
    eq /= eq.sum();

    return eq.matrix();
}

void iterate_generations(dvec& x, llong N, llong t, double s, double h, double u, double v, double alpha, bool verbose = false) {
    WF::Matrix wf = WF::Single(N, N, WF::NON_ABSORBING, s, h, u, v, true, alpha, verbose);
    wf.Q.multiply_inplace_rep(x, t, true);
    x /= x.sum();
}

dvec switch_population_size(dvec& x, llong Nx, llong Ny, double s, double h, double u, double v, double alpha, bool verbose = false) {
    WF::Matrix wf = WF::Single(Nx, Ny, WF::NON_ABSORBING, s, h, u, v, true, alpha, verbose);
    dvec next = wf.Q.multiply(x, true);
    next /= next.sum();
    return next;
}

int main(int argc, char const *argv[])
{
    args::ArgumentParser parser("WFAFLE - Wright-Fisher Allele Frequency calculation");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;

    args::ValueFlag<Eigen::Matrix<llong, Eigen::Dynamic, 1>, NumericVectorReader<llong>> population_sizes_f(parser, "k", "Population sizes", {'N', "pop-sizes"}, args::Options::Required);
    args::ValueFlag<Eigen::Matrix<llong, Eigen::Dynamic, 1>, NumericVectorReader<llong>> epoch_lengths_f(parser, "k", "Number of generations for each epoch", {'G', "generations"}, args::Options::Required);
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> selection_coefficient_f(parser, "float[k]", "Selection coefficients", {'s', "selection"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> dominance_f(parser, "float[k]", "Dominance coefficients", {'h', "dominance"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> backward_mutation_f(parser, "float[k]", "Backward mutation rates", {'u', "backward-mu"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> forward_mutation_f(parser, "float[k]", "Forward mutation rates", {'v', "forward-mu"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::ValueFlag<string> initial_f(parser, "path", "Path to initial probability distribution CSV", {'i', "initial"});
    args::ValueFlag<llong>  initial_count_f(parser, "int", "Initial allele count", {'p', "initial-count"});
    args::ValueFlag<llong>  n_threads_f(parser, "int", "Number of threads", {'t', "num-threads"});
    args::Flag verbose_f(parser, "verbose", "Verbose solver output", {"verbose"});

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

    time_point t_start, t_end;
    if (verbose_f) t_start = std::chrono::system_clock::now();

    lvec pop_sizes(args::get(population_sizes_f));
    lvec epoch_gens(args::get(epoch_lengths_f));

    if(pop_sizes.size() != epoch_gens.size()) throw runtime_error("Population size vector should be the same length as epoch lengths");
    llong k = pop_sizes.size();

    dvec s = selection_coefficient_f ? args::get(selection_coefficient_f) : dvec::Constant(k, 0);
    dvec h = dominance_f ? args::get(dominance_f) : dvec::Constant(k, 0.5);
    dvec u = backward_mutation_f ? args::get(backward_mutation_f) : dvec::Constant(k, 1e-9);
    dvec v = forward_mutation_f ? args::get(forward_mutation_f) : dvec::Constant(k, 1e-9);
    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    llong n_threads = n_threads_f ? args::get(n_threads_f) : 1;

    #ifdef OMP
        omp_set_num_threads(n_threads);
    #endif
    mkl_set_num_threads(n_threads);

    deque<dvec> d;
    dvec initial;
    if (initial_f) {
        initial = load_csv_vector(args::get(initial_f));
    } else if (initial_count_f) {
        llong p = args::get(initial_count_f);
        initial = dvec::Zero(2 * pop_sizes(0) + 1);
        initial[p] = 1;
    } else {
        initial = equilibrium(pop_sizes(0), s(0), h(0), u(0), v(0), a, verbose_f);
    }
    d.push_back(initial);

    for(llong i = 0; i < k - 1; i++) {
        iterate_generations(d[i], pop_sizes(i), epoch_gens(i), s(i), h(i), u(i), v(i), a, verbose_f);
        d.push_back(switch_population_size(d[i], pop_sizes(i), pop_sizes(i + 1), s(i + 1), h(i + 1), u(i + 1), v(i + 1), a));
    }

    iterate_generations(d[k - 1], pop_sizes(k - 1), epoch_gens(k - 1), s(k - 1), h(k - 1), u(k - 1), v(k - 1), a, verbose_f);

    write_vector_to_file(d[k - 1], "stdout");

    if (verbose_f) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Total runtime: " << dt.count() << " s" << std::endl;
    }

    return EXIT_SUCCESS;
}
