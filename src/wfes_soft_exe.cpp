#include "common.hpp"
#include "WrightFisher.hpp"
#include "PardisoSolver.hpp"
#include "parsing.hpp"
#include "args.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[])
{
    args::ArgumentParser parser("WFES-SOFT");
    args::ValueFlag<llong> population_size_f(parser, "int", "Size of the population", {'N', "pop-size"}, args::Options::Required);
    args::ValueFlag<double> switching_time_f(parser, "double", "Expected number of generations in non-absorbing model", {'l', "lambda"}, args::Options::Required);
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> backward_mutation_f(parser, "float[k]", "Backward mutation rates", {'u', "backward-mu"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> forward_mutation_f(parser, "float[k]", "Forward mutation rates", {'v', "forward-mu"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> selection_coefficient_f(parser, "k", "Selection coefficients", {'s', "selection"}, args::Options::Required);
    args::ValueFlag<string> output_Q_f(parser, "path", "Output Q matrix to file", {"output-Q"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
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

    dvec selection_coefficient(args::get(selection_coefficient_f));
    llong population_size = args::get(population_size_f);
    dvec u = backward_mutation_f ? args::get(backward_mutation_f) : dvec::Constant(2, 1e-9);
    dvec v = forward_mutation_f ? args::get(forward_mutation_f) : dvec::Constant(2, 1e-9);
    double lambda = 1 / args::get(switching_time_f);
    double a = alpha_f ? args::get(alpha_f) : 1e-20;

    if(selection_coefficient.size() < 2) throw runtime_error("Population size vector should be longer than 2");

    dmat switching(2, 2);
    switching << 1 - lambda, lambda, 0, 1;

    dvec dom(2); dom << 0.5, 0.5;
    WF::Matrix wf = WF::NonAbsorbingToFixationOnly(population_size, selection_coefficient, dom, u, v, switching, a);
    if(output_Q_f) wf.Q.save_market(args::get(output_Q_f));


    wf.Q.subtract_identity();

    llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    PardisoSolver solver(wf.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
    solver.analyze();
    dvec id = dvec::Zero(wf.Q.n_row);
    id(0) = 1;

    dvec N = solver.solve(id, true);

    cout << N << endl;

    return 0;
}
