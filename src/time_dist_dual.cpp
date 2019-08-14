#include "PardisoSolver.hpp"
#include "WrightFisher.hpp"
#include "args.hpp"
#include "common.hpp"
#include "parsing.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[]) {
    args::ArgumentParser parser("Distribution of time to fixation / extinction");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;

    args::ValueFlag<llong> population_size_f(parser, "int", "Size of the population",
                                             {'N', "pop-size"}, args::Options::Required);

    // Optional arguments
    args::ValueFlag<double> selection_coefficient_f(parser, "float", "Selection coefficient",
                                                    {'s', "selection"});
    args::ValueFlag<double> dominance_f(parser, "float", "Dominance coefficient",
                                        {'h', "dominance"});
    args::ValueFlag<double> backward_mutation_f(parser, "float", "Backward mutation rate",
                                                {'u', "backward-mu"});
    args::ValueFlag<double> forward_mutation_f(parser, "float", "Forward mutation rate",
                                               {'v', "forward-mu"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::ValueFlag<llong> block_size_f(parser, "int", "Block size", {'b', "block-size"});
    args::ValueFlag<llong> n_threads_f(parser, "int", "Number of threads", {'t', "num-threads"});
    args::ValueFlag<double> integration_cutoff_f(
        parser, "float", "Stop once this probability mass is reached", {'c', "integration-cutoff"});
    args::ValueFlag<llong> max_t_f(parser, "int", "Maximum number of generations", {'m', "max-t"});
    args::Flag no_recurrent_mutation_f(parser, "bool", "Exclude recurrent mutation",
                                       {'r', "no-recurrent-mu"});

    args::ValueFlag<string> output_Q_f(parser, "path", "Output Q matrix to file", {"output-Q"});
    args::ValueFlag<string> output_R_f(parser, "path", "Output R vectors to file", {"output-R"});
    args::ValueFlag<string> output_P_f(parser, "path", "Output phase-type distribution",
                                       {"output-P"});

    args::Flag verbose_f(parser, "verbose", "Verbose solver output", {"verbose"});

    args::HelpFlag help_f(parser, "help", "Display this help menu", {"help"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help &) {
        cerr << parser;
        return EXIT_FAILURE;
    } catch (args::Error &e) {
        cerr << e.what() << endl;
        cerr << parser;
        return EXIT_FAILURE;
    }

    llong population_size = args::get(population_size_f);
    // Set default values
    double s = selection_coefficient_f ? args::get(selection_coefficient_f) : 0;
    double h = dominance_f ? args::get(dominance_f) : 0.5;
    double u = backward_mutation_f ? args::get(backward_mutation_f) : 1e-9;
    double v = forward_mutation_f ? args::get(forward_mutation_f) : 1e-9;
    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    double b = block_size_f ? args::get(block_size_f) : 100;
    double n_threads = n_threads_f ? args::get(n_threads_f) : 1;
    llong max_t = max_t_f ? args::get(max_t_f) : 100000;
    double integration_cutoff = integration_cutoff_f ? args::get(integration_cutoff_f) : 1 - 1e-8;
    bool no_rem = no_recurrent_mutation_f ? args::get(no_recurrent_mutation_f) : false;
    bool rem = !no_rem;

#ifdef OMP
    omp_set_num_threads(n_threads);
#endif
    mkl_set_num_threads(n_threads);

    time_point t_start, t_end;
    if (verbose_f)
        t_start = std::chrono::system_clock::now();

    WF::Matrix wf =
        WF::DualMutation(population_size, population_size, s, h, u, v, rem, a, verbose_f, b);
    if (output_Q_f)
        wf.Q.save_market(args::get(output_Q_f));
    if (output_R_f)
        write_matrix_to_file(wf.R, args::get(output_R_f));

    dmat PH(max_t, 5);

    dvec c = dvec::Zero(2 * population_size);
    c(0) = 1;

    double cdf = 0;
    llong i;
    for (i = 0; cdf < integration_cutoff && i < max_t; i++) {

        double P_ext_t = wf.R.col(0).dot(c);
        double P_fix_t = wf.R.col(1).dot(c);
        cdf += P_fix_t + P_ext_t;

        PH(i, 0) = i + 1;
        PH(i, 1) = P_ext_t;
        PH(i, 2) = P_fix_t;
        PH(i, 3) = P_ext_t + P_fix_t;
        PH(i, 4) = cdf;

        c = wf.Q.multiply(c, true);
    }
    PH.conservativeResize(i, 5);

    if (output_P_f) {
        write_matrix_to_file(PH, args::get(output_P_f));
    }

    if (verbose_f) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Total runtime: " << dt.count() << " s" << std::endl;
    }

    return EXIT_SUCCESS;
}
