#include "PardisoSolver.hpp"
#include "WrightFisher.hpp"
#include "args.hpp"
#include "common.hpp"
#include "parsing.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[]) {
    args::ArgumentParser parser("Time distribution of fixation with standing genetics variation");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;

    args::ValueFlag<llong> population_size_f(parser, "int", "Size of the population",
                                             {'N', "pop-size"}, args::Options::Required);
    args::ValueFlag<double> lambda_f(parser, "float", "Transition probability", {'l', "lambda"}, args::Options::Required);
    // Optional arguments
    args::ValueFlag<dvec, NumericVectorReader<double>> selection_coefficient_f(parser, "float[k]", "Selection coefficients", {'s', "selection"}, args::Options::Required);
    args::ValueFlag<dvec, NumericVectorReader<double>> dominance_f(parser, "float[k]", "Dominance coefficients", {'h', "dominance"});
    args::ValueFlag<dvec, NumericVectorReader<double>> backward_mutation_f(parser, "float[k]", "Backward mutation rates", {'u', "backward-mu"});
    args::ValueFlag<dvec, NumericVectorReader<double>> forward_mutation_f(parser, "float[k]", "Forward mutation rates", {'v', "forward-mu"});
    
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
    args::Flag force_f(parser, "force", "Do not perform parameter checks", {"force"});
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
    dvec selection_coefficient(args::get(selection_coefficient_f));
    dvec h = dominance_f ? args::get(dominance_f) : dvec::Constant(2, 0.5);
    dvec u = backward_mutation_f ? args::get(backward_mutation_f) : dvec::Constant(2, 1e-9);
    dvec v = forward_mutation_f ? args::get(forward_mutation_f) : dvec::Constant(2, 1e-9);
    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    double b = block_size_f ? args::get(block_size_f) : 100;
    double n_threads = n_threads_f ? args::get(n_threads_f) : 1;
    llong max_t = max_t_f ? args::get(max_t_f) : 100000;
    double integration_cutoff = integration_cutoff_f ? args::get(integration_cutoff_f) : 1 - 1e-8;
    bool no_rem = no_recurrent_mutation_f ? args::get(no_recurrent_mutation_f) : false;
    bool rem = !no_rem;
    
    if(selection_coefficient.size() != 2) throw runtime_error("Selection coefficient vector should be longer than 2");

    if (!force_f) {
        if (population_size > 500000) {
            throw args::Error("Population size is quite large - the computations will take a long time. Use --force to ignore");   
        }
        double max_mu = max(u.maxCoeff(), v.maxCoeff());
        if (4 * population_size * max_mu > 1) {
            throw args::Error("The mutation rate might violate the Wright-Fisher assumptions. Use --force to ignore");
        }
        if (selection_coefficient.minCoeff() <= -1) {
            throw args::Error("The selection coefficient is quite negative. Fixations might be impossible. Use --force to ignore");
        }
        if (a > 1e-5) {
            throw args::Error("Zero cutoff value is quite high. This might produce inaccurate results. Use --force to ignore");   
        }
    }
#ifdef OMP
    omp_set_num_threads(n_threads);
#endif
    mkl_set_num_threads(n_threads);

    time_point t_start, t_end;
    if (verbose_f)
        t_start = std::chrono::system_clock::now();

    double l = args::get(lambda_f);
    dmat switching(2, 2); switching << 1 - l, l, 0, 1;

    WF::Matrix wf = WF::NonAbsorbingToFixationOnly(population_size, selection_coefficient, h, u, v, switching, a, verbose_f, b);
    if (output_Q_f)
        wf.Q.save_market(args::get(output_Q_f));
    if (output_R_f)
        write_matrix_to_file(wf.R, args::get(output_R_f));

    dmat PH(max_t, 3);

    dvec c = dvec::Zero(4 * population_size + 1);
    c(0) = 1;
    dvec R = wf.R.col(0);

    double cdf = 0;
    llong i;
    for (i = 0; cdf < integration_cutoff && i < max_t; i++) {

	double P_abs_t = R.dot(c);
        cdf += P_abs_t;

        PH(i, 0) = i + 1;
        PH(i, 1) = P_abs_t;
        PH(i, 2) = cdf;

        c = wf.Q.multiply(c, true);
    }
    PH.conservativeResize(i, 3);

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
