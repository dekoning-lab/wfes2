#include "common.hpp"
#include "WrightFisher.hpp"
#include "PardisoSolver.hpp"
#include "parsing.hpp"
#include "args.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[])
{
    args::ArgumentParser parser("Phase-type distribution - calculate moments of absorption times for a fixation-only model");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;
    
    args::ValueFlag<llong> population_size_f(parser, "int", "Size of the population", {'N', "pop-size"}, args::Options::Required);

    // Optional arguments
    args::ValueFlag<double> selection_coefficient_f(parser, "float", "Selection coefficient", {'s', "selection"});
    args::ValueFlag<double> dominance_f(parser, "float", "Dominance coefficient", {'h', "dominance"});
    args::ValueFlag<double> backward_mutation_f(parser, "float", "Backward mutation rate", {'u', "backward-mu"});
    args::ValueFlag<double> forward_mutation_f(parser, "float", "Forward mutation rate", {'v', "forward-mu"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::ValueFlag<llong>  block_size_f(parser, "int", "Block size", {'b', "block-size"});
    args::ValueFlag<llong>  n_threads_f(parser, "int", "Number of threads", {'t', "num-threads"});
    
    args::ValueFlag<llong> n_moments_f(parser, "int", "Numbers of moments to calculate", {'k', "n-moments"});

    args::Flag verbose_f(parser, "verbose", "Verbose solver output", {"verbose"});
    
    args::HelpFlag help_f(parser, "help", "Display this help menu", {"help"});
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

    llong population_size = args::get(population_size_f);
    // Set default values
    double s = selection_coefficient_f ? args::get(selection_coefficient_f) : 0;
    double h = dominance_f ? args::get(dominance_f) : 0.5;
    double u = backward_mutation_f ? args::get(backward_mutation_f) : 1e-9;
    double v = forward_mutation_f ? args::get(forward_mutation_f) : 1e-9;
    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    double b = block_size_f ? args::get(block_size_f) : 100;
    double n_threads = n_threads_f ? args::get(n_threads_f) : 1;
    llong k = n_moments_f ? args::get(n_moments_f) : 20;
    llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    llong size = 2 * population_size;

    #ifdef OMP
        omp_set_num_threads(n_threads);
    #endif
    mkl_set_num_threads(n_threads);

    time_point t_start, t_end;
    if (verbose_f) t_start = std::chrono::system_clock::now();

    WF::Matrix wf = WF::Single(population_size, population_size, WF::FIXATION_ONLY, s, h, u, v, true, a, verbose_f, b);
    wf.Q.subtract_identity();
    
    PardisoSolver solver(wf.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
    solver.analyze();
    
    dvec z = dvec::Zero(k+1);
    z(0) = 1; z(1) = -1;

    dvec rhs = dvec::Ones(size);
    dmat m = dmat::Zero(size, k+1);
    m.col(0) = rhs;
    m.col(1) = solver.solve(rhs, false);

    for(llong i = 1; i < k; i++) {
	z(i + 1) = -1;
	for (llong j = i; j > 0; j--) {
	    z(j) = z(j - 1) - z(j);
	}
	z(0) = -z(0);
	rhs.setZero();
	for (llong j = 0; j < i+1; j++) {
	    rhs += z(j) * m.col(j);
	}

	// note that we only need the first row of M - do we need to solve every time?
	m.col(i+1) = solver.solve(rhs, false);

	
        
    }
    double m1 = m(0, 1);
    double m2 = m(0, 2);
    std::cout << "Mean: " << m1 << std::endl;
    std::cout << "Standard deviation: " << sqrt(m2 - (m1 * m1)) << std::endl;
    std::cout << "Raw moments: " << std::endl;
    for (llong i = 1; i <= k; i++) {
	std::cout << i << "\t" << m(0, i) << std::endl;
    }


    if (verbose_f) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Total runtime: " << dt.count() << " s" << std::endl;
    }

    return EXIT_SUCCESS;
}
