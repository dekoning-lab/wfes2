
#include "common.hpp"
#include "WrightFisher.hpp"
#include "PardisoSolver.hpp"
#include "parsing.hpp"
#include "util.hpp"
#include "args.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[])
{

    // Arguments {{{

    // Parser setting{{{
    args::ArgumentParser parser("WFES-SWITCHING");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;
    // }}}

    // Required args {{{
    args::ValueFlag<lvec, NumericVectorReader<llong>> population_size_f(parser, "int[k]", "Sizes of the populations", {'N', "pop-sizes"}, args::Options::Required);
    args::ValueFlag<dvec, NumericVectorReader<double>> generations_f(parser, "float[k]", "Expected number of generations spent in each model", {'G', "generations"}, args::Options::Required);
    args::ValueFlag<double> factor_f(parser, "float", "Matrix approximation factor", {'f', "factor"}, args::Options::Required);
    // }}}

    // Optional arguments {{{
    args::ValueFlag<double> selection_coefficient_f(parser, "float", "2Ns", {'s', "selection"});
    args::ValueFlag<double> dominance_f(parser, "float", "Dominance coefficients", {'h', "dominance"});
    args::ValueFlag<double> backward_mutation_f(parser, "float", "4Nu", {'u', "backward-mu"});
    args::ValueFlag<double> forward_mutation_f(parser, "float", "4Nv", {'v', "forward-mu"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::ValueFlag<llong>  n_threads_f(parser, "int", "Number of threads", {'t', "num-threads"});

    args::ValueFlag<string> initial_f(parser, "path", "Path to initial probability distribution CSV", {'i', "initial"});
    args::ValueFlag<llong>  initial_count_f(parser, "int", "Initial allele count", {'p', "initial-count"});
    // }}}

    // Output options {{{
    args::ValueFlag<string> output_Q_f(parser, "path", "Output Q matrix to file", {"output-Q"});
    args::ValueFlag<string> output_R_f(parser, "path", "Output R vectors to file", {"output-R"});
    args::ValueFlag<string> output_N_f(parser, "path", "Output N matrix to file", {"output-N"});
    args::ValueFlag<string> output_B_f(parser, "path", "Output B vectors to file", {"output-B"});
    args::ValueFlag<string> output_N_ext_f(parser, "path", "Output extinction-conditional sojourn to file", {"output-N-ext"});
    args::ValueFlag<string> output_N_fix_f(parser, "path", "Output fixation-conditional sojourn to file", {"output-N-fix"});
    args::ValueFlag<string> output_N_tmo_f(parser, "path", "Output timeout-conditional sojourn to file", {"output-N-tmo"});
    // }}}

    // Flags{{{
    args::Flag csv_f(parser, "csv", "Output results in CSV format", {"csv"});
    args::Flag force_f(parser, "force", "Do not perform parameter checks", {"force"});
    args::Flag verbose_f(parser, "verbose", "Verbose solver output", {"verbose"});
    args::HelpFlag help_f(parser, "help", "Display this help menu", {"help"});
    // }}}

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
    // }}}

    // Start timer{{{
    time_point t_start, t_end;
    if (verbose_f) t_start = std::chrono::system_clock::now();
    // }}}

    // Set defaults {{{
    lvec population_sizes = args::get(population_size_f);
    llong n_models = population_sizes.size();
    dvec t = args::get(generations_f);
    double f = args::get(factor_f);

    // Set default values
    double gamma = selection_coefficient_f ? args::get(selection_coefficient_f) : 0;
    double dom = dominance_f ? args::get(dominance_f) : 0.5;
    double theta_u = backward_mutation_f ? args::get(backward_mutation_f) : 1e-6;
    double theta_v = forward_mutation_f ? args::get(forward_mutation_f) : 1e-6;

    population_sizes /= f;
    t /= f;

    // scale by population size
    dvec M = population_sizes.cast<double>();
    dvec s = gamma / (2 * M.array());
    dvec h = dvec::Constant(n_models, dom);
    dvec u = theta_u / (4 * M.array());
    dvec v = theta_v / (4 * M.array());

    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    llong n_threads = n_threads_f ? args::get(n_threads_f) : 1;

    std::cerr << population_sizes.transpose() << std::endl;
    std::cerr << t.transpose() << std::endl;
    std::cerr << s.transpose() << std::endl;
    std::cerr << u.transpose() << std::endl;
    std::cerr << v.transpose() << std::endl;


#ifdef OMP
    omp_set_num_threads(n_threads);
#endif
    mkl_set_num_threads(n_threads);

    llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    //}}}

    // Main calculation {{{ 
    llong size = (2 * population_sizes.sum()) + n_models;

    dmat switching = dmat::Zero(n_models, n_models);
    for(llong i = 0; i < n_models - 1; i++) {
        switching(i, i) = 1 - (1 / t(i));
        switching(i, i+1) = 1 / t(i);
    }
    switching(n_models - 1, n_models - 1) = 1 - (1 / t(n_models - 1));

    WF::Matrix W = WF::Switching(population_sizes, WF::NON_ABSORBING, 
            s, h, u, v, switching, a, verbose_f);


    // if(output_Q_f) W.Q.save_market(args::get(output_Q_f));

    W.Q.subtract_identity();

    dvec initial;
    if (initial_f) {
        initial = load_csv_vector(args::get(initial_f));
    } else if (initial_count_f) {
        llong p = args::get(initial_count_f);
        initial = dvec::Zero(2 * population_sizes(0) + 1);
        initial[p] = 1;
    } else {
        initial = WF::Equilibrium(population_sizes(0), s(0), h(0), u(0), v(0), a, verbose_f);
    }

    // llong n_rhs = 2 * population_sizes(0) + 1;
    llong n_rhs = 2 * population_sizes(0) + 1;
    // llong n_rhs = 2 * population_sizes(n_models - 1) + 1;

    PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level, n_rhs);
    solver.analyze();

    // dmat R = dmat::Identity(size, 2 * population_sizes(n_models - 1) + 1).reverse() * (1 / t(n_models - 1));
    llong nk = 2 * population_sizes(n_models - 1) + 1;
    // SparseMatrix R = SparseMatrix::LeftPaddedDiagonal(nk, 1 / t(n_models - 1), size - nk);
    // if(output_R_f) R.save_market(args::get(output_R_f));


    dmat id = dmat::Identity(n_rhs, size);

    dmat B = solver.solve_multiple(id, true);

    // if(output_N_f) write_matrix_to_file(Nt, args::get(output_N_f));

    B /= t(n_models - 1);

    // std::cout << eq.rows() << " x " << eq.cols() << std::endl;
    // std::cout << B.transpose().rightCols(nk).rows() << " x " << B.transpose().rightCols(nk).cols() << std::endl;
    dvec d = initial.transpose() * B.transpose().rightCols(nk);

    // dvec d = B.transpose().rightCols(nk).row(0);



    if (f != 0) {
        llong lt = n_models - 1;
        WF::Matrix sw_up = WF::Single(population_sizes(lt), population_sizes(lt) * f, WF::NON_ABSORBING, s(lt), h(lt), u(lt), v(lt), true, a, verbose_f);
        // WF::Matrix sw_down = WF::Single(population_sizes(lt) * f, population_sizes(lt), WF::NON_ABSORBING, s(0), h(0), u(0), v(0), true, a, verbose_f);

        dvec e = sw_up.Q.multiply(d, true);
        // dvec f = sw_down.Q.multiply(e, true);
        d = e;
    }

    std::cout << d << std::endl;
    // write_vector_to_file(d, "stdout");
    // if(output_B_f) write_matrix_to_file(Nt.transpose().rightCols(nk), args::get(output_B_f));

    // dmat Nt = dmat::Zero(size, n_rhs);
    // for(llong i = 0; i < n_rhs; i++) {
    //     dvec id_tmp = id.col(i);
    //     Nt.col(i) = solver.solve(id_tmp, true);
    // }

    // SparseMatrix B = R.multiply(Ntt);


    // if(output_N_f) write_matrix_to_file(Nt, args::get(output_N_f));
    // if(output_B_f) B.save_market(args::get(output_B_f));
    // if(output_B_f) write_matrix_to_file(B.transpose(), args::get(output_B_f));





// }}}

    // Print timing {{{ 
    if (verbose_f) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Total runtime: " << dt.count() << " s" << std::endl;
    }
    // }}}

    return EXIT_SUCCESS;

}
