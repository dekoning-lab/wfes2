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
    args::ArgumentParser parser("WFES-SWITCHING");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;

    args::Group model_f(parser, "Model type - specify one", args::Group::Validators::Xor, args::Options::Required);
    args::Flag absorption_f(model_f, "absorption", "Both fixation and extinction states are absorbing", {"absorption"});
    args::Flag fixation_f(model_f, "fixation", "Only fixation state is absorbing", {"fixation"});
    // args::Flag fundamental_f(model_f, "fundamental", "Calculate the entire fundamental matrix (slow)", {"fundamental"});
    // args::Flag equilibrium_f(model_f, "equilibrium", "Calculate the equilibrium distribtion of allele states", {"equilibrium"});

    args::ValueFlag<Eigen::Matrix<llong, Eigen::Dynamic, 1>, NumericVectorReader<llong>> population_size_f(parser, "int[k]", "Sizes of the populations", {'N', "pop-sizes"}, args::Options::Required);

    // Optional arguments
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> selection_coefficient_f(parser, "float[k]", "Selection coefficients", {'s', "selection"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> dominance_f(parser, "float[k]", "Dominance coefficients", {'h', "dominance"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> backward_mutation_f(parser, "float[k]", "Backward mutation rates", {'u', "backward-mu"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> forward_mutation_f(parser, "float[k]", "Forward mutation rates", {'v', "forward-mu"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> starting_prob_f(parser, "float[k]", "Starting probabilities", {'p', "starting-prob"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, NumericMatrixReader<double>> switching_f(parser, "float[k][k]", "Switching parameters over models", {'r', "switching"});
    args::ValueFlag<double> integration_cutoff_f(parser, "float", "Starting number of copies integration cutoff", {'c', "integration-cutoff"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::ValueFlag<llong>  n_threads_f(parser, "int", "Number of threads", {'t', "num-threads"});

    // TODO: add integration cutoff

    // Output options
    args::ValueFlag<string> output_Q_f(parser, "path", "Output Q matrix to file", {"output-Q"});
    args::ValueFlag<string> output_R_f(parser, "path", "Output R vectors to file", {"output-R"});
    args::ValueFlag<string> output_N_f(parser, "path", "Output N matrix to file", {"output-N"});
    args::ValueFlag<string> output_B_f(parser, "path", "Output B vectors to file", {"output-B"});
    args::ValueFlag<string> output_N_ext_f(parser, "path", "Output extinction-conditional sojourn to file", {"output-N-ext"});
    args::ValueFlag<string> output_N_fix_f(parser, "path", "Output fixation-conditional sojourn to file", {"output-N-fix"});
    // args::ValueFlag<string> output_E_f(parser, "path", "Output Equilibrium frequencies to file (--equilibrium only)", {"output-E"});

    args::Flag csv_f(parser, "csv", "Output results in CSV format", {"csv"});
    args::Flag force_f(parser, "force", "Do not perform parameter checks", {"force"});
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

    time_point t_start, t_end;
    if (verbose_f) t_start = std::chrono::system_clock::now();

    if(model_f.MatchedChildren() != 1) {
        throw args::Error("Should have exactly one of the 'Model type' options");
    }
    // }}}

    // Set defaults {{{
    lvec population_sizes = args::get(population_size_f);
    llong n_models = population_sizes.size();

    // Set default values
    dvec s = selection_coefficient_f ? args::get(selection_coefficient_f) : dvec::Constant(n_models, 0);
    dvec h = dominance_f ? args::get(dominance_f) : dvec::Constant(n_models, 0.5);
    dvec u = backward_mutation_f ? args::get(backward_mutation_f) : dvec::Constant(n_models, 1e-9);
    dvec v = forward_mutation_f ? args::get(forward_mutation_f) : dvec::Constant(n_models, 1e-9);
    dvec p = starting_prob_f ? args::get(starting_prob_f) : dvec::Constant(n_models, 1.0 / (double)(n_models));
    dmat switching = switching_f ? args::get(switching_f) : dmat::Ones(n_models, n_models);

    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    double integration_cutoff = integration_cutoff_f ? args::get(integration_cutoff_f) : 1e-10;
    llong n_threads = n_threads_f ? args::get(n_threads_f) : 1;

#ifdef OMP
    omp_set_num_threads(n_threads);
#endif
    mkl_set_num_threads(n_threads);

    llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    if (!force_f) {
        if (population_sizes.maxCoeff() > 500000) {
            throw args::Error("Population size is quite large - the computations will take a long time. Use --force to ignore");   
        }
        dvec N = population_sizes.cast<double>();
        dvec theta_f = dvec::Constant(n_models, 4).array() * N.array() * v.array();
        dvec theta_b = dvec::Constant(n_models, 4).array() * N.array() * u.array();
        double max_theta = max(theta_b.maxCoeff(), theta_f.maxCoeff());
        if (max_theta > 1) {
            throw args::Error("The mutation rate might violate the Wright-Fisher assumptions. Use --force to ignore");
        }
        dvec gamma = dvec::Constant(n_models, 2).array() * N.array() * s.array();
        if (2 * N.maxCoeff() * s.minCoeff() <= -100) {
            throw args::Error("The selection coefficient is quite negative. Fixations might be impossible. Use --force to ignore");
        }
        if (a > 1e-5) {
            throw args::Error("Zero cutoff value is quite high. This might produce inaccurate results. Use --force to ignore");   
        }
    }

    // dmat switching = GTR::Matrix(p, r);

    // switching.diagonal() = dvec::Zero(n_models);
    // switching /= switching.sum();
    dvec row_sums = switching.rowwise().sum();
    for (llong i = 0; i < n_models; i++) {
        for (llong j = 0; j < n_models; j++) {
            switching(i, j) /= row_sums(i);
        }
    }
    //}}}

    if(fixation_f) // BEGIN SWITCHING FIXATION {{{
    {
        WF::Matrix W = WF::Switching(population_sizes, WF::FIXATION_ONLY, s, h, u, v, switching, a, verbose_f);

        if(output_Q_f) W.Q.save_market(args::get(output_Q_f));
        if(output_R_f) write_matrix_to_file(W.R, args::get(output_R_f));

        W.Q.subtract_identity();

        llong size = (2 * population_sizes.sum());
        dvec id(size);

        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();

        dmat N(n_models, size);
        lvec start_state_index(n_models);
        start_state_index(0) = 0;
        for(llong i = 1; i < n_models; i++) {
            start_state_index(i) = (2 * population_sizes(i - 1)) + start_state_index(i - 1);
        }

        for(llong i = 0; i < n_models; i++) {
            id.setZero();
            id(start_state_index(i)) = 1;
            N.row(i) = solver.solve(id, true);
            N.row(i) *= p(i);
        }

        double T_fix = N.sum();
        double rate = 1.0 / T_fix;

        dmat B(size, n_models);

        for(llong i = 0; i < n_models; i++) {
            dvec R_col = W.R.col(i);
            B.col(i) = solver.solve(R_col, false);
        }

        if(output_N_f) write_matrix_to_file(N, args::get(output_N_f));
        if(output_B_f) {
            write_matrix_to_file(B, args::get(output_B_f));
        }


        if (csv_f) {
            print_vector(population_sizes, "", ", ");
            print_vector(s, "", ", ");
            print_vector(h, "", ", ");
            print_vector(u, "", ", ");
            print_vector(v, "", ", ");
            print_vector(p, "", ", ");
            // print_vector(r, "", ", ");
            printf(DPF ", ", a);
            printf(DPF ", ", T_fix);
            printf(DPF "\n", rate);
        } else {
            print_vector(population_sizes, "N = ", "\n");

            print_vector(s, "s = ", "\n");
            print_vector(h, "h = ", "\n");
            print_vector(u, "u = ", "\n");
            print_vector(v, "v = ", "\n");
            print_vector(p, "p = ", "\n");
            // print_vector(r, "r = ", "\n");

            printf("a = " DPF "\n", a);
            printf("T_fix = " DPF "\n", T_fix);
            printf("Rate = " DPF "\n", rate);
        }

    } // END SWITCHING FIXATION }}}

    if(absorption_f) // BEGIN SWITCHING ABSORPTION {{{
    {
        WF::Matrix W = WF::Switching(population_sizes, WF::BOTH_ABSORBING, 
                s, h, u, v, switching, a, verbose_f);

        // SI are start indeces - a vector of size n_models
        lvec si = start_indeces(2 * population_sizes - lvec::Ones(n_models));

        if(output_Q_f) W.Q.save_market(args::get(output_Q_f));
        if(output_R_f) write_matrix_to_file(W.R, args::get(output_R_f));

        W.Q.subtract_identity();

        llong size = (2 * population_sizes.sum()) - n_models;

        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();

        // Get initial probabilities of mu within each model
        lvec nnz_p0(n_models);
        vector<dvec> p0(n_models);
        for (llong i = 0; i < n_models; i++) {
            llong pop_size = population_sizes(i);
            dvec first_row = WF::binom_row(2 * pop_size, WF::psi_diploid(0, pop_size, s(i), h(i), u(i), v(i)), a).Q;
            p0[i] = first_row.tail(first_row.size() - 1) / (1 - first_row(0)); // renormalize
            nnz_p0[i] = (p0[i].array() > integration_cutoff).count();
        }

        // extinction and fixation column for each submodel
        dmat B(size, n_models * 2);
        for (llong i = 0; i < n_models * 2; i++) {
            dvec R_col = W.R.col(i);
            B.col(i) = solver.solve(R_col, false);
        }

        map<llong, dvec> N_rows;
        map<llong, dvec> N2_rows;
        dvec id(size);
        for (llong i_ = 0; i_ < si.size(); i_++) {
            llong i = si[i_];
            for(llong o_ = 0; o_ < nnz_p0[i_]; o_++) {
                llong idx = i + o_;
                id.setZero();
                id(idx) = 1;
                N_rows[idx] = solver.solve(id, true);
                N2_rows[idx] = solver.solve(N_rows[idx], true);
            }
        }

        // absorbing extinction columns of B
        lvec ke = range_step(0, 2*n_models, 2);
        // absorbing fixation columns of B
        lvec kf = range_step(1, 2*n_models, 2);

        double P_ext = 0;
        double P_fix = 0;
        dvec P_cond_fix = dvec::Zero(n_models);
        dvec P_cond_ext = dvec::Zero(n_models);

        for (llong i_ = 0; i_ < si.size(); i_++) {
            llong i = si[i_];
            for(llong o_ = 0; o_ < nnz_p0[i_]; o_++) {
                double o = p0[i_](o_);
                llong idx = i + o_;

                double P_ext_i = 0;
                for (llong k_ = 0; k_ < ke.size(); k_++) {
                    P_ext_i += B(idx, ke[k_]);
                    P_cond_ext[k_] += B(idx, ke[k_]) * o * p[i_];
                }
                // P_cond_ext[i_] += P_ext_i * o;
                P_ext += P_ext_i * o * p[i_];

                double P_fix_i = 0;
                for (llong k_ = 0; k_ < kf.size(); k_++) {
                    P_fix_i += B(idx, kf[k_]);
                    P_cond_fix[k_] += B(idx, kf[k_]) * o * p[i_];
                }
                // P_cond_fix[i_] += P_fix_i * o;
                P_fix += P_fix_i * o * p[i_];
            }
        }

        // Summarize extinction and fixation absorption vectors
        dvec B_fix = dvec::Zero(size);
        dvec B_ext = dvec::Zero(size);
        for(llong k_ = 0; k_ < ke.size(); k_++) { B_ext += B.col(ke[k_]); }
        for(llong k_ = 0; k_ < kf.size(); k_++) { B_fix += B.col(kf[k_]); }

        double T_ext = 0;
        double T_fix = 0;
        double T_ext_var = 0;
        double T_fix_var = 0;

        dvec E_ext = dvec::Zero(size);
        dvec E_fix = dvec::Zero(size);
        dvec E_uncond = dvec::Zero(size);

        for (llong i_ = 0; i_ < si.size(); i_++) {
            llong i = si[i_];
            for(llong o_ = 0; o_ < nnz_p0[i_]; o_++) {
                double o = p0[i_](o_);
                llong idx = i + o_;
                double iw = o * p[i_]; // integration weight

                E_uncond += iw * N_rows[idx];

                dvec E_ext_i = B_ext.array() * N_rows[idx].array() / B_ext[idx];
                dvec E_ext_var_i = B_ext.array() * N2_rows[idx].array() / B_ext[idx];
                T_ext += iw * E_ext_i.sum();
                T_ext_var += (2 * E_ext_var_i.sum() - E_ext_i.sum() - pow(E_ext_i.sum(), 2)) * iw;
                E_ext += iw * E_ext_i;

                dvec E_fix_i = B_fix.array() * N_rows[idx].array() / B_fix[idx];
                dvec E_fix_var_i = B_fix.array() * N2_rows[idx].array() / B_fix[idx];
                T_fix += iw * E_fix_i.sum();
                T_fix_var += (2 * E_fix_var_i.sum() - E_fix_i.sum() - pow(E_fix_i.sum(), 2)) * iw;
                E_fix += iw * E_fix_i;
            }
        }

        // time spent in each model conditional on absorbing in a particuar state
        dvec T_cond_fix = dvec::Zero(n_models);
        dvec T_cond_ext = dvec::Zero(n_models);
        dvec T_uncond   = dvec::Zero(n_models);
        for(llong i_ = 0; i_ < si.size(); i_++){
            T_cond_ext[i_] = E_ext.segment(si[i_], (2 * population_sizes[i_]) - 1).sum();
            T_cond_fix[i_] = E_fix.segment(si[i_], (2 * population_sizes[i_]) - 1).sum();
            T_uncond[i_]   = E_uncond.segment(si[i_], (2 * population_sizes[i_]) - 1).sum();
        }

        double T_ext_std = sqrt(T_ext_var);
        double T_fix_std = sqrt(T_fix_var);

        // Output {{{
        if(output_N_ext_f) write_vector_to_file(E_ext, args::get(output_N_ext_f));
        if(output_N_fix_f) write_vector_to_file(E_fix, args::get(output_N_fix_f));
        if(output_N_f) write_vector_map_to_file(N_rows, args::get(output_N_f));
        if(output_B_f) write_matrix_to_file(B, args::get(output_B_f));

        if (csv_f) {
            print_vector(population_sizes, "", ", ");
            print_vector(s, "", ", ");
            print_vector(h, "", ", ");
            print_vector(u, "", ", ");
            print_vector(v, "", ", ");
            print_vector(p, "", ", ");
            // print_vector(r, "", ", ");
            printf(DPF ", ", a);
            printf(DPF ", ", P_ext);
            printf(DPF ", ", P_fix);
            printf(DPF ", ", T_ext);
            printf(DPF ", ", T_ext_std);
            printf(DPF ", ", T_fix);
            printf(DPF ", ", T_fix_std);
            print_vector(P_cond_ext, "", ", ");
            print_vector(P_cond_fix, "", ", ");
            print_vector(T_uncond, "", ", ");
            print_vector(T_cond_ext, "", ", ");
            print_vector(T_cond_fix, "", "\n");
        } else {
            print_vector(population_sizes, "N = ", "\n");
            print_vector(s, "s = ", "\n");
            print_vector(h, "h = ", "\n");
            print_vector(u, "u = ", "\n");
            print_vector(v, "v = ", "\n");
            print_vector(p, "p = ", "\n");
            // print_vector(r, "r = ", "\n");

            printf("a = " DPF "\n", a);
            printf("P_ext = " DPF "\n", P_ext);
            printf("P_fix = " DPF "\n", P_fix);
            printf("T_ext = " DPF "\n", T_ext);
            printf("T_ext_std = " DPF "\n", T_ext_std);
            printf("T_fix = " DPF "\n", T_fix);
            printf("T_fix_std = " DPF "\n", T_fix_std);
            print_vector(P_cond_ext, "P_cond_ext = ", "\n");
            print_vector(P_cond_fix, "P_cond_fix = ", "\n");
            print_vector(T_uncond, "T_uncond = ", "\n");
            print_vector(T_cond_ext, "T_cond_ext = ", "\n");
            print_vector(T_cond_fix, "T_cond_fix = ", "\n");
        } // }}}

    } // END SWITCHING ABSORPTION }}}

    if (verbose_f) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Total runtime: " << dt.count() << " s" << std::endl;
    }

    return EXIT_SUCCESS;
}
