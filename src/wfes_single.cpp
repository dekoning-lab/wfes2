#include "PardisoSolver.hpp"
#include "WrightFisher.hpp"
#include "args.hpp"
#include "common.hpp"
#include "parsing.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
using namespace std;

int main(int argc, char const *argv[]) {

    args::ArgumentParser parser("WFES-SINGLE");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;

    args::Group model_f(parser, "Model type - specify one", args::Group::Validators::Xor,
                        args::Options::Required);
    args::Flag absorption_f(model_f, "absorption",
                            "Both fixation and extinction states are absorbing", {"absorption"});
    args::Flag fixation_f(model_f, "fixation", "Only fixation state is absorbing", {"fixation"});
    args::Flag establishment_f(model_f, "establishment", "Calculate establishment properties",
                               {"establishment"});
    args::Flag fundamental_f(model_f, "fundamental",
                             "Calculate the entire fundamental matrix (slow)", {"fundamental"});
    args::Flag equilibrium_f(model_f, "equilibrium",
                             "Calculate the equilibrium distribtion of allele states",
                             {"equilibrium"});
    args::Flag non_absorbing_f(model_f, "non-absorbing", "Build a non-absorbing WF matrix",
                               {"non-absorbing"});
    args::Flag allele_age_f(model_f, "allele-age", "Calculate age of an allele", {"allele-age"});

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
    args::Flag no_recurrent_mutation_f(parser, "bool", "Exclude recurrent mutation",
                                       {'m', "no-recurrent-mu"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::ValueFlag<llong> block_size_f(parser, "int", "Block size", {'b', "block-size"});
    args::ValueFlag<llong> n_threads_f(parser, "int", "Number of threads", {'t', "num-threads"});
    args::ValueFlag<double> integration_cutoff_f(parser, "float",
                                                 "Starting number of copies integration cutoff",
                                                 {'c', "integration-cutoff"});
    args::ValueFlag<string> initial_distributon_csv_f(
        parser, "path", "Path to initial probability distribution CSV", {'i', "initial"});
    args::ValueFlag<llong> starting_copies_f(
        parser, "int", "Starting number of copies - no integration", {'p', "starting-copies"});
    args::ValueFlag<llong> observed_copies_f(
        parser, "int", "Observed number of copies (--allele-age only)", {'x', "observed-copies"});
    args::ValueFlag<double> odds_ratio_f(parser, "float", "Odds ratio (--establishment only)",
                                         {'k', "odds-ratio"});

    // Output options
    args::ValueFlag<string> output_Q_f(parser, "path", "Output Q matrix to file", {"output-Q"});
    args::ValueFlag<string> output_R_f(parser, "path", "Output R vectors to file", {"output-R"});
    args::ValueFlag<string> output_N_f(parser, "path", "Output N matrix to file", {"output-N"});
    args::ValueFlag<string> output_N_ext_f(
        parser, "path", "Output extinction-conditional sojourn to file", {"output-N-ext"});
    args::ValueFlag<string> output_N_fix_f(
        parser, "path", "Output fixation-conditional sojourn to file", {"output-N-fix"});
    args::ValueFlag<string> output_B_f(parser, "path", "Output B vectors to file", {"output-B"});
    args::ValueFlag<string> output_I_f(parser, "path", "Output Initial probability distribution",
                                       {"output-I"});
    args::ValueFlag<string> output_E_f(
        parser, "path", "Output Equilibrium frequencies to file (--equilibrium only)",
        {"output-E"});
    args::ValueFlag<string> output_V_f(
        parser, "path", "Output Variance time matrix to file (--fundamental only)", {"output-V"});

    args::Flag csv_f(parser, "csv", "Output results in CSV format", {"csv"});
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

    time_point t_start, t_end;
    if (verbose_f)
        t_start = std::chrono::system_clock::now();

    if (model_f.MatchedChildren() != 1) {
        throw args::Error("Should have exactly one of the 'Model type' options");
    }

    llong population_size = args::get(population_size_f);
    // Set default values
    double s = selection_coefficient_f ? args::get(selection_coefficient_f) : 0;
    double h = dominance_f ? args::get(dominance_f) : 0.5;
    double u = backward_mutation_f ? args::get(backward_mutation_f) : 1e-9;
    double v = forward_mutation_f ? args::get(forward_mutation_f) : 1e-9;
    bool no_rem = no_recurrent_mutation_f ? args::get(no_recurrent_mutation_f) : false;
    bool rem = !no_rem;
    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    double odds_ratio = odds_ratio_f ? args::get(odds_ratio_f) : 1.0;
    double b = block_size_f ? args::get(block_size_f) : 100;
    double n_threads = n_threads_f ? args::get(n_threads_f) : 1;
    double integration_cutoff = integration_cutoff_f ? args::get(integration_cutoff_f) : 1e-10;

    // translate starting number of copies into model state (p - 1)
    llong starting_copies = starting_copies_f ? (args::get(starting_copies_f) - 1) : 0;

    if (!force_f) {
        if (population_size > 500000) {
            throw args::Error("Population size is quite large - the computations will take a long "
                              "time. Use --force to ignore");
        }
        double max_mu = max(u, v);
        if ((4 * population_size * max_mu) > 1) {
            throw args::Error("The mutation rate might violate the Wright-Fisher assumptions. Use "
                              "--force to ignore");
        }
        if ((2 * population_size * s) <= -100) {
            throw args::Error("The selection coefficient is quite negative. Fixations might be "
                              "impossible. Use --force to ignore");
        }
        if (a > 1e-5) {
            throw args::Error("Zero cutoff value is quite high. This might produce inaccurate "
                              "results. Use --force to ignore");
        }
    }

    llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

#ifdef OMP
    omp_set_num_threads(n_threads);
#endif
    mkl_set_num_threads(n_threads);

    dvec starting_copies_p;
    if (initial_distributon_csv_f) {
        // cout << "Reading initial from file" << args::get(initial_distributon_csv_f) << "" <<
        // endl;
        starting_copies_p = load_csv_col_vector(args::get(initial_distributon_csv_f));
    } else {
        dvec first_row =
            WF::binom_row(2 * population_size, WF::psi_diploid(0, population_size, s, h, u, v), a)
                .Q;
        starting_copies_p = first_row.tail(first_row.size() - 1); // renormalize
        starting_copies_p /= 1 - first_row(0);
    }

    if (output_I_f)
        write_vector_to_file(starting_copies_p, args::get(output_I_f));

    llong z = 0;

    if (initial_distributon_csv_f) {
        z = starting_copies_p.size();
    } else if (integration_cutoff <= 0 || v == 0) { // no integration
        z = 1;
        starting_copies_p[0] = 1;
    } else {
        for (llong i = 0; starting_copies_p(i) > integration_cutoff; i++) {
            z++;
        }
    }
    if (starting_copies_f)
        z = 1;

    if (fixation_f) // BEGIN SINGLE FIXATION
    {
        WF::Matrix W = WF::Single(population_size, population_size, WF::FIXATION_ONLY, s, h, u, v,
                                  rem, a, verbose_f, b);

        if (output_Q_f)
            W.Q.save_market(args::get(output_Q_f));
        if (output_R_f)
            write_matrix_to_file(W.R, args::get(output_R_f));

        W.Q.subtract_identity();

        llong size = (2 * population_size);

        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();

        dvec id(size);
        dmat N_mat(1, size);

        id.setZero();
        id(starting_copies) = 1;
        N_mat.row(0) = solver.solve(id, true);
        dvec N1 = N_mat.row(0);
        dvec N2 = solver.solve(N1, true);
        double T_fix = N1.sum();
        double T_var = ((2 * N2.sum()) - N1.sum()) - pow(N1.sum(), 2);

        double rate = 1.0 / T_fix;
        double T_std = sqrt(T_var);

        if (output_N_f)
            write_matrix_to_file(N_mat, args::get(output_N_f));
        if (output_B_f) {
            dvec B = dvec::Ones(size);
            write_vector_to_file(B, args::get(output_B_f));
        }

        if (csv_f) {
            printf("%lld, " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF "\n",
                   population_size, s, h, u, v, a, T_fix, T_std, rate);
        } else {
            printf("N = " LPF "\n", population_size);
            printf("s = " DPF "\n", s);
            printf("h = " DPF "\n", h);
            printf("u = " DPF "\n", u);
            printf("v = " DPF "\n", v);
            printf("a = " DPF "\n", a);
            printf("T_fix = " DPF "\n", T_fix);
            printf("T_std = " DPF "\n", T_std);
            printf("Rate = " DPF "\n", rate);
        }
    } // END SINGLE FIXATION

    if (absorption_f) // BEGIN SINGLE ABSORPTION
    {
        WF::Matrix W = WF::Single(population_size, population_size, WF::BOTH_ABSORBING, s, h, u, v,
                                  rem, a, verbose_f, b);

        if (output_Q_f)
            W.Q.save_market(args::get(output_Q_f));
        if (output_R_f)
            write_matrix_to_file(W.R, args::get(output_R_f));

        W.Q.subtract_identity();

        llong size = (2 * population_size) - 1;

        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();

        dvec R_ext = W.R.col(0);
        dvec B_ext = solver.solve(R_ext, false);
        dvec B_fix = dvec::Ones(size) - B_ext;
        dvec id(size);

        // integrate over starting number of copies
        double P_ext = 0;
        double P_fix = 0;
        double T_abs = 0;
        double T_ext = 0;
        double T_fix = 0;
        double T_abs_var = 0;
        double T_ext_var = 0;
        double T_fix_var = 0;
        double N_ext = 0;

        dmat N_mat(z, size);
        dmat E_ext_mat(z, size);
        dmat E_fix_mat(z, size);
        dmat N2_mat(z, size);
        if (!starting_copies_f) {
            for (llong i = 0; i < z; i++) {
                double p_i = starting_copies_p(i);
                id.setZero();
                id(i) = 1;

                N_mat.row(i) = solver.solve(id, true);
                dvec N1 = N_mat.row(i);
                N2_mat.row(i) = solver.solve(N1, true);
                dvec N2 = N2_mat.row(i);

                T_abs += N1.sum() * p_i;
                T_abs_var += (2 * N2.sum() - N1.sum()) - pow(N1.sum(), 2) * p_i;

                P_ext += B_ext(i) * p_i;
                dvec E_ext = B_ext.array() * N1.array() / B_ext(i);
                E_ext_mat.row(i) = E_ext;
                dvec E_ext_var = B_ext.array() * N2.array() / B_ext(i);
                T_ext += E_ext.sum() * p_i;
                T_ext_var += ((2 * E_ext_var.sum() - E_ext.sum()) - pow(E_ext.sum(), 2)) * p_i;

                dvec C_ext = E_ext.array() * dvec::LinSpaced(size, 1, size).array();
                N_ext += p_i * C_ext.sum();

                P_fix += B_fix(i) * p_i;
                dvec E_fix = B_fix.array() * N1.array() / B_fix(i);
                E_fix_mat.row(i) = E_fix;
                dvec E_fix_var = B_fix.array() * N2.array() / B_fix(i);
                T_fix += E_fix.sum() * p_i;
                T_fix_var += ((2 * E_fix_var.sum() - E_fix.sum()) - pow(E_fix.sum(), 2)) * p_i;

            }
        } else {
            // TODO: combine this with the previous clause
            id.setZero();
            id(starting_copies) = 1;
            N_mat.row(0) = solver.solve(id, true);
            dvec N1 = N_mat.row(0);
            N2_mat.row(0) = solver.solve(N1, true);
            dvec N2 = N2_mat.row(0);

            T_abs = N1.sum();
            T_abs_var = (2 * N2.sum() - N1.sum()) - pow(N1.sum(), 2);

            P_ext = B_ext(starting_copies);
            dvec E_ext = B_ext.array() * N1.array() / B_ext(starting_copies);
            E_ext_mat.row(0) = E_ext;
            dvec E_ext_var = B_ext.array() * N2.array() / B_ext(starting_copies);
            T_ext = E_ext.sum();
            T_ext_var = (2 * E_ext_var.sum() - E_ext.sum()) - pow(E_ext.sum(), 2);

            P_fix = B_fix(starting_copies);
            dvec E_fix = B_fix.array() * N1.array() / B_fix(starting_copies);
            E_fix_mat.row(0) = E_fix;
            dvec E_fix_var = B_fix.array() * N2.array() / B_fix(starting_copies);
            T_fix = E_fix.sum();
            T_fix_var = (2 * E_fix_var.sum() - E_fix.sum()) - pow(E_fix.sum(), 2);

            dvec C_ext = E_ext.array() * dvec::LinSpaced(size, 1, size).array();
            N_ext = C_ext.sum();
            // N_ext = (N_mat.row(0) * B_ext * dvec::LinSpaced(size, 1, size)).sum() /
            // B_ext(starting_copies);
        }

        double T_abs_std = sqrt(T_abs_var);
        double T_ext_std = sqrt(T_ext_var);
        double T_fix_std = sqrt(T_fix_var);

        N_ext /= (1 / (2 * population_size * v)) + T_ext;

        if (output_N_f)
            write_matrix_to_file(N_mat, args::get(output_N_f));
        if (output_N_ext_f)
            write_matrix_to_file(E_ext_mat, args::get(output_N_ext_f));
        if (output_N_fix_f)
            write_matrix_to_file(E_fix_mat, args::get(output_N_fix_f));
        if (output_B_f) {
            dmat B(size, 2);
            B.col(0) = B_ext;
            B.col(1) = B_fix;
            write_matrix_to_file(B, args::get(output_B_f));
        }

        if (csv_f) {
            printf("%lld, " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF
                   ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF "\n",
                   population_size, s, h, u, v, a, P_ext, P_fix, T_abs, T_abs_std, T_ext, T_ext_std, N_ext,
                   T_fix, T_fix_std);

        } else {
            printf("N = " LPF "\n", population_size);
            printf("s = " DPF "\n", s);
            printf("h = " DPF "\n", h);
            printf("u = " DPF "\n", u);
            printf("v = " DPF "\n", v);
            printf("a = " DPF "\n", a);
            printf("P_ext = " DPF "\n", P_ext);
            printf("P_fix = " DPF "\n", P_fix);
            printf("T_abs = " DPF "\n", T_abs);
            printf("T_abs_std = " DPF "\n", T_abs_std);
            printf("T_ext = " DPF "\n", T_ext);
            printf("T_ext_std = " DPF "\n", T_ext_std);
            printf("N_ext = " DPF "\n", N_ext);
            printf("T_fix = " DPF "\n", T_fix);
            printf("T_fix_std = " DPF "\n", T_fix_std);
            // printf("N_ext = " DPF "\n", N_ext);
        }
    } // END SINGLE ABSORPTION

    if (fundamental_f) {
        llong size = (2 * population_size) - 1;
        WF::Matrix W = WF::Single(population_size, population_size, WF::BOTH_ABSORBING, s, h, u, v,
                                  rem, a, verbose_f, b);
        if (output_Q_f)
            W.Q.save_market(args::get(output_Q_f));
        if (output_R_f)
            write_matrix_to_file(W.R, args::get(output_R_f));

        W.Q.subtract_identity();
        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();
        dmat N(size, size);
        dvec id(size);
        for (llong i = 0; i < size; i++) {
            id.setZero();
            id(i) = 1;
            N.row(i) = solver.solve(id, true);
        }
        if (output_N_f)
            write_matrix_to_file(N, args::get(output_N_f));

        if (output_V_f) {
            dvec Ndg = (2 * N.diagonal().array()) - 1;
            dmat Nsq = N.array().square();
            dmat V = (N * diagmat(Ndg)) - Nsq;

            write_matrix_to_file(V, args::get(output_V_f));
        }
    }

    if (equilibrium_f) {
        llong size = (2 * population_size) + 1;
        WF::Matrix W = WF::EquilibriumSolvingMatrix(population_size, s, h, u, v, a, verbose_f, b);
        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();
        dvec O = dvec::Zero(size);
        O(size - 1) = 1;

        dvec pi = solver.solve(O, true);
        write_vector_to_file(pi, args::get(output_E_f));

        // Calculate expected frequency
        double e_freq = 0.0;
        for (llong i = 0; i < size; i++) {
            e_freq += i * pi[i];
        }
        e_freq /= (size - 1);

        if (csv_f) {
            printf("%lld, " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF "\n",
                   population_size, s, h, u, v, a, e_freq, 1 - e_freq);

        } else {
            printf("N = " LPF "\n", population_size);
            printf("s = " DPF "\n", s);
            printf("h = " DPF "\n", h);
            printf("u = " DPF "\n", u);
            printf("v = " DPF "\n", v);
            printf("a = " DPF "\n", a);
            printf("E[freq mut] = " DPF "\n", e_freq);
            printf("E[freq  wt] = " DPF "\n", (1.0 - e_freq));
        }
    }

    if (establishment_f) {

        // Full Wright-Fisher
        WF::Matrix W_full = WF::Single(population_size, population_size, WF::BOTH_ABSORBING, s, h,
                                       u, v, rem, a, verbose_f, b);

        W_full.Q.subtract_identity();

        llong size = (2 * population_size) - 1;

        PardisoSolver solver_full(W_full.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver_full.analyze();

        dvec R_full_fix = W_full.R.col(1);
        dvec B_full_fix = solver_full.solve(R_full_fix, false);
        dvec B_full_ext = dvec::Constant(size, 1) - B_full_fix;
        dvec id_full(size);

        // establishment
        llong est_idx = 0;
        // find closest to k / (1 + k)
        // cout << odds_ratio / (1 + odds_ratio) << endl;
        dvec B_est_closest = B_full_fix - dvec::Constant(size, odds_ratio / (1 + odds_ratio));
        B_est_closest.array().abs().minCoeff(&est_idx);

        if (est_idx == 1) {
            throw args::Error("Establishment is near-certain: establishment-count is 1");
        }
        if (z >= est_idx) {
            throw args::Error("Establishment can be reached by mutation alone");
        }

        // Since the B indexes begin at 1
        est_idx++;
        // cout << est_idx << endl;
        double est_freq = (double)(est_idx) / (2 * population_size);

        // post-establishment time before absorption
        id_full.setZero();
        id_full(est_idx) = 1;
        dvec N1_aft_est = solver_full.solve(id_full, true);
        dvec N2_aft_est = solver_full.solve(N1_aft_est, true);

        // Segregation
        double T_seg = N1_aft_est.sum();
        double T_seg_var = (2 * N2_aft_est.sum() - N1_aft_est.sum()) - pow(N1_aft_est.sum(), 2);
        double T_seg_std = sqrt(T_seg_var);

        dvec E_seg_ext = B_full_ext.array() * N1_aft_est.array() / B_full_ext(est_idx);
        dvec E_seg_ext_var = B_full_ext.array() * N2_aft_est.array() / B_full_ext(est_idx);
        double T_seg_ext = E_seg_ext.sum();
        double T_seg_ext_var =
            (2 * E_seg_ext_var.sum() - E_seg_ext.sum()) - pow(E_seg_ext.sum(), 2);
        double T_seg_ext_std = sqrt(T_seg_ext_var);

        dvec E_seg_fix = B_full_fix.array() * N1_aft_est.array() / B_full_fix(est_idx);
        dvec E_seg_fix_var = B_full_fix.array() * N2_aft_est.array() / B_full_fix(est_idx);
        double T_seg_fix = E_seg_fix.sum();
        double T_seg_fix_var =
            (2 * E_seg_fix_var.sum() - E_seg_fix.sum()) - pow(E_seg_fix.sum(), 2);
        double T_seg_fix_std = sqrt(T_seg_fix_var);

        // Truncated model
        WF::Matrix W_tr = WF::Truncated(population_size, population_size, est_idx, s, h, u, v, rem,
                                        a, verbose_f, b);
        if (output_Q_f)
            W_tr.Q.save_market(args::get(output_Q_f));
        if (output_R_f)
            write_matrix_to_file(W_tr.R, args::get(output_R_f));

        // To test
        // cout << W_tr.R.col(0) + W_tr.Q.dense().rowwise().sum() + W_tr.R.col(1) << endl;
        // cout << W.Q.dense() << endl;

        W_tr.Q.subtract_identity();

        PardisoSolver solver_tr(W_tr.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver_tr.analyze();

        dvec R_est = W_tr.R.col(1);
        dvec B_est = solver_tr.solve(R_est, false);
        dvec B_ext = dvec::Ones(est_idx - 1) - B_est;

        // integrate over starting number of copies
        double P_ext = 0;
        double P_est = 0;
        double T_ext = 0;
        double T_est = 0;
        double T_ext_var = 0;
        double T_est_var = 0;

        dmat N_mat(z, est_idx - 1);
        dmat E_ext_mat(z, est_idx - 1);
        dmat E_est_mat(z, est_idx - 1);
        dmat N2_mat(z, est_idx - 1);

        dvec id(est_idx);
        if (!starting_copies_f) {
            for (llong i = 0; i < z; i++) {
                double p_i = starting_copies_p(i);
                id.setZero();
                id(i) = 1;

                N_mat.row(i) = solver_tr.solve(id, true);
                dvec N1 = N_mat.row(i);
                N2_mat.row(i) = solver_tr.solve(N1, true);
                dvec N2 = N2_mat.row(i);

                P_ext += B_ext(i) * p_i;
                dvec E_ext = B_ext.array() * N1.array() / B_ext(i);
                E_ext_mat.row(i) = E_ext;
                dvec E_ext_var = B_ext.array() * N2.array() / B_ext(i);
                T_ext += E_ext.sum() * p_i;
                T_ext_var +=
                    (((2 * E_ext_var.sum() - E_ext.sum()) * p_i) - pow(E_ext.sum() * p_i, 2));

                P_est += B_est(i) * p_i;
                dvec E_est = B_est.array() * N1.array() / B_est(i);
                E_est_mat.row(i) = E_est;
                dvec E_est_var = B_est.array() * N2.array() / B_est(i);
                T_est += E_est.sum() * p_i;
                T_est_var +=
                    (((2 * E_est_var.sum() - E_est.sum()) * p_i) - pow(E_est.sum() * p_i, 2));
            }
        } else {
            // TODO: combine this with the previous clause
            id.setZero();
            id(starting_copies) = 1;
            N_mat.row(0) = solver_tr.solve(id, true);
            dvec N1 = N_mat.row(0);
            N2_mat.row(0) = solver_tr.solve(N1, true);
            dvec N2 = N2_mat.row(0);

            P_ext = B_ext(starting_copies);
            dvec E_ext = B_ext.array() * N1.array() / B_ext(starting_copies);
            E_ext_mat.row(0) = E_ext;
            dvec E_ext_var = B_ext.array() * N2.array() / B_ext(starting_copies);
            T_ext = E_ext.sum();
            T_ext_var = (2 * E_ext_var.sum() - E_ext.sum()) - pow(E_ext.sum(), 2);

            P_est = B_est(starting_copies);
            dvec E_est = B_est.array() * N1.array() / B_est(starting_copies);
            E_est_mat.row(0) = E_est;
            dvec E_est_var = B_est.array() * N2.array() / B_est(starting_copies);
            T_est = E_est.sum();
            T_est_var = (2 * E_est_var.sum() - E_est.sum()) - pow(E_est.sum(), 2);
        }
        double T_est_std = sqrt(T_est_var);

        if (csv_f) {
            printf("%lld, " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF
                   "," DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF "\n",
                   population_size, s, h, u, v, odds_ratio, a, est_freq, P_est, T_seg, T_seg_std,
                   T_seg_ext, T_seg_ext_std, T_seg_fix, T_seg_fix_std, T_est, T_est_std);

        } else {
            printf("N = " LPF "\n", population_size);
            printf("s = " DPF "\n", s);
            printf("h = " DPF "\n", h);
            printf("u = " DPF "\n", u);
            printf("v = " DPF "\n", v);
            printf("odds_ratio = " DPF "\n", odds_ratio);
            printf("a = " DPF "\n", a);
            printf("F_est = " DPF "\n", est_freq);
            printf("P_est = " DPF "\n", P_est);
            printf("T_seg = " DPF "\n", T_seg);
            printf("T_seg_std = " DPF "\n", T_seg_std);
            printf("T_seg_ext = " DPF "\n", T_seg_ext);
            printf("T_seg_ext_std = " DPF "\n", T_seg_ext_std);
            printf("T_seg_fix = " DPF "\n", T_seg_fix);
            printf("T_seg_fix_std = " DPF "\n", T_seg_fix_std);
            printf("T_est = " DPF "\n", T_est);
            printf("T_est_std = " DPF "\n", T_est_std);
        }
    }

    if (allele_age_f) // BEGIN SINGLE ALLELE AGE
    {
        if (!observed_copies_f) {
            throw args::Error("-x | --observed-copies required for allele are calculation");
        }
        llong x = args::get(observed_copies_f) - 1;

        llong size = (2 * population_size) - 1;
        WF::Matrix W = WF::Single(population_size, population_size, WF::BOTH_ABSORBING, s, h, u, v,
                                  rem, a, verbose_f, b);
        if (output_Q_f)
            W.Q.save_market(args::get(output_Q_f));
        if (output_R_f)
            write_matrix_to_file(W.R, args::get(output_R_f));
        dvec Q_x = W.Q.col_copy(x);
        W.Q.subtract_identity();
        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();

        W.Q.subtract_identity();
        dvec Q_I_x = W.Q.col_copy(x);
        Q_I_x(x) += 1;
        dvec A_x = W.Q.multiply(Q_I_x);

        double E_allele_age = 0;
        double S_allele_age = 0;
        if (!starting_copies_f) {
            // Iterate over starting states
            for (llong i = 0; i < z; i++) {
                dvec e_p = dvec::Zero(size);
                e_p(i) = 1;

                dvec M1 = solver.solve(e_p, true);
                dvec M2 = solver.solve(M1, true);

                double mu1 = M2.dot(Q_x) / M1(x);

                dvec M3 = solver.solve(M2, true);

                double mu2 = sqrt((M3.dot(A_x) / M1(x)) - pow(mu1, 2));

                E_allele_age += mu1 * starting_copies_p(i);
                S_allele_age += mu2 * starting_copies_p(i);
            }
        } else {
            dvec e_p = dvec::Zero(size);
            e_p(starting_copies) = 1;

            dvec M1 = solver.solve(e_p, true);
            dvec M2 = solver.solve(M1, true);

            E_allele_age = M2.dot(Q_x) / M1(x);

            dvec M3 = solver.solve(M2, true);

            S_allele_age = sqrt((M3.dot(A_x) / M1(x)) - pow(E_allele_age, 2));
        }

        if (csv_f) {
            printf("%lld, " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", " DPF "\n",
                   population_size, s, h, u, v, a, E_allele_age, S_allele_age);

        } else {
            printf("N = " LPF "\n", population_size);
            printf("s = " DPF "\n", s);
            printf("h = " DPF "\n", h);
            printf("u = " DPF "\n", u);
            printf("v = " DPF "\n", v);
            printf("a = " DPF "\n", a);
            printf("E(A) = " DPF "\n", E_allele_age);
            printf("S(A) = " DPF "\n", S_allele_age);
        }
    } // END SINGLE ALLELE AGE

    if (non_absorbing_f) {
        WF::Matrix W = WF::Single(population_size, population_size, WF::NON_ABSORBING, s, h, u, v,
                                  rem, a, verbose_f, b);
        if (output_Q_f)
            W.Q.save_market(args::get(output_Q_f));
    }

    if (verbose_f) {
        t_end = std::chrono::system_clock::now();
        time_diff dt = t_end - t_start;
        std::cout << "Total runtime: " << dt.count() << " s" << std::endl;
    }

    return EXIT_SUCCESS;
}
