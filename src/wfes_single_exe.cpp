#include "common.hpp"
#include "WrightFisher.hpp"
#include "PardisoSolver.hpp"
#include "parsing.hpp"
#include "args.hpp"

namespace WF = WrightFisher;
using namespace std;


int main(int argc, char const *argv[])
{
    
    args::ArgumentParser parser("WFES-SINGLE");
    
    args::Group model_f(parser, "Model type - specify one", args::Group::Validators::Xor, args::Options::Required);
    args::Flag absorption_f(model_f, "absorption", "Both fixation and extinction states are absorbing", {"absorption"});
    args::Flag fixation_f(model_f, "fixation", "Only fixation state is absorbing", {"fixation"});
    args::Flag allele_age_f(model_f, "allele-age", "Calculate age of an allele", {"allele-age"});

    args::ValueFlag<llong> population_size_f(parser, "int", "Size of the population", {'N', "pop-size"}, args::Options::Required);

    // Optional arguments
    args::ValueFlag<double> selection_coefficient_f(parser, "float", "Selection coefficient", {'s', "selection"});
    args::ValueFlag<double> dominance_f(parser, "float", "Dominance coefficient", {'h', "dominance"});
    args::ValueFlag<double> backward_mutation_f(parser, "float", "Backward mutation rate", {'u', "backward-mu"});
    args::ValueFlag<double> forward_mutation_f(parser, "float", "Forward mutation rate", {'v', "forward-mu"});
    args::ValueFlag<double> observed_count_f(parser, "float", "Observed allele count (allele-age only)", {'i', "observed-count"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::ValueFlag<double> integration_cutoff_f(parser, "float", "Starting number of copies integration cutoff", {'c', "integration-cutoff"});

    // Output options
    args::ValueFlag<string> output_Q_f(parser, "path", "Output Q matrix to file", {"output-Q"});
    args::ValueFlag<string> output_R_f(parser, "path", "Output R vectors to file", {"output-R"});
    args::ValueFlag<string> output_N_f(parser, "path", "Output N matrix to file", {"output-N"});
    args::ValueFlag<string> output_B_f(parser, "path", "Output B vectors to file", {"output-B"});
    args::ValueFlag<string> output_E_f(parser, "path", "Output Equilibrium frequencies to file (eigen only)", {"output-E"});
    
    args::Flag csv_f(parser, "csv", "Output results in CSV format", {"csv"});
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
    
    if(model_f.MatchedChildren() != 1) {
        throw args::Error("Should have exactly one of the 'Model type' options");
    }

    llong population_size = args::get(population_size_f);
    // Set default values
    double s = selection_coefficient_f ? args::get(selection_coefficient_f) : 0;
    double h = dominance_f ? args::get(dominance_f) : 0.5;
    double u = backward_mutation_f ? args::get(backward_mutation_f) : 1e-9;
    double v = forward_mutation_f ? args::get(forward_mutation_f) : 1e-9;
    double a = alpha_f ? args::get(alpha_f) : 1e-20;
    double integration_cutoff = integration_cutoff_f ? args::get(integration_cutoff_f) : 1e-10;

    llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    dvec first_row = WF::binom_row(0, population_size, population_size, s, h, u, v, a).Q;
    dvec starting_copies_p = first_row.tail(first_row.size() - 1); // renormalize
    starting_copies_p /= 1 - first_row(0);

    llong z = 0;

    if(integration_cutoff == -1) { // no integration
        z = 1;
        starting_copies_p[0] = 1;
    } else {
        for(llong i = 0; starting_copies_p(i) > integration_cutoff; i++, z++);
    }

    if(fixation_f) // BEGIN SINGLE FIXATION
    {
        WF::Matrix W = WF::SingleAlt(population_size, population_size, WF::FIXATION_ONLY, s, h, u, v, a);

        if(output_Q_f) W.Q.save_market(args::get(output_Q_f));
        if(output_R_f) write_matrix_to_file(W.R, args::get(output_R_f));

        W.Q.subtract_identity();

        llong size = (2 * population_size);

        PardisoSolver solver(W.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
        solver.analyze();

        dvec id(size);
        dmat N_mat(z, size);
        double T_fix = 0;

        // Iterate over starting states
        for(llong i = 0; i < z; i++) {
            id.setZero();
            id(i) = 1;

            N_mat.row(i) = solver.solve(id, true);

            T_fix += N_mat.row(i).sum();
            T_fix *= starting_copies_p(i);

        }
        double rate = 1.0 / T_fix;

        if(output_N_f) write_matrix_to_file(N_mat, args::get(output_N_f));
        if(output_B_f) {
            dvec B = dvec::Ones(size);
            write_vector_to_file(B, args::get(output_B_f));
        }

        if (csv_f) {
            printf("%lld, " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", "
                           DPF ", " DPF "\n", population_size, s, h, u, v, a, T_fix, rate);
        } else {
            printf("N = " LPF "\n", population_size);
            printf("s = " DPF "\n", s);
            printf("h = " DPF "\n", h);
            printf("u = " DPF "\n", u);
            printf("v = " DPF "\n", v);
            printf("a = " DPF "\n", a);
            printf("T_fix = " DPF "\n", T_fix);
            printf("Rate = " DPF "\n", rate);
        }
    } // END SINGLE FIXATION

    if(absorption_f) // BEGIN SINGLE ABSORPTION
    {
        WF::Matrix W = WF::SingleAlt(population_size, population_size, WF::BOTH_ABSORBING, s, h, u, v, a);

        if(output_Q_f) W.Q.save_market(args::get(output_Q_f));
        if(output_R_f) write_matrix_to_file(W.R, args::get(output_R_f));

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
        double T_ext = 0;
        double T_fix = 0;

        dmat N_mat(z, size);
        for(llong i = 0; i < z; i++) {
            id.setZero();
            id(i) = 1;

            N_mat.row(i) = solver.solve(id, true);
            // N_mat.row(i) *= starting_copies_p(i);

            P_ext += B_ext(i) * starting_copies_p(i);
            dvec E_ext = B_ext.transpose() * N_mat.row(i).transpose() / B_ext(i);
            T_ext += E_ext.sum() * starting_copies_p(i);

            P_fix += B_fix(i) * starting_copies_p(i);
            dvec E_fix = B_fix.transpose() * N_mat.row(i).transpose() / B_fix(i);
            T_fix += E_fix.sum() * starting_copies_p(i);
        }

        if(output_N_f) write_matrix_to_file(N_mat, args::get(output_N_f));
        if(output_B_f) {
            dmat B(size, 2);
            B.col(0) = B_ext;
            B.col(1) = B_fix;
            write_matrix_to_file(B, args::get(output_B_f));
        }

        if (csv_f) {
            printf("%lld, " DPF ", " DPF ", " DPF ", " DPF ", " DPF ", "
                           DPF ", " DPF ", " DPF ", " DPF "\n",
                   population_size, s, h, u, v, a, B_ext(0), B_fix(0), T_ext, T_fix);

        } else {
            printf("N = " LPF "\n", population_size);
            printf("s = " DPF "\n", s);
            printf("h = " DPF "\n", h);
            printf("u = " DPF "\n", u);
            printf("v = " DPF "\n", v);
            printf("a = " DPF "\n", a);
            printf("P_ext = " DPF "\n", B_ext(0));
            printf("P_fix = " DPF "\n", B_fix(0));
            printf("T_ext = " DPF "\n", T_ext);
            printf("T_fix = " DPF "\n", T_fix);
        }
    } // END SINGLE ABSORPTION

    if(allele_age_f) // BEGIN SINGLE ALLELE AGE
    {
        throw runtime_error("Unimplemented");
    } // END SINGLE ALLELE AGE

    return EXIT_SUCCESS;
}
