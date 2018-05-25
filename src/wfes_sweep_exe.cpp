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
    args::ArgumentParser parser("WFES-SOFT");
    parser.helpParams.width = 120;
    parser.helpParams.helpindent = 50;
    parser.helpParams.flagindent = 2;

    args::ValueFlag<llong> population_size_f(parser, "int", "Size of the population", {'N', "pop-size"}, args::Options::Required);
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> backward_mutation_f(parser, "float[k]", "Backward mutation rates", {'u', "backward-mu"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> forward_mutation_f(parser, "float[k]", "Forward mutation rates", {'v', "forward-mu"});
    args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, 1>, NumericVectorReader<double>> selection_coefficient_f(parser, "float[k]", "Selection coefficients", {'s', "selection"}, args::Options::Required);
    // args::ValueFlag<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, NumericMatrixReader<double>> switching_f(parser, "float[k][k]", "Switching parameters over models", {'r', "switching"}, args::Options::Required);
    args::ValueFlag<double> lambda_f(parser, "float", "Transition probability", {'l', "lambda"}, args::Options::Required);
    args::ValueFlag<string> output_Q_f(parser, "path", "Output Q matrix to file", {"output-Q"});
    args::ValueFlag<string> output_N_f(parser, "path", "Output N matrix to file", {"output-N"});
    args::ValueFlag<double> alpha_f(parser, "float", "Tail truncation weight", {'a', "alpha"});
    args::Flag force_f(parser, "force", "Do not perform parameter checks", {"force"});
    args::Flag verbose_f(parser, "verbose", "Verbose solver output", {"verbose"});
    args::Flag csv_f(parser, "csv", "Output results in CSV format", {"csv"});

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
    double a = alpha_f ? args::get(alpha_f) : 1e-20;

    if(selection_coefficient.size() < 2) throw runtime_error("Population size vector should be longer than 2");

    if (!force_f) {
        if (population_size > 500000) {
            throw args::Error("Population size is quite large - the computations will take a long time. Use --force to ignore");   
        }
        double max_mu = max(u.maxCoeff(), v.maxCoeff());
        if (4 * population_size * max_mu > 1) {
            throw args::Error("The mutation rate might violate the Wright-Fisher assumptions. Use --force to ignore");
        }
        if (2 * population_size * selection_coefficient.minCoeff() < -10) {
            throw args::Error("The selection coefficient is quite negative. Fixations might be impossible. Use --force to ignore");
        }
        if (a > 1e-5) {
            throw args::Error("Zero cutoff value is quite high. This might produce inaccurate results. Use --force to ignore");   
        }
    }

    // dmat switching = switching_f ? args::get(switching_f) : dmat::Ones(1, 1);

    // dvec row_sums = switching.rowwise().sum();
    // for (llong i = 0; i < 2; i++) {
    //     for (llong j = 0; j < 2; j++) {
    //         switching(i, j) /= row_sums(i);
    //     }
    // }

    double l = args::get(lambda_f);
    dmat switching(2, 2); switching << 1 - l, l, 0, 1;


    dvec dom(2); dom << 0.5, 0.5;
    WF::Matrix wf = WF::NonAbsorbingToFixationOnly(population_size, selection_coefficient, dom, u, v, switching, a);
    if(output_Q_f) wf.Q.save_market(args::get(output_Q_f));


    wf.Q.subtract_identity();

    llong msg_level = verbose_f ? MKL_PARDISO_MSG_VERBOSE : MKL_PARDISO_MSG_QUIET;

    PardisoSolver solver(wf.Q, MKL_PARDISO_MATRIX_TYPE_REAL_UNSYMMETRIC, msg_level);
    solver.analyze();
    dvec id = dvec::Zero(wf.Q.n_row);
    id(0) = 1;

    dmat N(1, wf.Q.n_row);
    N.row(0) = solver.solve(id, true);

    double T_fix = N.sum();
    double rate = 1.0 / T_fix;

    if(output_N_f) write_matrix_to_file(N, args::get(output_N_f));

    if (csv_f) {
        printf("%lld, ", population_size);
        print_vector(selection_coefficient, "", ", ");
        print_vector(dom, "", ", ");
        print_vector(u, "", ", ");
        print_vector(v, "", ", ");
        printf(DPF ", ", a);
        printf(DPF ", ", T_fix);
        printf(DPF "\n", rate);
    } else {
        printf("N = %lld\n", population_size);
        print_vector(selection_coefficient, "s = ", "\n");
        print_vector(dom, "h = ", "\n");
        print_vector(u, "u = ", "\n");
        print_vector(v, "v = ", "\n");
        printf("a = " DPF "\n", a);
        printf("T_fix = " DPF "\n", T_fix);
        printf("Rate = " DPF "\n", rate);
    }


    return 0;
}
