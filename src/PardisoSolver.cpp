#include "PardisoSolver.hpp"


PardisoSolver::PardisoSolver(SparseMatrix& A, llong matrix_type, llong message_level, llong n_rhs):
        size(A.n_row),
        n_right_hand_sides(n_rhs),
        max_factors(1),
        matrix_type(matrix_type),
        matrix_number(1),
        error(0),
        message_level(message_level),
        m(A),
        control(MKL_IFS_SIZE),
        internal(MKL_IFS_SIZE),
        workspace((ulong)(size * n_rhs))
{
    control(MKL_PARDISO_DEFAULT_SETTINGS) = MKL_PARDISO_FALSE;
    control(MKL_PARDISO_INDEXING_OPTION) = MKL_PARDISO_INDEXING_ZERO;
    control(MKL_PARDISO_FILL_IN_REDUCING_ORDERING_OPTION) = MKL_PARDISO_FILL_IN_REDUCING_ORDERING_NESTED_DISSECTION_OMP;
    control(MKL_PARDISO_ITERATIVE_REFINEMENT_MAX) = 0;
    control(MKL_PARDISO_PIVOTING_PERTURBATION) = 20; // Perturb the pivot elements with 1E-20
    control(MKL_PARDISO_SCALING_OPTION) = MKL_PARDISO_SCALING_ENABLE; 
    control(MKL_PARDISO_SOLVE_OPTION) = MKL_PARDISO_DEFAULT;
    control(MKL_PARDISO_WEIGHTED_MATCHING_OPTION) = MKL_PARDISO_WEIGHTED_MATCHING_ENABLE;
    control(MKL_PARDISO_PRECISION_OPTION) = MKL_PARDISO_PRECISION_DOUBLE;
    control(MKL_PARDISO_INDEXING_OPTION) = MKL_PARDISO_INDEXING_ZERO;
    control(MKL_PARDISO_OOC_OPTION) = MKL_PARDISO_OOC_OVERFLOW;
    control(MKL_PARDISO_REPORT_NNZ_FACTORS) = MKL_PARDISO_REPORT_ENABLE;
    control(MKL_PARDISO_REPORT_FLOP_FACTOR_PHASE) = MKL_PARDISO_REPORT_ENABLE;
    control(MKL_PARDISO_REPORT_CGS_CG_DIAGNOSTIC) = MKL_PARDISO_REPORT_ENABLE;
    control(MKL_PARDISO_MATRIX_CHECK_OPTION) = MKL_PARDISO_MATRIX_CHECK_ENABLE;
    control(MKL_PARDISO_PARALLEL_FACTORIZATION_OPTION) = MKL_PARDISO_PARALLEL_FACTORIZATION_TWO_LEVEL;
    control(MKL_PARDISO_PIVOT_OPTION) = MKL_PARDISO_PIVOT_CALLBACK;
}

void PardisoSolver::analyze()
{
    phase = MKL_PARDISO_SOLVER_PHASE_ANALYSIS;

    pardiso_64(internal.values, &max_factors, &matrix_number,
               &matrix_type, &phase, &size,
               m.values, m.row_index, m.columns,
               nullptr, &n_right_hand_sides, control.values,
               &message_level, nullptr, nullptr, &error);

    if(error != 0) throw runtime_error("PardisoSolver::analyze(): Symbolic factorization error: " + to_string(error));

    phase = MKL_PARDISO_SOLVER_PHASE_NUMERICAL_FACTORIZATION;
    pardiso_64(internal.values, &max_factors, &matrix_number,
               &matrix_type, &phase, &size,
               m.values, m.row_index, m.columns,
               nullptr, &n_right_hand_sides, control.values,
               &message_level, nullptr, nullptr, &error);

    if(error != 0) throw runtime_error("PardisoSolver::analyze(): Numerical factorization error: " + to_string(error));
}

dvec PardisoSolver::solve(dvec& b, bool transpose)
{
    phase = MKL_PARDISO_SOLVER_PHASE_SOLVE_ITERATIVE_REFINEMENT;
    if(transpose) control(MKL_PARDISO_SOLVE_OPTION) = MKL_PARDISO_SOLVE_TRANSPOSED;
    else control(MKL_PARDISO_SOLVE_OPTION) = MKL_PARDISO_DEFAULT;

    pardiso_64(internal.values, &max_factors, &matrix_number,
               &matrix_type, &phase, &size,
               m.values, m.row_index, m.columns,
               nullptr, &n_right_hand_sides, control.values,
               &message_level, b.values, workspace.values, &error);

    if(error != 0) throw runtime_error("PardisoSolver::solve(): Solution error: " + to_string(error));

    dvec x(size);
    for(llong i = 0; i < size; i++) x(i) = workspace(i);
    return x;
}

dvec PardisoSolver::get_diagonal() 
{
  dvec df(size);
  dvec da(size);

  pardiso_getdiag(internal.values, df.values, da.values, &matrix_number, &error);

  if(error == 1) throw runtime_error("PardisoSolver::get_diagonal(): Diagonal information not turned on before pardiso main loop: " + to_string(error));

  return df;
}

PardisoSolver::~PardisoSolver()
{
    phase = MKL_PARDISO_SOLVER_PHASE_RELEASE_MEMORY_ALL;
    pardiso_64(internal.values, &max_factors, &matrix_number,
               &matrix_type, &phase, &size,
               nullptr, m.row_index, m.columns,
               nullptr, &n_right_hand_sides, control.values,
               &message_level, nullptr, nullptr, &error);

}
