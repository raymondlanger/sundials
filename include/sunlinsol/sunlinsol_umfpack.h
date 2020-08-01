/*
 * -----------------------------------------------------------------
 * Programmer:  Raymond Langer @ RWTH Aachen University
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 * This is the header file for the UMFPACK implementation of the
 * SUNLINSOL module.
 *
 * UMFPACK is a set of routines for solving unsymmetric sparse 
 * linear systems.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file 
 *     include/sundials/sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_UMFPACK_H
#define _SUNLINSOL_UMFPACK_H

#include <sundials/sundials_config.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "umfpack.h"

/* Interfaces to match 'sunindextype' with the correct UMFPACK types/functions
 */
#if defined(SUNDIALS_INT64_T)
#error INT64 sunindextype no implemented for UMFPACK 
#elif defined(SUNDIALS_INT32_T)
#else /* incompatible sunindextype for UMFPACK */
#error Incompatible sunindextype for UMFPACK
#endif

#if defined(SUNDIALS_DOUBLE_PRECISION)
#else
#error Incompatible realtype for UMFPACK
#endif

struct _SUNLinearSolverContent_UMFPACK {
  /* indicates whether the last sundials linear solver routine succeeded */ 
  long int last_flag;

  int first_factorize;
  /*  Control [UMFPACK_CONTROL] is an array controling the symbolic analysis,
   *  the numerical factorization, and the solve step from the umfpack library.
   *  
   *  Input argument, not modified by the linear solver
   *
   *  Only the primary parameters are listed below:
   *
   *  GENERAL CONTROLS:
   *  Control [UMFPACK_PRL]: 
   *     - umfpack * report status:
   *     No output if the print level is 0 or less, even when an error occurs. 
   *     If 1, then error messages are printed, and nothing is printed if the 
   *     status is UMFPACK OK. A warning message is printedif the matrix is 
   *     singular. If 2 or more, then the status is always printed. If 4 or 
   *     more, then the UMFPACK Copyright is printed. If 6 or more, then the 
   *     UMFPACK License is printed.
   *     See also the first page of this User Guide for the Copyright and 
   *     License.
   *     - umfpack * report control:
   *     No output if the print level is 1 or less. If 2 or more, all of Control 
   *     is printed.
   *     - umfpack * report info:
   *     No output if the print level is 1 or less. If 2 or more, all of Info 
   *     is printed.
   *     - all other umfpack * report * routines:
   *     If the print level is 2 or less, then these routines return silently
   *     without checking their inputs.
   *     If 3 or more, the inputs are fully verified and a short status summary
   *     is printed.
   *     If 4, then the first few  entries of the input arguments are printed. 
   *     If 5, then all of the input arguments are printed.
   *
   *
   *  SYMBOLIC FACTORIZATION: 
   *
   *  Control [UMFPACK_STRATEGY]: 
   *     This is the most important control parameter. It determines what kind 
   *     of ordering and pivoting strategy that UMFPACK should use. It is new 
   *     to Version 4.1. There are 4 options:
   *     1) UMFPACK_STRATEGY_AUTO: This is the default. The input matrix is
   *        analyzed to determine how symmetric the nonzero pattern is, and
   *        how many entries there are on the diagonal. It then selects one
   *        of the following strategies. Refer to the User Guide for a
   *        description of how the strategy is automatically selected.
   *
   *     2) UMFPACK_STRATEGY_UNSYMMETRIC: Use the unsymmetric strategy. COLAMD
   *        is used to order the columns of A, followed by a postorder of
   *        the column elimination tree. No attempt is made to perform diagonal
   *        pivoting. The column ordering is refined during factorization. This 
   *        strategy was the only one provided with UMFPACK V4.0.
   *
   *        In the numerical factorization, the 
   *        Control[UMFPACK_SYM_PIVOT_TOLERANCE] parameter is ignored. A
   *        pivot is selected if its magnitude is >=
   *        Control[UMFPACK_PIVOT_TOLERANCE] (default 0.1) times the
   *        largest entry in its column.
   *
   *     3) UMFPACK_STRATEGY_SYMMETRIC: Use the symmetric strategy (new to
   *        Version 4.1). In this method, the approximate minimum degree
   *        ordering (AMD) is applied to A+A’, followed by a postorder of
   *        the elimination tree of A+A’. UMFPACK attempts to perform
   *        diagonal pivoting during numerical factorization. No refinement
   *        of the column preordering is performed during factorization.
   *        In the numerical factorization, a nonzero entry on the diagonal
   *        is selected as the pivot if its magnitude is 
   *        >= Control[UMFPACK_SYM_PIVOT_TOLERANCE] (default 0.001) times the
   *        largest entry in its column. If this is not acceptable, then an
   *        off-diagonal pivot is selected with magnitude >= Control
   *        [UMFPACK_PIVOT_TOLERANCE] (default 0.1) times the largest entry
   *        in its column.
   *
   *     4) Not provided in the manuals I've copied...
   *
   *  Control [UMFPACK_ORDERING]: The ordering method to use:
   *     UMFPACK_ORDERING_CHOLMOD try AMD/COLAMD, then METIS if needed
   *     UMFPACK_ORDERING_AMD just AMD or COLAMD
   *     UMFPACK_ORDERING_GIVEN just Qinit (umfpack_*_qsymbolic only)
   *     UMFPACK_ORDERING_NONE no fill-reducing ordering
   *     UMFPACK_ORDERING_METIS just METIS(A+A’) or METIS(A’A)
   *     UMFPACK_ORDERING_BEST try AMD/COLAMD, METIS, and NESDIS
   *     UMFPACK_ORDERING_USER just user function (*_fsymbolic only)
   *
   *  Control [UMFPACK_SCALE]: This parameter is new to V4.1. See
   *        umfpack_numeric.h for a description. Only affects the 2-by-2
   *        strategy. Default: UMFPACK_SCALE_SUM. 
   *
   *
   *  NUMERICAL FACTORIZATION: 
   *
   *  Control[UMFPACK_PIVOT_TOLERANCE]: relative pivot tolerance for
   *        threshold partial pivoting with row interchanges. In any given
   *        column, an entry is numerically acceptable if its absolute value is
   *        greater than or equal to Control [UMFPACK_PIVOT_TOLERANCE] times
   *        the largest absolute value in the column. A value of 1.0 gives true
   *        partial pivoting. If less than or equal to zero, then any nonzero
   *        entry is numerically acceptable as a pivot (this is changed from
   *        Version 4.0). Default: 0.1.
   *
   *  Control[UMFPACK_SYM_PIVOT_TOLERANCE]: This parameter is new to V4.1.
   *        If diagonal pivoting is attempted (the symmetric
   *        strategy is used) then this parameter is used to control when the
   *        diagonal entry is selected in a given pivot column. The absolute
   *        value of the entry must be >= Control [UMFPACK_SYM_PIVOT_TOLERANCE]
   *        times the largest absolute value in the column. A value of zero
   *        will ensure that no off-diagonal pivoting is performed, except that
   *        zero diagonal entries are not selected if there are any off-diagonal
   *        nonzero entries.
   *
   *  Control[UMFPACK_SCALE]: This parameter is new to V4.1. Version 4.0
   *        did not scale the matrix. Note that the user’s input matrix is
   *        never modified, only an internal copy is scaled.
   *        1) UMFPACK_SCALE_NONE: no scaling is performed.
   *        2) UMFPACK_SCALE_SUM: each row of the input matrix A is divided by
   *           the sum of the absolute values of the entries in that row.
   *           The scaled matrix has an infinity norm of 1.
   *        3) UMFPACK_SCALE_MAX: each row of the input matrix A is divided by
   *           the maximum the absolute values of the entries in that row.
   *           In the scaled matrix the largest entry in each row has
   *           a magnitude exactly equal to 1.
   *        Scaling is very important for the "symmetric" strategy when
   *        diagonal pivoting is attempted. It also improves the performance
   *        of the "unsymmetric" strategy.
   *        Default: UMFPACK_SCALE_SUM.
   *
   *
   *  SOLVE STEP:
   *
   *  Control [UMFPACK_IRSTEP]: 
   *        The maximum number of iterative refinement
   *        steps to attempt. A value less than zero is treated as zero. If
   *        less than 1, or if Ax=b or A’x=b is not being solved, or
   *        if A is singular, then the Ap, Ai, and Ax arguments are not
   *        accessed. 
   *        Default: 2.
   *
   */

   double Control[UMFPACK_CONTROL];

  /*  Info[UMFPACK_INFO] is an array containing statistics about the 
   *  symbolic analysis, the numeric factorization, or the solve step
   *
   *
   *  SYMBOLIC FACTORIZATION: 
   *
   *  Info[UMFPACK_STATUS]: 
   *     status code. This is also the return value of a symbolic 
   *     factorization
   *  Info [UMFPACK_SIZE_OF_UNIT]: 
   *     the number of bytes in a Unit, for memory usage statistics
   *     below 
   *  Info [UMFPACK_SYMBOLIC_PEAK_MEMORY]: 
   *     the amount of memory (in Units) required for umfpack_di_symbolic 
   *     to complete. This count includes the size of the Symbolic object 
   *     itself, which is also reported in Info [UMFPACK_SYMBOLIC_SIZE].
   *  Info [UMFPACK_NUMERIC_SIZE_ESTIMATE]: 
   *     an estimate of the final size (in Units) of the entire Numeric 
   *     object (both fixed-size and variablesized parts), which holds 
   *     the LU factorization (including the L, U, P and Q matrices).
   *      The count includes the size of both the Symbolic and Numeric 
   *      objects themselves. It can be a very loose upper bound,
   *      particularly when the symmetric or 2-by-2 strategies are used.
   *  Info [UMFPACK_PEAK_MEMORY_ESTIMATE], was computed by
   *      umfpack_di_symbolic. The estimate is normally an upper bound on the
   *      actual peak usage, but this is not guaranteed. With testing on
   *      hundreds of matrix arising in real applications, I have never
   *      observed a matrix where this estimate or the Numeric size estimate
   *      was less than the actual result, but this is theoretically possible.
   *      Please send me one if you find such a matrix.
   *  Info [UMFPACK_FLOPS_ESTIMATE]: 
   *      an estimate of the total floating-point operations required to 
   *      factorize the matrix. This is a "true" theoretical estimate of the 
   *      number of flops that would be performed by a flop-parsimonious sparse 
   *      LU algorithm. It assumes that no extra flops are performed except for 
   *      what is strictly required to compute the LU factorization. It ignores,
   *      for example, the flops performed by umfpack_di_numeric to add 
   *      contribution blocks of frontal matrices together.
   *      The actual "true flop" count found by umfpack_di_numeric will be less
   *      than this estimate.
   *  Info [UMFPACK_LNZ_ESTIMATE]: 
   *      an estimate of the number of nonzeros in L, including the diagonal. 
   *      Since L is unit-diagonal, the diagonal of L is not stored. This 
   *      estimate is a strict upper bound on the actual nonzeros in L to be 
   *      computed by umfpack_di_numeric.
   *  Info [UMFPACK_UNZ_ESTIMATE]: 
   *      an estimate of the number of nonzeros in U, including the diagonal. 
   *      This estimate is a strict upper bound on the actual nonzeros in U to 
   *      be computed by umfpack_di_numeric.
   *  Info [UMFPACK_SYMBOLIC_TIME]: 
   *      The CPU time taken, in seconds.
   *  Info [UMFPACK_STRATEGY_USED]: 
   *      The ordering strategy used: 
   *      UMFPACK_STRATEGY_SYMMETRIC or UMFPACK_STRATEGY_UNSYMMETRIC
   *
   *
   *  NUMERICAL FACTORIZATION: 
   *
   *  Info[UMFPACK_NUMERIC_SIZE]: the actual final size (in Units) of the
   *      entire Numeric object, including the final size of the variable
   *      part of the object. Info [UMFPACK_NUMERIC_SIZE_ESTIMATE],
   *      an estimate, was computed by umfpack_di_symbolic. The estimate is
   *      normally an upper bound on the actual final size, but this is not
   *      guaranteed.
   *  Info [UMFPACK_PEAK_MEMORY]: the actual peak memory usage (in Units) of
   *      both umfpack_di_symbolic and umfpack_di_numeric. An estimate,
   *  Info [UMFPACK_FLOPS]: the actual count of the (useful) floating-point
   *      operations performed. An estimate, Info [UMFPACK_FLOPS_ESTIMATE],
   *      was computed by umfpack_di_symbolic. The estimate is guaranteed to
   *      be an upper bound on this flop count. The flop count excludes
   *      "useless" flops on zero values, flops performed during the pivot
   *      search (for tentative updates and assembly of candidate columns),
   *      and flops performed to add frontal matrices together.
   *  Info [UMFPACK_LNZ]: the actual nonzero entries in final factor L,
   *      including the diagonal. This excludes any zero entries in L,
   *      although some of these are stored in the Numeric object. The
   *      Info [UMFPACK_LU_ENTRIES] statistic does account for all
   *      explicitly stored zeros, however. Info [UMFPACK_LNZ_ESTIMATE]
   *      an estimate, was computed by umfpack_di_symbolic. The estimate is
   *      guaranteed to be an upper bound on Info [UMFPACK_LNZ].
   *  Info [UMFPACK_UNZ]: the actual nonzero entries in final factor U,
   *      including the diagonal. This excludes any zero entries in U,
   *      although some of these are stored in the Numeric object. The
   *      Info [UMFPACK_LU_ENTRIES] statistic does account for all
   *      explicitly stored zeros, however. Info [UMFPACK_UNZ_ESTIMATE],
   *      an estimate, was computed by umfpack_di_symbolic. The estimate is
   *      guaranteed to be an upper bound on Info [UMFPACK_UNZ].
   *  Info [UMFPACK_NUMERIC_TIME]: The CPU time taken, in seconds.
   *
   *
   *  SOLVE STEP:
   *
   *  Info [UMFPACK_SOLVE_FLOPS]:
   *      the number of floating point operations
   *      performed to solve the linear system. This includes the work
   *      taken for all iterative refinement steps, including the backtrack
   *      (if any).
   *   Info [UMFPACK_SOLVE_TIME]:
   *      The time taken, in seconds.
   *   */

  double Info[UMFPACK_INFO];

  /*  Symbolic is a (void *) pointer variable. The address of Symbolic is 
   *  provided to the symbolic factorization routine umfpack_di_symbolic. On 
   *  input to umfpack_di_symbolic, the contents of this variable are not 
   *  defined. On output, this variable holds a (void *) pointer to the 
   *  Symbolic object (if successful), or (void *) NULL if a failure occurred.
   *
   *  After calling umfpack_di_symbolic, the Symbolic object holds the 
   *  symbolic factorization. The Symbolic object is then provided to
   *  umfpack_di_numeric as a const input.
   *   */
   void *Symbolic;

  /*  Numeric is a (void *) pointer variable. The address of Numeric is 
   *  provided to the numerical factorization routine umfpack_di_numerical. On 
   *  input to umfpack_di_numeric, the contents of this variable are not 
   *  defined. On output, this variable holds a (void *) pointer to the 
   *  Numeric object (if successful), or (void *) NULL if a failure occurred.
   *
   *  After calling umfpack_di_numeric, the Numeric object holds the 
   *  numerical factorization. The Numeric object is then provided to
   *  umfpack_di_solve as a const input.
   *   */
   void *Numeric;


  /* The error indicator: Meaning depends on whether last step was 
   *                      a symbolic, numeric, or solve step
   *                      use it together with last_umfpack_step
   *
   *  ################################################
   *  SYMBOLIC FACTORIZATION (last_umfpack_step == 1):
   *
   *  UMFPACK_OK
   *      Each column of the input matrix contained row indices
   *      in increasing order, with no duplicates. Only in this case
   *      does umfpack_*_symbolic compute a valid symbolic factorization.
   *      For the other cases below, no Symbolic object is created
   *      (*Symbolic is (void *) NULL).
   *  UMFPACK_ERROR_n_nonpositive
   *      n is less than or equal to zero.
   *  UMFPACK_ERROR_invalid_matrix
   *      Number of entries in the matrix is negative, Ap [0] is nonzero,
   *      a column has a negative number of entries, a row index is out of
   *      bounds, or the columns of input matrix were jumbled (unsorted
   *      columns or duplicate entries).
   *  UMFPACK_ERROR_out_of_memory
   *      Insufficient memory to perform the symbolic analysis. If the
   *      analysis requires more than 2GB of memory and you are using
   *      the 32-bit ("int") version of UMFPACK, then you are guaranteed
   *      to run out of memory. Try using the 64-bit version of UMFPACK.
   *  UMFPACK_ERROR_argument_missing
   *      One or more required arguments is missing.
   *  UMFPACK_ERROR_internal_error
   *      Something very serious went wrong. This is a bug.
   *      Please contact the author (DrTimothyAldenDavis@gmail.com).
   *
   *  ################################################
   *  NUMERIC FACTORIZATION (last_umfpack_step == 2):
   *
   *  UMFPACK_OK
   *      Numeric factorization was successful. umfpack_di_numeric
   *      computed a valid numeric factorization.
   *      UMFPACK_WARNING_singular_matrix
   *      Numeric factorization was successful, but the matrix is
   *      singular. umfpack_di_numeric computed a valid numeric
   *      factorization, but you will get a divide by zero in
   *      umfpack_di_solve. For the other cases below, no Numeric object
   *      is created (*Numeric is (void *) NULL).
   *  UMFPACK_ERROR_out_of_memory
   *      Insufficient memory to complete the numeric factorization.
   *  UMFPACK_ERROR_argument_missing
   *      One or more required arguments are missing.
   *  UMFPACK_ERROR_invalid_Symbolic_object
   *      Symbolic object provided as input is invalid.
   *  UMFPACK_ERROR_different_pattern
   *      The pattern (Ap and/or Ai) has changed since the call to
   *      umfpack_di_symbolic which produced the Symbolic object.
   *
   *  ################################################
   *  SOLVE STEP (last_umfpack_step == 3):
   * 
   *  UMFPACK_OK
   *      The linear system was successfully solved.
   *  UMFPACK_WARNING_singular_matrix
   *       A divide-by-zero occurred. Your solution will contain Inf’s
   *       and/or NaN’s. Some parts of the solution may be valid. For
   *       example, solving Ax=b with
   *       A = [2 0] b = [ 1 ] returns x = [ 0.5 ]
   *           [0 0]     [ 0 ]             [ Inf ]
   *  UMFPACK_ERROR_out_of_memory
   *       Insufficient memory to solve the linear system.
   *  UMFPACK_ERROR_argument_missing
   *       One or more required arguments are missing. The B, X, (or
   *       Bx and Xx for the complex versions) arguments
   *       are always required. Info and Control are not required. Ap,
   *       Ai, Ax are required if Ax=b,
   *       A’x=b, A.’x=b is to be solved, the (default) iterative
   *       refinement is requested, and the matrix A is nonsingular.
   *  UMFPACK_ERROR_invalid_system
   *       The sys argument is not valid, or the matrix A is not square.
   *  UMFPACK_ERROR_invalid_Numeric_object
   *       The Numeric object is not valid.
   */
  int error;
  /* Last umfpack routine:
   *    0   ||  initial value, no umfpack routines called so far
   *    1   ||  last step was a symbolic factorization
   *    2   ||  last step was a numeric factorization
   *    3   ||  last step was a solve step
   */
  int last_umfpack_step;

  /* Maximal number of factors with identical nonzero sparsity structure that
   * the user would like to keep at the same time in memory. It is possible to
   * store several different factorizations with the same nonzero structure at
   * the same time in the internal data management of the solver. In many
   * applications this value is equal to 1.
   * Note: Matrices with different sparsity structure can be kept in memory
   * with different memory address pointers pt. */
  int maxfct;
  /* Actual matrix for the solution phase. With this scalar the user
   * can define the matrix that he would like to factorize.
   * The value must be: 1 ≤ mnum ≤ maxfct. In many applications this
   * value is equal to 1. */
  int mnum;

  /* Number of equations in the sparse linear systems of equations */
  int n;
  /* ai is an integer array of size n+1. ai(i) points to the first column index
   * of row i in the array ja in compressed sparse row format. */
  int *ia;
  /* Nonzero values of the coefficient matrix a. The array a contains the values
   * corresponding to the indices in ja. The size and order of a is the same as
   * that of ja. The coefficients can be either real or complex. The matrix must
   * be stored in compressed sparse row format with increasing values of ja for
   * each row. */
  double *a;
  /* The integer array ja contains the column indices of the sparse matrix A
   * stored in compress sparse row format. The indices in each row must be
   * sorted in increasing order. */
  int *ja;
  /* nrhs is the number of right-hand sides that need to be solved for. */
  int nrhs;

};
typedef struct _SUNLinearSolverContent_UMFPACK *SUNLinearSolverContent_UMFPACK;

/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_umfpack
 *
 * CONSTRUCTOR:
 *    SUNUMFPACK creates and allocates memory for a UMFPACK sparse-direct
 *      linear solver
 *
 * OTHER:
 *        NOTE: NOT UP-TO-DATE, NOT SURE WHETHER REQUIRED
 *    SUNUMFPACKReInit reinitializes memory and flags for a new
 *      factorization (symbolic and numeric) to be conducted at the
 *      next solver setup call.  This routine is useful in the
 *      cases where the number of nonzeroes has changed or if the
 *      structure of the linear system has changed which would
 *      require a new symbolic (and numeric factorization).
 *
 *      The reinit_type argument governs the level of
 *      reinitialization:
 *
 *      reinit_type = 1: The Jacobian matrix will be destroyed and
 *                       a new one will be allocated based on the
 *                       nnz value passed to this call. New
 *                       symbolic and numeric factorizations will
 *                       be completed at the next solver setup.
 *
 *      reinit_type = 2: Only symbolic and numeric factorizations
 *                       will be completed.  It is assumed that the
 *                       Jacobian size has not exceeded the size of
 *                       nnz given in the sparse matrix provided to
 *                       the original constructor routine (or the
 *                       previous SUNUMFPACKReInit call)
 *
 *      This routine assumes no other changes to solver use are
 *      necessary.
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_UMFPACK(N_Vector y, SUNMatrix A);

SUNDIALS_EXPORT int SUNLinSol_UMFPACKReInit(SUNLinearSolver S, SUNMatrix A,
                                     sunindextype nnz, int reinit_type);

SUNDIALS_EXPORT int SUNLinSol_UMFPACKSetOrdering(SUNLinearSolver S,
                                          int ordering_choice);

/*
 * -----------------------------------------------------------------
 * UMFPACK implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type
SUNLinSolGetType_UMFPACK(SUNLinearSolver S);
SUNDIALS_EXPORT SUNLinearSolver_ID SUNLinSolGetID_UMFPACK(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_UMFPACK(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_UMFPACK(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_UMFPACK(SUNLinearSolver S, SUNMatrix A,
                                           N_Vector x, N_Vector b,
                                           realtype tol);
SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_UMFPACK(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_UMFPACK(SUNLinearSolver S, long int *lenrwLS,
                                           long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_UMFPACK(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif  // _SUNLINSOL_UMFPACK_H
