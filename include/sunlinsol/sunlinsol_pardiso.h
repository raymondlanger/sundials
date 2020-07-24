/*
 * -----------------------------------------------------------------
 * Programmer:  Raymond Langer @ RWTH Aachen University
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 * This is the header file for the PARDISO implementation of the
 * SUNLINSOL module.
 *
 * Notes:
 *
 *   - The definition of the generic SUNLinearSolver structure can
 *     be found in the header file sundials_linearsolver.h.
 *
 * -----------------------------------------------------------------
 */

#ifndef _SUNLINSOL_PARDISO_H
#define _SUNLINSOL_PARDISO_H

#include <sundials/sundials_config.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>

#ifdef __cplusplus /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifdef MKL_PARDISO
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

#define PARDISOINIT PARDISO
#define PARDISO PARDISO

#else

#define PARDISOINIT pardisoinit
#define PARDISO pardiso
/* pardisoinit() checks the current license in the file pardiso.lic and
 * initializes the internal timer and the address pointer pt. It sets the
 * solver default values according to the matrix type. */
void pardisoinit(void * /* pt */, int * /* mtype */, int * /* solver */,
                 int * /* iparm */, double * /* dparm */, int * /* error */);

void pardiso(void * /* pt */, int * /* maxfct */, int * /* mnum */,
             int * /* mtype */, int * /* phase */, int * /* n */,
             double * /* a */, int * /* ia */, int * /* ja */, int * /* perm */,
             int * /* nrhs */, int * /* iparm */, int * /* msglvl */,
             double * /* b */, double * /* x */, int * /* error */,
             double * /* dparm */);

void pardiso_chkmatrix(int * /* mtype */, int * /* n */, double * /* a */,
                       int * /* ia */, int * /* ja */, int * /* error */);
void pardiso_chkvec(int * /* n */, int * /* nrhs */, double * /* b */,
                    int * /* error */);
void pardiso_printstats(int *, int *, double *, int *, int *, int *, double *,
                        int *);
#endif // MKL_PARDISO

/* Interfaces to match 'sunindextype' with the correct PARDISO types/functions
 */
#if defined(SUNDIALS_INT64_T)
#elif defined(SUNDIALS_INT32_T)
#else /* incompatible sunindextype for PARDISO */
#error Incompatible sunindextype for PARDISO
#endif

#if defined(SUNDIALS_DOUBLE_PRECISION)
#else
#error Incompatible realtype for PARDISO
#endif

#ifdef MKL_PARDISO
struct _SUNLinearSolverContent_PARDISO {
  // for documentation cf. below
  long int last_flag;
  int first_factorize;
  void *pt[64];
  MKL_INT mtype;
  int solver;
  MKL_INT iparm[64];
  double dparm[64];
  MKL_INT error;
  MKL_INT maxfct;
  MKL_INT mnum;
  MKL_INT phase;
  MKL_INT n;
  MKL_INT *ia;
  double *a;
  MKL_INT *ja;
  MKL_INT *perm;
  MKL_INT nrhs;
  MKL_INT msglvl;
  double *b;
  double *x;
};
#else
struct _SUNLinearSolverContent_PARDISO {
  long int last_flag;
  int first_factorize;
  /* Internal solver memory pointer pt,                  */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
  /* or void *pt[64] should be OK on both architectures  */
  void *pt[64];

  /* mtype is the matrix type; the following matrix types are supported by
   * pardiso */
  /*   1  ||  real and structurally symmetric               */
  /*   2  ||  real and symmetric positive definite          */
  /*  -2  ||  real and symmetric indefinite                 */
  /*   3  ||  complex and structurally symmetric            */
  /*   4  ||  complex and Hermitian positive definite       */
  /*  -4  ||  complex and Hermitian indefinite              */
  /*   6  ||  complex and symmetric                         */
  /*  11  ||  real and nonsymmetric                         */
  /*  13  ||  complex and nonsymmetric                      */
  int mtype; /* real unsymmetric matrix */

  /* Pardiso control parameters. */
  /* the solver method
     0   || sparse direct solver
     1   || multi-recursive iterative solver */
  int solver;
  /*  iparm is an integer array of size 64 that is used to pass various
   *  parameters to PARDISO and to return some useful information after
   *  the execution of the solver */
  int iparm[64];
  /*  dparm is a double-precision array of size 64 that is used to pass various
   *  parameters to PARDISO and to return some useful information after the
   *  execution of the solver. The array will only be used to initialize the
   *  multi-recursive iterative linear solver in PARDISO. */
  double dparm[64];

  /* The error indicator:
   *    0   ||  No error.
   *   -1   ||  Input inconsistent.
   *   -2   ||  Not enough memory.
   *   -3   ||  Reordering problem.
   *   -4   ||  Zero pivot, numerical fact. or iterative refinement problem.
   *   -5   ||  Unclassified (internal) error.
   *   -6   ||  Preordering failed (matrix types 11, 13 only).
   *   -7   ||  Diagonal matrix problem.
   *   -8   ||  32-bit integer overflow problem.
   *  -10   ||  No license file pardiso.lic found.
   *  -11   ||  License is expired.
   *  -12   ||  Wrong username or hostname.
   * -100   ||  Reached maximum number of Krylov-subspace iteration in iterative
   *            solver. 
   * -101   ||  No sufficient convergence in Krylov-subspace iteration
   *            within 25 iterations. 
   * -102   ||  Error in Krylov-subspace iteration. 
   * -103   ||  Break-Down in Krylov-subspace iteration.
   */
  int error;

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

  /* Solving a linear system is split into four tasks:
   * 1. analysis and symbolic factorization
   * 2. numerical factorization
   * 3. forward and backward substitution including iterative refinement
   * 4. termination to release all internal solver memory. phase ≤ 0
   *
   * phase controls the execution of the solver. It is a two-digit integer ij
   * (10i + j, 1 ≤ i ≤ 3, i ≤ j ≤ 3 for normal execution modes). The i digit
   * indicates the starting phase of execution, and j indicates the ending
   * phase. */
  int phase;
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
  /* The user can supply his own fill-in reducing ordering to the solver. This
   * permutation vector is only accessed if iparm[4] = 1.
   * The permutation iperm must be a vector of size n.
   */
  int *perm;
  /* nrhs is the number of right-hand sides that need to be solved for. */
  int nrhs;

  /* Message level information. If msglvl = 0 then PARDISO generates no output.
   * If msglvl = 1 the direct solver prints statistical information to the
   * screen. The multi-recursive iterative solver prints statistical information
   * in the file pardiso-ml.out. */
  int msglvl;

  /* Right-hand side vector/matrix.
   * The array is replaced by the solution if iparm[5] = 1. */
  double *b;
  /* Solution if iparm[5]= 0.
   * Note: x is only accessed in the solution phase. */
  double *x;
};
#endif // MKL_PARDISO
typedef struct _SUNLinearSolverContent_PARDISO *SUNLinearSolverContent_PARDISO;

/*
 * -----------------------------------------------------------------
 * PART II: functions exported by sunlinsol_klu
 *
 * CONSTRUCTOR:
 *    SUNPARDISO creates and allocates memory for a PARDISO sparse-direct
 *      linear solver
 *
 * OTHER:
 *        NOTE: NOT UP-TO-DATE, NOT SURE WHETHER REQUIRED
 *    SUNPARDISOReInit reinitializes memory and flags for a new
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
 *                       previous SUNPARDISOReInit call)
 *
 *      This routine assumes no other changes to solver use are
 *      necessary.
 *
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver SUNPARDISO(N_Vector y, SUNMatrix A);

SUNDIALS_EXPORT int SUNPARDISOReInit(SUNLinearSolver S, SUNMatrix A,
                                     sunindextype nnz, int reinit_type);

SUNDIALS_EXPORT int SUNPARDISOSetOrdering(SUNLinearSolver S,
                                          int ordering_choice);

/*
 * -----------------------------------------------------------------
 * PARDISO implementations of various useful linear solver operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT SUNLinearSolver_Type
SUNLinSolGetType_PARDISO(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolInitialize_PARDISO(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSetup_PARDISO(SUNLinearSolver S, SUNMatrix A);
SUNDIALS_EXPORT int SUNLinSolSolve_PARDISO(SUNLinearSolver S, SUNMatrix A,
                                           N_Vector x, N_Vector b,
                                           realtype tol);
SUNDIALS_EXPORT long int SUNLinSolLastFlag_PARDISO(SUNLinearSolver S);
SUNDIALS_EXPORT int SUNLinSolSpace_PARDISO(SUNLinearSolver S, long int *lenrwLS,
                                           long int *leniwLS);
SUNDIALS_EXPORT int SUNLinSolFree_PARDISO(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif
