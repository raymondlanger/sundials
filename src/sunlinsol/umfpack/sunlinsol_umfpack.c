/*
 * -----------------------------------------------------------------
 * Programmer:  Raymond Langer @ RWTH Aachen University
 * -----------------------------------------------------------------
 * This is the implementation file for the UMFPACK implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_umfpack.h>

#define ONE RCONST(1.0)
#define TWOTHIRDS RCONST(0.666666666666666666666666666666667)

/*
 * -----------------------------------------------------------------
 * UMFPACK solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define UMFPACK_CONTENT(S) ((SUNLinearSolverContent_UMFPACK)(S->content))
#define CONTROL(S) (UMFPACK_CONTENT(S)->Control)
#define INFO(S) (UMFPACK_CONTENT(S)->Info)

#define SYMBOLIC(S) (UMFPACK_CONTENT(S)->Symbolic)
#define NUMERIC(S) (UMFPACK_CONTENT(S)->Numeric)

#define ERROR(S) (UMFPACK_CONTENT(S)->error)
#define LAST_UMFPACK_STEP(S) (UMFPACK_CONTENT(S)->last_umfpack_step)

#define MAXFCT(S) (UMFPACK_CONTENT(S)->maxfct)
#define MNUN(S) (UMFPACK_CONTENT(S)->mnum)
#define N(S) (UMFPACK_CONTENT(S)->n)
#define IA(S) (UMFPACK_CONTENT(S)->ia)
#define A(S) (UMFPACK_CONTENT(S)->a)
#define JA(S) (UMFPACK_CONTENT(S)->ja)
#define NRHS(S) (UMFPACK_CONTENT(S)->nrhs)
#define B(S) (UMFPACK_CONTENT(S)->b)
#define X(S) (UMFPACK_CONTENT(S)->x)

#define LASTFLAG(S) (UMFPACK_CONTENT(S)->last_flag)
#define FIRSTFACTORIZE(S) (UMFPACK_CONTENT(S)->first_factorize)
/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new UMFPACK linear solver
 */

SUNLinearSolver SUNLinSol_UMFPACK(N_Vector y, SUNMatrix A) {
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_UMFPACK content;
  int flag;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) {
    printf("Matrix type is inappropriate for UMFPACK\n");
    return (NULL);
  }
  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) {
    return (NULL);
  }
  if ((N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS))
    return (NULL);

  if (SUNSparseMatrix_Rows(A) != N_VGetLength(y))
    return (NULL);

  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver)malloc(sizeof *S);
  if (S == NULL)
    return (NULL);

  /* Create linear solver operation structure */
  ops = NULL;
  ops =
      (SUNLinearSolver_Ops)malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) {
    free(S);
    return (NULL);
  }
  /* Attach operations */
  ops->gettype = SUNLinSolGetType_UMFPACK;
  ops->getid = SUNLinSolGetID_UMFPACK;
  ops->initialize = SUNLinSolInitialize_UMFPACK;
  ops->setup = SUNLinSolSetup_UMFPACK;
  ops->solve = SUNLinSolSolve_UMFPACK;
  ops->lastflag = SUNLinSolLastFlag_UMFPACK;
  ops->space = SUNLinSolSpace_UMFPACK;
  ops->free = SUNLinSolFree_UMFPACK;
  ops->setatimes = NULL;
  ops->setpreconditioner = NULL;
  ops->setscalingvectors = NULL;
  ops->numiters = NULL;
  ops->resnorm = NULL;
  ops->resid = NULL;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_UMFPACK)malloc(
      sizeof(struct _SUNLinearSolverContent_UMFPACK));
  if (content == NULL) {
    free(ops);
    free(S);
    return (NULL);
  }

  S->content = content;
  /* Fill content */
  LASTFLAG(S) = 0;
  FIRSTFACTORIZE(S) = 1;
  LAST_UMFPACK_STEP(S) = 0;

  S->ops = ops;
#if defined(SUNDIALS_INT64_T)
#error Incompatible sunindextype. Currently (08/01/2020);
#elif defined(SUNDIALS_INT32_T)
#else
#error sunindextype not implemented for UMFPACK
#endif

  ERROR(S) = UMFPACK_OK;

  umfpack_di_defaults(CONTROL(S));
  int i = 0;
  for (i = 0; i < UMFPACK_INFO; ++i)
    INFO(S)[i] = 0.0;
  //}

  MAXFCT(S) = 1; /* Maximum number of numerical factorizations.  */
  MNUN(S) = 1;   /* Which factorization to use. */
  NRHS(S) = 1;

  SYMBOLIC(S) = NULL;
  NUMERIC(S) = NULL;

  N(S) = SUNSparseMatrix_NP(A);
  IA(S) = SUNSparseMatrix_IndexPointers(A);
  JA(S) = SUNSparseMatrix_IndexValues(A);
  A(S) = SUNSparseMatrix_Data(A);

  return (S);
}

/*
/----------------------------------------------------------------------------
* Function to reinitialize a UMFPACK linear solver
*/

int SUNLinSol_UMFPACKReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz,
                            int reinit_type) {
  sunindextype n;
  int type;

  /* Check for non-NULL SUNLinearSolver */
  if ((S == NULL) || (A == NULL))
    return (SUNLS_MEM_NULL);

  /* Check for valid SUNMatrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
    return (SUNLS_ILL_INPUT);

  /* Check for valid reinit_type */
  if ((reinit_type != 1) && (reinit_type != 2))
    return (SUNLS_ILL_INPUT);

  /* Perform re-initialization */
  if (reinit_type == 1) {
    /* Get size/type of current matrix */
    n = SUNSparseMatrix_Rows(A);
    type = SUNSparseMatrix_SparseType(A);

    /* Destroy previous matrix */
    SUNMatDestroy(A);

    /* Create new sparse matrix */
    A = SUNSparseMatrix(n, n, nnz, type);
    if (A == NULL)
      return (SUNLS_MEM_FAIL);
  }

  printf("REVISE REINIT METHOD\n");
  exit(2);
  /* Free the prior factorazation and reset for first factorization */
  if (SYMBOLIC(S) != NULL) {
    umfpack_di_free_symbolic(&SYMBOLIC(S));
    SYMBOLIC(S) = NULL;
  }
  if (NUMERIC(S) != NULL) {
    umfpack_di_free_numeric(&NUMERIC(S));
    NUMERIC(S) = NULL;
  }
  FIRSTFACTORIZE(S) = 1;
  N(S) = SUNSparseMatrix_NP(A);
  IA(S) = SUNSparseMatrix_IndexPointers(A);
  JA(S) = SUNSparseMatrix_IndexValues(A);
  A(S) = SUNSparseMatrix_Data(A);

  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a UMFPACK linear solver
 */

// int SUNLinSol_UMFPACKSetOrdering(SUNLinearSolver S, int ordering_choice)
//{
//  /* Check for legal ordering_choice */
//  if ((ordering_choice < 0) || (ordering_choice > 2))
//    return(SUNLS_ILL_INPUT);
//
//  /* Check for non-NULL SUNLinearSolver */
//  if (S == NULL) return(SUNLS_MEM_NULL);
//
//  /* Set ordering_choice */
//  COMMON(S).ordering = ordering_choice;
//
//  LASTFLAG(S) = SUNLS_SUCCESS;
//  return(LASTFLAG(S));
//}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_UMFPACK(SUNLinearSolver S) {
  return (SUNLINEARSOLVER_DIRECT);
}

SUNLinearSolver_ID SUNLinSolGetID_UMFPACK(SUNLinearSolver S) {
  return (SUNLINEARSOLVER_UMFPACK);
}

int SUNLinSolInitialize_UMFPACK(SUNLinearSolver S) {
  /* Force symbolic factorization */
  FIRSTFACTORIZE(S) = 1;

  LASTFLAG(S) = SUNLS_SUCCESS;
  LAST_UMFPACK_STEP(S) = 0;
  return (LASTFLAG(S));
}

int SUNLinSolSetup_UMFPACK(SUNLinearSolver S, SUNMatrix A) {
  int retval;
  realtype uround_twothirds = SUNRpowerR(UNIT_ROUNDOFF, TWOTHIRDS);

  /* Ensure that A is a sparse matrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return (LASTFLAG(S));
  }
  //if (N(S) != SUNSparseMatrix_NP(A)) {
  //  printf(
  //      "WARNING Dimensions of A are different from the stored dimensions\n");
  //}
  //if (IA(S) != SUNSparseMatrix_IndexPointers(A)) {
  //  printf("WARNING Index pointer of A is different from the stored pointer\n");
  //}
  //if (JA(S) != SUNSparseMatrix_IndexValues(A)) {
  //  printf("WARNING Index values pointer of A is different from the stored "
  //         "pointer\n");
  //}
  //if (A(S) != SUNSparseMatrix_Data(A)) {
  //  printf("WARNING Data pointer of A is different from the stored pointer\n");
  //}

  /* On first decomposition, get the symbolic factorization */
  if (FIRSTFACTORIZE(S)) {
    if (SYMBOLIC(S)) {
      umfpack_di_free_symbolic(&SYMBOLIC(S));
      SYMBOLIC(S) = NULL;
    }

    ERROR(S) = umfpack_di_symbolic(N(S), N(S), IA(S), JA(S), A(S), &SYMBOLIC(S),
                                   CONTROL(S), INFO(S));
    if (ERROR(S) != 0) {
      printf("\nERROR during symbolic and numerical factorization: %d",
             ERROR(S));
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return (LASTFLAG(S));
    }
    // printf("\nReordering completed ... ");
    // printf("\nNumber of nonzeros in factors  = %d", CONTROL(S)[17]);
    // printf("\nNumber of factorization MFLOPS = %d", CONTROL(S)[18]);

    FIRSTFACTORIZE(S) = 0;
    LAST_UMFPACK_STEP(S) = 1;
  }

  ERROR(S) = umfpack_di_numeric(IA(S), JA(S), A(S), SYMBOLIC(S), &NUMERIC(S),
                                CONTROL(S), INFO(S));

  if (ERROR(S) != 0) {
    printf("\nERROR during numerical factorization: %d", ERROR(S));
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    return (LASTFLAG(S));
    // printf("\nNumerical factorization completed ...\n ");
    //
    //   /*-----------------------------------------------------------
    //     Check if a cheap estimate of the reciprocal of the condition
    //     number is getting too small.  If so, delete
    //     the prior numeric factorization and recompute it.
    //     -----------------------------------------------------------*/
    //
    //   retval = sun_umfpack_rcond(SYMBOLIC(S), NUMERIC(S), &COMMON(S));
    //   if (retval == 0) {
    //     LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    //     return(LASTFLAG(S));
    //   }
    //
    //   if ( COMMON(S).rcond < uround_twothirds ) {
    //
    //     /* Condition number may be getting large.
    //	 Compute more accurate estimate */
    //     retval = sun_umfpack_condest(SUNSparseMatrix_IndexPointers(A),
    //                              SUNSparseMatrix_Data(A),
    //                              SYMBOLIC(S),
    //                              NUMERIC(S),
    //                              &COMMON(S));
    // if (retval == 0) {
    // LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    // return(LASTFLAG(S));
    // }
    //
    // if ( COMMON(S).condest > (ONE/uround_twothirds) ) {
    //
    // /* More accurate estimate also says condition number is
    // large, so recompute the numeric factorization */
    // sun_umfpack_free_numeric(&NUMERIC(S), &COMMON(S));
    // NUMERIC(S) = sun_umfpack_factor(SUNSparseMatrix_IndexPointers(A),
    // SUNSparseMatrix_IndexValues(A),
    // SUNSparseMatrix_Data(A),
    // SYMBOLIC(S),
    // &COMMON(S));
    // if (NUMERIC(S) == NULL) {
    // LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
    // return(LASTFLAG(S));
    // }
    // }
    //
    // }
  }

  LAST_UMFPACK_STEP(S) = 2;
  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

int SUNLinSolSolve_UMFPACK(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                           N_Vector b, realtype tol) {
  int flag;
  realtype *xdata;
  /* check for valid inputs */
  if ((A == NULL) || (S == NULL) || (x == NULL) || (b == NULL))
    return (SUNLS_MEM_NULL);

  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return (LASTFLAG(S));
  }

  //if (N(S) != SUNSparseMatrix_NP(A)) {
  //  printf(
  //      "WARNING Dimensions of A are different from the stored dimensions\n");
  //}
  //if (IA(S) != SUNSparseMatrix_IndexPointers(A)) {
  //  printf("WARNING Index pointer of A is different from the stored pointer\n");
  //}
  //if (JA(S) != SUNSparseMatrix_IndexValues(A)) {
  //  printf("WARNING Index values pointer of A is different from the stored "
  //         "pointer\n");
  //}
  //if (A(S) != SUNSparseMatrix_Data(A)) {
  //  printf("WARNING Data pointer of A is different from the stored pointer\n");
  //}

  /* Call UMFPACK to solve the linear system */
  /* -------------------------------------------------------------------- */
  /* ..  Back substitution and iterative refinement.                      */
  /* -------------------------------------------------------------------- */
  CONTROL(S)
  [UMFPACK_IRSTEP] = 1; /* Max numbers of iterative refinement steps. */
  if (SUNSparseMatrix_SparseType(A) == CSC_MAT) {
    // UMFPACK expect CSR format, solve transposed system J^T x = b instead
    // printf("[UMFPACK]:  solve transposed system J^T x = b (CSC_MAT used)\n");
    umfpack_di_solve(UMFPACK_A, IA(S), JA(S), A(S), xdata,
                     N_VGetArrayPointer(b), NUMERIC(S), CONTROL(S), INFO(S));
  } else {
    // printf("[UMFPACK]:  solve regular system J x = b (CSR_MAT used)\n");

    umfpack_di_solve(UMFPACK_At, IA(S), JA(S), A(S), xdata,
                     N_VGetArrayPointer(b), NUMERIC(S), CONTROL(S), INFO(S));
  }

  ERROR(S) = INFO(S)[UMFPACK_STATUS];
  if (ERROR(S) != UMFPACK_OK) {
    printf("\nERROR during solve step: %d", ERROR(S));
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    return (LASTFLAG(S));
  }

  LAST_UMFPACK_STEP(S) = 3;
  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

sunindextype SUNLinSolLastFlag_UMFPACK(SUNLinearSolver S) {
  /* return the stored 'last_flag' value */
  return (LASTFLAG(S));
}

int SUNLinSolSpace_UMFPACK(SUNLinearSolver S, long int *lenrwLS,
                           long int *leniwLS) {
  // code below is just copied from klu and might not be correct (rlanger
  // 08/01/2020) 
  // printf("%s %s %d:%s()\n", "Fix this:\n", __FILE__, __LINE__,
  // __func__); 
  // exit(99);

  /* since the umfpack structures are opaque objects, we
     omit those from these results */
  *leniwLS = 2;
  *lenrwLS = 0;
  return (SUNLS_SUCCESS);
}

int SUNLinSolFree_UMFPACK(SUNLinearSolver S) {
  /* return with success if already freed */
  if (S == NULL)
    return (SUNLS_SUCCESS);

  /* delete items from the contents structure (if it exists) */
  if (S->content) {
    if (SYMBOLIC(S) != NULL) {
      umfpack_di_free_symbolic(&SYMBOLIC(S));
      SYMBOLIC(S) = NULL;
    }
    if (NUMERIC(S) != NULL) {
      umfpack_di_free_numeric(&NUMERIC(S));
      NUMERIC(S) = NULL;
    }
    free(S->content);
    S->content = NULL;
  }

  /* delete generic structures */
  if (S->ops) {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return (SUNLS_SUCCESS);
}
