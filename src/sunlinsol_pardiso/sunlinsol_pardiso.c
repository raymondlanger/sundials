/*
 * -----------------------------------------------------------------
 * Programmer:  Raymond Langer @ RWTH Aachen University
 * -----------------------------------------------------------------
 * This is the implementation file for the PARDISO implementation of
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_pardiso.h>

#define ZERO RCONST(0.0)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)
#define TWOTHIRDS RCONST(0.666666666666666666666666666666667)

/* Private function prototypes */
static sunindextype GlobalVectorLength_PARDISO(N_Vector y);

/*
 * -----------------------------------------------------------------
 * PARDISO solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define PARDISO_CONTENT(S) ((SUNLinearSolverContent_PARDISO)(S->content))
#define PT(S) (PARDISO_CONTENT(S)->pt)
#define MTYPE(S) (PARDISO_CONTENT(S)->mtype)
#define SOLVER(S) (PARDISO_CONTENT(S)->solver)
#define IPARM(S) (PARDISO_CONTENT(S)->iparm)
#define DPARM(S) (PARDISO_CONTENT(S)->dparm)
#define ERROR(S) (PARDISO_CONTENT(S)->error)

#define MAXFCT(S) (PARDISO_CONTENT(S)->maxfct)
#define MNUN(S) (PARDISO_CONTENT(S)->mnum)
#define PHASE(S) (PARDISO_CONTENT(S)->phase)
#define N(S) (PARDISO_CONTENT(S)->n)
#define IA(S) (PARDISO_CONTENT(S)->ia)
#define A(S) (PARDISO_CONTENT(S)->a)
#define JA(S) (PARDISO_CONTENT(S)->ja)
#define PERM(S) (PARDISO_CONTENT(S)->perm)
#define NRHS(S) (PARDISO_CONTENT(S)->nrhs)
#define MSGLVL(S) (PARDISO_CONTENT(S)->msglvl)
#define B(S) (PARDISO_CONTENT(S)->b)
#define X(S) (PARDISO_CONTENT(S)->x)

#define LASTFLAG(S) (PARDISO_CONTENT(S)->last_flag)
#define FIRSTFACTORIZE(S) (PARDISO_CONTENT(S)->first_factorize)
/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new PARDISO linear solver
 */

SUNLinearSolver SUNPARDISO(N_Vector y, SUNMatrix A) {
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_PARDISO content;
  sunindextype MatrixRows, VecLength;
  int flag;
  int i = 0;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) {
    printf("Matrix type is inappropriate for PARDISO\n");
    return (NULL);
  }
  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) {
    return (NULL);
  }
  MatrixRows = SUNSparseMatrix_Rows(A);
  if ((N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
      (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS))
    return (NULL);

  /* optimally this function would be replaced with a generic N_Vector routine
   */
  VecLength = GlobalVectorLength_PARDISO(y);
  if (MatrixRows != VecLength) return (NULL);

  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver)malloc(sizeof *S);
  if (S == NULL) return (NULL);

  /* Create linear solver operation structure */
  ops = NULL;
  ops =
      (SUNLinearSolver_Ops)malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) {
    free(S);
    return (NULL);
  }
  /* Attach operations */
  ops->gettype = SUNLinSolGetType_PARDISO;
  ops->initialize = SUNLinSolInitialize_PARDISO;
  ops->setup = SUNLinSolSetup_PARDISO;
  ops->solve = SUNLinSolSolve_PARDISO;
  ops->lastflag = SUNLinSolLastFlag_PARDISO;
  ops->space = SUNLinSolSpace_PARDISO;
  ops->free = SUNLinSolFree_PARDISO;
  ops->setatimes = NULL;
  ops->setpreconditioner = NULL;
  ops->setscalingvectors = NULL;
  ops->numiters = NULL;
  ops->resnorm = NULL;
  ops->resid = NULL;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_PARDISO)malloc(
      sizeof(struct _SUNLinearSolverContent_PARDISO));
  if (content == NULL) {
    free(ops);
    free(S);
    return (NULL);
  }
  /* Fill content */
  content->last_flag = 0;
  content->first_factorize = 1;

  S->content = content;
  S->ops = ops;

  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  for (i = 0; i < 64; ++i) IPARM(S)[i] = 0;
  IPARM(S)[0] = 1;  /* No solver default */
  IPARM(S)[1] = 2;  /* Fill-in reordering from METIS */
  IPARM(S)[3] = 0;  /* No iterative-direct algorithm */
  IPARM(S)[4] = 0;  /* No user fill-in reducing permutation */
  IPARM(S)[5] = 0;  /* Write solution into x */
  IPARM(S)[6] = 0;  /* Not in use */
  IPARM(S)[7] = 2;  /* Max numbers of iterative refinement steps */
  IPARM(S)[8] = 0;  /* Not in use */
  IPARM(S)[9] = 13; /* Perturb the pivot elements with 1E-13 */
  IPARM(S)[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
  // the distinction between CSC_MAT format and CSR_MAT format is done in 
  IPARM(S)[11] = 0;
  IPARM(S)
  [12] = 1; /* Maximum weighted matching algorithm is switched-on (default for
               non-symmetric) */
  IPARM(S)[13] = 0;  /* Output: Number of perturbed pivots */
  IPARM(S)[14] = 0;  /* Not in use */
  IPARM(S)[15] = 0;  /* Not in use */
  IPARM(S)[16] = 0;  /* Not in use */
  IPARM(S)[17] = -1; /* Output: Number of nonzeros in the factor LU */
  IPARM(S)[18] = -1; /* Output: Mflops for LU factorization */
  IPARM(S)[19] = 0;  /* Output: Numbers of CG Iterations */
  IPARM(S)
  [23] = 0; /* 1 selects a two-level scheduling algorithms which should result
               in a better parallel efficiency */

#if defined(SUNDIALS_INT64_T)
#error Incompatible sunindextype. Currently (03/25/2018); not sure if I only have to implement this or if this is not supported by PARDISO
#elif defined(SUNDIALS_INT32_T)
#else
#error Incompatible sunindextype for PARDISO
#endif

  SOLVER(S) = 0; /* use sparse direct solver */
  ERROR(S) = 0;
  MTYPE(S) = 11;

#ifdef MKL_PARDISO
  for (i = 0; i < 64; ++i)
        PT(S)[i] = 0;
  //PARDISOINIT(PT(S), &MTYPE(S), IPARM(S));
#else
  PARDISOINIT(PT(S), &MTYPE(S), &SOLVER(S), IPARM(S), DPARM(S), &ERROR(S));

  if (ERROR(S) != 0) {
    if (ERROR(S) == -10) printf("No license file found \n");
    if (ERROR(S) == -11) printf("License is expired \n");
    if (ERROR(S) == -12) printf("Wrong username or hostname \n");
    exit(ERROR(S));
  } else
    printf("[PARDISO]: License check was successful ... \n");
#endif
  /* Numbers of processors, value of OMP_NUM_THREADS */
  char *var = getenv("OMP_NUM_THREADS");
  if (var != NULL) {
    sscanf(var, "%d", &IPARM(S)[2]);
    printf("[PARDISO]: Using %d threads\n", IPARM(S)[2]);
  } else {
    IPARM(S)[2] = 1;
    printf(
        "[PARDISO]: Environment variable OMP_NUM_THREADS not set, using 1\n");
  }

  MAXFCT(S) = 1; /* Maximum number of numerical factorizations.  */
  MNUN(S) = 1;   /* Which factorization to use. */
  MSGLVL(S) = 0; /* Print statistical information  */
  NRHS(S) = 1;

  N(S) = SUNSparseMatrix_NP(A);
  IA(S) = SUNSparseMatrix_IndexPointers(A);
  JA(S) = SUNSparseMatrix_IndexValues(A);
  A(S) = SUNSparseMatrix_Data(A);

  return (S);
}

/*
/----------------------------------------------------------------------------
* Function to reinitialize a PARDISO linear solver
*/

int SUNPARDISOReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz,
                     int reinit_type) {
  sunindextype n;
  int type;

  /* Check for non-NULL SUNLinearSolver */
  if ((S == NULL) || (A == NULL)) return (SUNLS_MEM_NULL);

  /* Check for valid SUNMatrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return (SUNLS_ILL_INPUT);

  /* Check for valid reinit_type */
  if ((reinit_type != 1) && (reinit_type != 2)) return (SUNLS_ILL_INPUT);

  /* Perform re-initialization */
  if (reinit_type == 1) {
    /* Get size/type of current matrix */
    n = SUNSparseMatrix_Rows(A);
    type = SUNSparseMatrix_SparseType(A);

    /* Destroy previous matrix */
    SUNMatDestroy(A);

    /* Create new sparse matrix */
    A = SUNSparseMatrix(n, n, nnz, type);
    if (A == NULL) return (SUNLS_MEM_FAIL);
  }

  printf("REVISE REINIT METHOD\n");
  exit(2);
  /* Free the prior factorazation and reset for first factorization */
  // if( SYMBOLIC(S) != NULL)
  //  sun_pardiso_free_symbolic(&SYMBOLIC(S), &COMMON(S));
  // if( NUMERIC(S) != NULL)
  //  sun_pardiso_free_numeric(&NUMERIC(S), &COMMON(S));
  FIRSTFACTORIZE(S) = 1;
  N(S) = SUNSparseMatrix_NP(A);
  IA(S) = SUNSparseMatrix_IndexPointers(A);
  JA(S) = SUNSparseMatrix_IndexValues(A);
  A(S) = SUNSparseMatrix_Data(A);

  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a PARDISO linear solver
 */

// int SUNPARDISOSetOrdering(SUNLinearSolver S, int ordering_choice)
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

SUNLinearSolver_Type SUNLinSolGetType_PARDISO(SUNLinearSolver S) {
  return (SUNLINEARSOLVER_DIRECT);
}

int SUNLinSolInitialize_PARDISO(SUNLinearSolver S) {
  /* Force factorization */
  FIRSTFACTORIZE(S) = 1;

  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

int SUNLinSolSetup_PARDISO(SUNLinearSolver S, SUNMatrix A) {
  int retval;
  realtype uround_twothirds = SUNRpowerR(UNIT_ROUNDOFF, TWOTHIRDS);

  /* Ensure that A is a sparse matrix */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) {
    LASTFLAG(S) = SUNLS_ILL_INPUT;
    return (LASTFLAG(S));
  }
  //if(N(S) != SUNSparseMatrix_NP(A)){
  //  printf("WARNING Dimensions of A are different from the stored dimensions\n");
  //}
  //if(IA(S) != SUNSparseMatrix_IndexPointers(A)){
  //  printf("WARNING Index pointer of A is different from the stored pointer\n");
  //}
  //if(JA(S) != SUNSparseMatrix_IndexValues(A)){
  //  printf("WARNING Index values pointer of A is different from the stored pointer\n");
  //}
  //if(A(S) != SUNSparseMatrix_Data(A)){
  //  printf("WARNING Data pointer of A is different from the stored pointer\n");
  //}
  /* -------------------------------------------------------------------- */
  /*  .. pardiso_chk_matrix(...)                                          */
  /*     Checks the consistency of the given matrix.                      */
  /*     Use this functionality only for debugging purposes               */
  /* -------------------------------------------------------------------- */
  if (SUNSparseMatrix_IndexValues(A)[0] == 0) {
    int i;
    for (i = 0; i < SM_NNZ_S(A); ++i) {
      ++SUNSparseMatrix_IndexValues(A)[i];
    }
    for (i = 0; i < SUNSparseMatrix_NP(A) + 1; ++i) {
      ++SUNSparseMatrix_IndexPointers(A)[i];
    }
  }
#ifndef MKL_PARDISO
  pardiso_chkmatrix(&MTYPE(S), &N(S), SUNSparseMatrix_Data(A), SUNSparseMatrix_IndexPointers(A), SUNSparseMatrix_IndexValues(A), &ERROR(S));
#endif
  if (ERROR(S) != 0) {
    printf("\nERROR in SUNLinSolSetup_PARDISO in consistency of matrix: %d",
           ERROR(S));
    exit(1);
  }

  /* On first decomposition, get the symbolic factorization */
  if (FIRSTFACTORIZE(S)) {
    /* -------------------------------------------------------------------- */
    /* ..  Reordering, Symbolic Factorization, Numerical factorization. This
     *     step also allocates all memory that is necessary for the 
     *     factorization.                                                   */
    /* -------------------------------------------------------------------- */
    PHASE(S) = 12;

    PARDISO(PT(S), &MAXFCT(S), &MNUN(S), &MTYPE(S), &PHASE(S), &N(S), SUNSparseMatrix_Data(A),
            SUNSparseMatrix_IndexPointers(A), SUNSparseMatrix_IndexValues(A), NULL, &NRHS(S), IPARM(S), &MSGLVL(S), NULL, NULL,
            &ERROR(S)
#ifndef MKL_PARDISO
,DPARM(S)
#endif
);

    if (ERROR(S) != 0) {
      printf("\nERROR during symbolic and numerical factorization: %d",
             ERROR(S));
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
      return (LASTFLAG(S));
    }
    //printf("\nReordering completed ... ");
    //printf("\nNumber of nonzeros in factors  = %d", IPARM(S)[17]);
    //printf("\nNumber of factorization MFLOPS = %d", IPARM(S)[18]);

    FIRSTFACTORIZE(S) = 0;

  } else { /* not the first decomposition, so just refactor */
    PHASE(S) = 22;
    PARDISO(PT(S), &MAXFCT(S), &MNUN(S), &MTYPE(S), &PHASE(S), &N(S), SUNSparseMatrix_Data(A),
            SUNSparseMatrix_IndexPointers(A), SUNSparseMatrix_IndexValues(A), NULL, &NRHS(S), IPARM(S), &MSGLVL(S), NULL, NULL,
            &ERROR(S)
#ifndef MKL_PARDISO
, DPARM(S)
#endif
);

    if (ERROR(S) != 0) {
      printf("\nERROR during numerical factorization: %d", ERROR(S));
      LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
      return (LASTFLAG(S));
    }
    //printf("\nNumerical factorization completed ...\n ");
    //
    //   /*-----------------------------------------------------------
    //     Check if a cheap estimate of the reciprocal of the condition
    //     number is getting too small.  If so, delete
    //     the prior numeric factorization and recompute it.
    //     -----------------------------------------------------------*/
    //
    //   retval = sun_pardiso_rcond(SYMBOLIC(S), NUMERIC(S), &COMMON(S));
    //   if (retval == 0) {
    //     LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    //     return(LASTFLAG(S));
    //   }
    //
    //   if ( COMMON(S).rcond < uround_twothirds ) {
    //
    //     /* Condition number may be getting large.
    //	 Compute more accurate estimate */
    //     retval = sun_pardiso_condest(SUNSparseMatrix_IndexPointers(A),
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
    // sun_pardiso_free_numeric(&NUMERIC(S), &COMMON(S));
    // NUMERIC(S) = sun_pardiso_factor(SUNSparseMatrix_IndexPointers(A),
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

  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

int SUNLinSolSolve_PARDISO(SUNLinearSolver S, SUNMatrix A, N_Vector x,
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

  //if(N(S) != SUNSparseMatrix_NP(A)){
  //  printf("WARNING Dimensions of A are different from the stored dimensions\n");
  //}
  //if(IA(S) != SUNSparseMatrix_IndexPointers(A)){ 
  //  printf("WARNING Index pointer of A is different from the stored pointer\n");
  //}
  //if(JA(S) != SUNSparseMatrix_IndexValues(A)){
  //  printf("WARNING Index values pointer of A is different from the stored pointer\n");
  //}
  //if(A(S) != SUNSparseMatrix_Data(A)){
  //  printf("WARNING Data pointer of A is different from the stored pointer\n");
  //}
  /* Call PARDISO to solve the linear system */
  /* -------------------------------------------------------------------- */
  /* ..  Back substitution and iterative refinement.                      */
  /* -------------------------------------------------------------------- */

  IPARM(S)[7] = 1; /* Max numbers of iterative refinement steps. */
  PHASE(S) = 33;
  if (SUNSparseMatrix_SparseType(A) == CSC_MAT) {
   // PARDISO expect CSR format, solve transposed system J^T x = b instead
   //printf("[PARDISO]:  solve transposed system J^T x = b (CSC_MAT used)\n");
   IPARM(S)[11] = 1;
  } else {
    //printf("[PARDISO]:  solve regular system J x = b (CSR_MAT used)\n");
    IPARM(S)[11] = 0;
  }
  PARDISO(PT(S), &MAXFCT(S), &MNUN(S), &MTYPE(S), &PHASE(S), &N(S), SUNSparseMatrix_Data(A), SUNSparseMatrix_IndexPointers(A),
          SUNSparseMatrix_IndexValues(A), NULL, &NRHS(S), IPARM(S), &MSGLVL(S), N_VGetArrayPointer(b),
          xdata, &ERROR(S)
#ifndef MKL_PARDISO
, DPARM(S)
#endif
);

  if (ERROR(S) != 0) {
    printf("\nERROR during solution: %d", ERROR(S));
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_REC;
    return (LASTFLAG(S));
  }

  LASTFLAG(S) = SUNLS_SUCCESS;
  return (LASTFLAG(S));
}

long int SUNLinSolLastFlag_PARDISO(SUNLinearSolver S) {
  /* return the stored 'last_flag' value */
  return (LASTFLAG(S));
}

int SUNLinSolSpace_PARDISO(SUNLinearSolver S, long int *lenrwLS,
                           long int *leniwLS) {
  // this is might not be up to date (rlagner 26.03.2018)
  /* since the pardiso structures are opaque objects, we
     omit those from these results */
  *leniwLS = 2;
  *lenrwLS = 0;
  return (SUNLS_SUCCESS);
}

int SUNLinSolFree_PARDISO(SUNLinearSolver S) {
  /* return with success if already freed */
  if (S == NULL) return (SUNLS_SUCCESS);

  /* -------------------------------------------------------------------- */
  /* ..  Termination and release of memory.                               */
  /* -------------------------------------------------------------------- */
  if (PHASE(S) == -1) printf("Warning: pardiso seems to be freed twice\n");
  PHASE(S) = -1;
  PARDISO(PT(S), &MAXFCT(S), &MNUN(S), &MTYPE(S), &PHASE(S), &N(S), NULL, IA(S),
          JA(S), NULL, &NRHS(S), IPARM(S), &MSGLVL(S), NULL, NULL, &ERROR(S)
#ifndef MKL_PARDISO
,DPARM(S)
#endif
);

  /* delete items from the contents structure (if it exists) */
  if (S->content) {
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

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

/* Inefficient kludge for determining the number of entries in a N_Vector
  object (replace if such a routine is ever added to the N_Vector API).

  Returns "-1" on an error. */
static sunindextype GlobalVectorLength_PARDISO(N_Vector y) {
  realtype len;
  N_Vector tmp = NULL;
  tmp = N_VClone(y);
  if (tmp == NULL) return (-1);
  N_VConst(ONE, tmp);
  len = N_VDotProd(tmp, tmp);
  N_VDestroy(tmp);
  return ((sunindextype)len);
}
