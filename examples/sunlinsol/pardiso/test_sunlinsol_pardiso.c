/*/*
 * -----------------------------------------------------------------
 * Programmer:  Raymond Langer @ RWTH Aachen University
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol PARDISO module 
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_pardiso.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunlinsol.h"

/* ----------------------------------------------------------------------
 * SUNPARDISO Linear Solver Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int             fails = 0;          /* counter for test failures  */
  sunindextype    N;                  /* matrix columns, rows       */
  SUNLinearSolver LS;                 /* linear solver object       */
  SUNMatrix       A, B;               /* test matrices              */
  N_Vector        x, y, b;            /* test vectors               */
  realtype        *matdata, *xdata;
  int             mattype, print_timing;

  /* check input and set matrix dimensions */
  if (argc < 4){
    printf("ERROR: THREE (3) Inputs required: matrix size, matrix type (0/1), print timing \n");
    return(-1);
  }

  N = atol(argv[1]); 
  if (N <= 0) {
    printf("ERROR: matrix size must be a positive integer \n");
    return(-1); 
  }

  mattype = atoi(argv[2]);
  if ((mattype != 0) && (mattype != 1)) {
    printf("ERROR: matrix type must be 0 or 1 \n");
    return(-1); 
  }
  mattype = (mattype == 0) ? CSC_MAT : CSR_MAT;
  
  print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  printf("\nPARDISO linear solver test: size %ld, type %i\n\n",
         (long int) N, mattype);

  /* Create matrices and vectors */
  B = SUNDenseMatrix(N, N);
  x = N_VNew_Serial(N);
  y = N_VNew_Serial(N);
  b = N_VNew_Serial(N);

  /* Fill matrix with uniform random data in [0,1/N] */
  int k =0;
  for (k=0; k<5*N; k++) {
    int i = rand() % N;
    int j = rand() % N;
    matdata = SUNDenseMatrix_Column(B,j);
    matdata[i] = (realtype) rand() / (realtype) RAND_MAX / N;
  }

  /* Add identity to matrix */
  fails = SUNMatScaleAddI(ONE, B);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatScaleAddI failure\n");
    return(1);
  }

  /* Fill x vector with uniform random data in [0,1] */
  xdata = N_VGetArrayPointer(x);
  int i;
  for (i=0; i<N; i++)
   xdata[i] = (realtype) rand() / (realtype) RAND_MAX;

  /* Create sparse matrix from dense, and destroy B */
  A = SUNSparseFromDenseMatrix(B, ZERO, mattype);
  SUNMatDestroy(B);

  // copy x into y to print in case of solver failure 
  N_VScale(ONE, x, y);

  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  if (fails) {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");
    return(1);
  }
  

  // D O N E :

//    /* Matrix data. */
//    int    n = 8;
//    int    ia[ 9] = { 0, 4, 7, 9, 11, 12, 15, 17, 20 };
//    int    ja[20] = { 0,    2,       5, 6, 
//                         1, 2,    4,
//                            2,             7,
//                               3,       6,
//                         1,
//                            2,       5,    7,
//                         1,             6,
//                            2,          6, 7 };
//    double  a[20] = { 7.0,      1.0,           2.0, 7.0, 
//                          -4.0, 8.0,      2.0,
//                                1.0,                     5.0,
//                                     7.0,           9.0,
//                          -4.0,
//                                7.0,           3.0,      8.0,
//                           1.0,                    11.0,
//                               -3.0,                2.0, 5.0 };
//
//    int      nnz = ia[n];
//    int      mtype = 11;        /* Real unsymmetric matrix */
//
//    /* RHS and solution vectors. */
//    double   bb[8], xx[8];
//    int      nrhs = 1;          /* Number of right hand sides. */
//
//    /* Internal solver memory pointer pt,                  */
//    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
//    /* or void *pt[64] should be OK on both architectures  */ 
//    void    *pt[64];
//
//    /* Pardiso control parameters. */
//    int      iparm[64];
//    double   dparm[64];
//    int      solver;
//    int      maxfct, mnum, phase, error, msglvl;
//
//    /* Number of processors. */
//    int      num_procs;
//
//    /* Auxiliary variables. */
//    char    *var;
//
//    error = 0;
//    solver = 0; /* use sparse direct solver */
//    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);
//
//    if (error != 0)
//    {
//        if (error == -10 )
//           printf("No license file found \n");
//        if (error == -11 )
//           printf("License is expired \n");
//        if (error == -12 )
//           printf("Wrong username or hostname \n");
//         return 1;
//    }
//    else
//        printf("[PARDISO]: License check was successful ... \n");
// 
//
//    /* Numbers of processors, value of OMP_NUM_THREADS */
//    var = getenv("OMP_NUM_THREADS");
//    if(var != NULL)
//        sscanf( var, "%d", &num_procs );
//    else {
//      num_procs = 1;
//      printf("Environment variable OMP_NUM_THREADS not set, using 1\n");
//    }
//    iparm[2]  = num_procs;
   
    

  // ________ NOT ______________ D O N E :
//    maxfct = 1;         /* Maximum number of numerical factorizations.  */
//    mnum   = 1;         /* Which factorization to use. */
//    
//    msglvl = 1;         /* Print statistical information  */
//    error  = 0;         /* Initialize error flag */
//
//
///* -------------------------------------------------------------------- */    
///* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
///*     notation.                                                        */
///* -------------------------------------------------------------------- */ 
//    for (int i = 0; i < n+1; i++) {
//        ia[i] += 1;
//    }
//    for (int i = 0; i < nnz; i++) {
//        ja[i] += 1;
//    }
//
//    /* Set right hand side to one. */
//    for (int i = 0; i < n; i++) {
//        bb[i] = i;
//    }
//
///* -------------------------------------------------------------------- */
///*  .. pardiso_chk_matrix(...)                                          */
///*     Checks the consistency of the given matrix.                      */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */
//    
//    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
//    if (error != 0) {
//        printf("\nERROR in consistency of matrix: %d", error);
//        exit(1);
//    }
//
///* -------------------------------------------------------------------- */
///* ..  pardiso_chkvec(...)                                              */
///*     Checks the given vectors for infinite and NaN values             */
///*     Input parameters (see PARDISO user manual for a description):    */
///*     Use this functionality only for debugging purposes               */
///* -------------------------------------------------------------------- */
//
//    pardiso_chkvec (&n, &nrhs, bb, &error);
//    if (error != 0) {
//        printf("\nERROR  in right hand side: %d", error);
//        exit(1);
//    }
//
///* -------------------------------------------------------------------- */
///* .. pardiso_printstats(...)                                           */
///*    prints information on the matrix to STDOUT.                       */
///*    Use this functionality only for debugging purposes                */
///* -------------------------------------------------------------------- */
//
//    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, bb, &error);
//    if (error != 0) {
//        printf("\nERROR right hand side: %d", error);
//        exit(1);
//    }
// 
///* -------------------------------------------------------------------- */    
///* ..  Reordering and Symbolic Factorization.  This step also allocates */
///*     all memory that is necessary for the factorization.              */
///* -------------------------------------------------------------------- */ 
//    phase = 11; 
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, NULL, &nrhs,
//             iparm, &msglvl, NULL, NULL, &error,  dparm);
//    
//    if (error != 0) {
//        printf("\nERROR during symbolic factorization: %d", error);
//        exit(1);
//    }
//    printf("\nReordering completed ... ");
//    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
//    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
//   
///* -------------------------------------------------------------------- */    
///* ..  Numerical factorization.                                         */
///* -------------------------------------------------------------------- */    
//    phase = 22;
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, NULL, &nrhs,
//             iparm, &msglvl, NULL, NULL, &error, dparm);
//   
//    if (error != 0) {
//        printf("\nERROR during numerical factorization: %d", error);
//        exit(2);
//    }
//    printf("\nFactorization completed ...\n ");
//
///* -------------------------------------------------------------------- */    
///* ..  Back substitution and iterative refinement.                      */
///* -------------------------------------------------------------------- */    
//    phase = 33;
//
//    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
//
//   
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, NULL, &nrhs,
//             iparm, &msglvl, bb, xx, &error,  dparm);
//   
//    if (error != 0) {
//        printf("\nERROR during solution: %d", error);
//        exit(3);
//    }
//
//    printf("\nSolve completed ... ");
//    printf("\nThe solution of the system is: ");
//    for (int i = 0; i < n; i++) {
//        printf("\n x [%d] = % f", i, xx[i] );
//    }
//    printf ("\n");
//
///* -------------------------------------------------------------------- */
///* ..  Back substitution with tranposed matrix A^t x=b                   */
///* -------------------------------------------------------------------- */
//    phase = 33;
//
//    iparm[7]  = 1;       /* Max numbers of iterative refinement steps. */
//    iparm[11] = 1;       /* Solving with transpose matrix. */
//
//    /* Set right hand side to one. */
//    for (int i = 0; i < n; i++) {
//        bb[i] = 1;
//    }
//  
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, a, ia, ja, NULL, &nrhs,
//             iparm, &msglvl, bb, xx, &error,  dparm);
//  
//    if (error != 0) {
//        printf("\nERROR during solution: %d", error);
//        exit(3);
//    }
//
//    printf("\nSolve completed ... ");
//    printf("\nThe solution of the system is: ");
//    for (int i = 0; i < n; i++) {
//        printf("\n x [%d] = % f", i, xx[i] );
//    }
//    printf ("\n");
//
///* -------------------------------------------------------------------- */    
///* ..  Convert matrix back to 0-based C-notation.                       */
///* -------------------------------------------------------------------- */ 
//    for (int i = 0; i < n+1; i++) {
//        ia[i] -= 1;
//    }
//    for (int i = 0; i < nnz; i++) {
//        ja[i] -= 1;
//    }
//
///* -------------------------------------------------------------------- */    
///* ..  Termination and release of memory.                               */
///* -------------------------------------------------------------------- */ 
//    phase = -1;                 /* Release internal memory. */
//
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, NULL, ia, ja, NULL, &nrhs,
//             iparm, &msglvl, NULL, NULL, &error,  dparm);

  printf("Creating linear solver\n");
   // Create PARDISO linear solver 
  LS = SUNPARDISO(x, A);
  
  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, RCONST(1.0e-13), 0);
 
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    printf("\nA =\n");
    SUNSparseMatrix_Print(A,stdout);
    printf("\nx (original) =\n");
    N_VPrint_Serial(y);
    printf("\nb =\n");
    N_VPrint_Serial(b);
    printf("\nx (computed) =\n");
    N_VPrint_Serial(x);
  } else {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  /* Free solver, matrix and vectors */
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector X, N_Vector Y, realtype tol)
{
  int failure = 0;
  sunindextype i, local_length, maxloc;
  realtype *Xdata, *Ydata, maxerr;
  
  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);
  local_length = N_VGetLength_Serial(X);
  
  /* check vector data */
  for(i=0; i < local_length; i++)
    failure += FNEQ(Xdata[i], Ydata[i], tol);

  if (failure > ZERO) {
    maxerr = ZERO;
    maxloc = -1;
    for(i=0; i < local_length; i++) {
      if (SUNRabs(Xdata[i]-Ydata[i]) >  maxerr) {
        maxerr = SUNRabs(Xdata[i]-Ydata[i]);
        maxloc = i;
      }
    }
    printf("check err failure: maxerr = %g at loc %li (tol = %g)\n",
	   maxerr, (long int) maxloc, tol);
    return(1);
  }
  else
    return(0);
}
