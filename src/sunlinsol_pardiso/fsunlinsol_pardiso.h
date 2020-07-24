/*
 * -----------------------------------------------------------------
 * Programmer:  Raymond Langer @ RWTH Aachen University
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 * This file (companion of fsunlinsol_pardiso.c) contains the
 * definitions needed for the initialization of pardiso
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_PARDISO_H
#define _FSUNLINSOL_PARDISO_H

#include <sunlinsol/sunlinsol_pardiso.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNPARDISO_INIT            SUNDIALS_F77_FUNC(fsunpardisoinit,        FSUNPARDISOINIT)
#define FSUNPARDISO_REINIT          SUNDIALS_F77_FUNC(fsunpardisoreinit,      FSUNPARDISOREINIT)
#define FSUNPARDISO_SETORDERING     SUNDIALS_F77_FUNC(fsunpardisosetordering, FSUNPARDISOSETORDERING)
#define FSUNMASSPARDISO_INIT        SUNDIALS_F77_FUNC(fsunmasspardisoinit,        FSUNMASSPARDISOINIT)
#define FSUNMASSPARDISO_REINIT      SUNDIALS_F77_FUNC(fsunmasspardisoreinit,      FSUNMASSPARDISOREINIT)
#define FSUNMASSPARDISO_SETORDERING SUNDIALS_F77_FUNC(fsunmasspardisosetordering, FSUNMASSPARDISOSETORDERING)
#else
#define FSUNPARDISO_INIT            fsunpardisoinit_
#define FSUNPARDISO_REINIT          fsunpardisoreinit_
#define FSUNPARDISO_SETORDERING     fsunpardisosetordering_
#define FSUNMASSPARDISO_INIT        fsunmasspardisoinit_
#define FSUNMASSPARDISO_REINIT      fsunmasspardisoreinit_
#define FSUNMASSPARDISO_SETORDERING fsunmasspardisosetordering_
#endif


/* Declarations of global variables */

extern SUNLinearSolver F2C_CVODE_linsol;
extern SUNLinearSolver F2C_IDA_linsol;
extern SUNLinearSolver F2C_KINSOL_linsol;
extern SUNLinearSolver F2C_ARKODE_linsol;
extern SUNLinearSolver F2C_ARKODE_mass_sol;

/* 
 * Prototypes of exported functions 
 *
 * FSUNPARDISO_INIT - initializes pardiso linear solver for main problem
 * FSUNPARDISO_REINIT - reinitializes pardiso linear solver for main problem
 * FSUNPARDISO_SETORDERING - sets the ordering choice used by PARDISO for main problem
 * FSUNMASSPARDISO_INIT - initializes pardiso linear solver for mass matrix solve
 * FSUNMASSPARDISO_REINIT - reinitializes pardiso linear solver for mass matrix solve
 * FSUNMASSPARDISO_SETORDERING - sets the ordering choice used by PARDISO for mass matrix solve
 */

void FSUNPARDISO_INIT(int *code, int *ier);
void FSUNPARDISO_REINIT(int *code, long int *NNZ, 
                    int *reinit_type, int *ier);
void FSUNPARDISO_SETORDERING(int *code, int *ordering,
                         int *ier);
void FSUNMASSPARDISO_INIT(int *ier);
void FSUNMASSPARDISO_REINIT(long int *NNZ, 
                        int *reinit_type, int *ier);
void FSUNMASSPARDISO_SETORDERING(int *ordering, int *ier);

#ifdef __cplusplus
}
#endif

#endif
