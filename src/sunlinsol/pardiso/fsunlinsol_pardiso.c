/*
 * -----------------------------------------------------------------
 * Programmer:  Raymond Langer @ RWTH Aachen University
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 * This file (companion of fsunlinsol_pardiso.h) contains the
 * implementation needed for the Fortran initialization of pardiso
 * linear solver operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunlinsol_pardiso.h"

/* Define global linsol variables */

SUNLinearSolver F2C_CVODE_linsol;
SUNLinearSolver F2C_IDA_linsol;
SUNLinearSolver F2C_KINSOL_linsol;
SUNLinearSolver F2C_ARKODE_linsol;
SUNLinearSolver F2C_ARKODE_mass_sol;

/* Declarations of external global variables */

extern SUNMatrix F2C_CVODE_matrix;
extern SUNMatrix F2C_IDA_matrix;
extern SUNMatrix F2C_KINSOL_matrix;
extern SUNMatrix F2C_ARKODE_matrix;
extern SUNMatrix F2C_ARKODE_mass_matrix;

extern N_Vector F2C_CVODE_vec;
extern N_Vector F2C_IDA_vec;
extern N_Vector F2C_KINSOL_vec;
extern N_Vector F2C_ARKODE_vec;

/* Fortran callable interfaces */

void FSUNPARDISO_INIT(int *code, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (F2C_CVODE_linsol)  SUNLinSolFree(F2C_CVODE_linsol);
    F2C_CVODE_linsol = NULL;
    F2C_CVODE_linsol = SUNPARDISO(F2C_CVODE_vec, F2C_CVODE_matrix);
    if (F2C_CVODE_linsol == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    if (F2C_IDA_linsol)  SUNLinSolFree(F2C_IDA_linsol);
    F2C_IDA_linsol = NULL;
    F2C_IDA_linsol = SUNPARDISO(F2C_IDA_vec, F2C_IDA_matrix);
    if (F2C_IDA_linsol == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    if (F2C_KINSOL_linsol)  SUNLinSolFree(F2C_KINSOL_linsol);
    F2C_KINSOL_linsol = NULL;
    F2C_KINSOL_linsol = SUNPARDISO(F2C_KINSOL_vec, F2C_KINSOL_matrix);
    if (F2C_KINSOL_linsol == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    if (F2C_ARKODE_linsol)  SUNLinSolFree(F2C_ARKODE_linsol);
    F2C_ARKODE_linsol = NULL;
    F2C_ARKODE_linsol = SUNPARDISO(F2C_ARKODE_vec, F2C_ARKODE_matrix);
    if (F2C_ARKODE_linsol == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNPARDISO_REINIT(int *code, long int *NNZ, int *reinit_type, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNPARDISOReInit(F2C_CVODE_linsol, F2C_CVODE_matrix,
                        *NNZ, *reinit_type);
    break;
  case FCMIX_IDA:
    *ier = SUNPARDISOReInit(F2C_IDA_linsol, F2C_IDA_matrix,
                        *NNZ, *reinit_type);
    break;
  case FCMIX_KINSOL:
    *ier = SUNPARDISOReInit(F2C_KINSOL_linsol, F2C_KINSOL_matrix,
                        *NNZ, *reinit_type);
    break;
  case FCMIX_ARKODE:
    *ier = SUNPARDISOReInit(F2C_ARKODE_linsol, F2C_ARKODE_matrix,
                        *NNZ, *reinit_type);
    break;
  default:
    *ier = -1;
  }
}


void FSUNPARDISO_SETORDERING(int *code, int *ordering_choice, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNPARDISOSetOrdering(F2C_CVODE_linsol, *ordering_choice);
    break;
  case FCMIX_IDA:
    *ier = SUNPARDISOSetOrdering(F2C_IDA_linsol, *ordering_choice);
    break;
  case FCMIX_KINSOL:
    *ier = SUNPARDISOSetOrdering(F2C_KINSOL_linsol, *ordering_choice);
    break;
  case FCMIX_ARKODE:
    *ier = SUNPARDISOSetOrdering(F2C_ARKODE_linsol, *ordering_choice);
    break;
  default:
    *ier = -1;
  }
}


void FSUNMASSPARDISO_INIT(int *ier)
{
  *ier = 0;
  if (F2C_ARKODE_mass_sol)  SUNLinSolFree(F2C_ARKODE_mass_sol);
  F2C_ARKODE_mass_sol = NULL;
  F2C_ARKODE_mass_sol = SUNPARDISO(F2C_ARKODE_vec, 
                               F2C_ARKODE_mass_matrix);
  if (F2C_ARKODE_mass_sol == NULL) *ier = -1;
}


void FSUNMASSPARDISO_REINIT(long int *NNZ, int *reinit_type, int *ier)
{
  *ier = 0;
  *ier = SUNPARDISOReInit(F2C_ARKODE_mass_sol, F2C_ARKODE_mass_matrix,
                      *NNZ, *reinit_type);
}


void FSUNMASSPARDISO_SETORDERING(int *ordering_choice, int *ier)
{
  *ier = 0;
  *ier = SUNPARDISOSetOrdering(F2C_ARKODE_mass_sol, *ordering_choice);
}
