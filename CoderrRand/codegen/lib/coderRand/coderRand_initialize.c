/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * coderRand_initialize.c
 *
 * Code generation for function 'coderRand_initialize'
 *
 */

/* Include files */
#include "coderRand.h"
#include "coderRand_initialize.h"
#include "eml_rand_shr3cong_stateful.h"
#include "eml_rand_mcg16807_stateful.h"
#include "eml_rand.h"
#include "eml_rand_mt19937ar_stateful.h"

/* Function Definitions */
void coderRand_initialize(void)
{
  state_not_empty_init();
  eml_rand_init();
  eml_rand_mcg16807_stateful_init();
  eml_rand_shr3cong_stateful_init();
}

/* End of code generation (coderRand_initialize.c) */
