/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_coderRand_api.h
 *
 * Code generation for function '_coder_coderRand_api'
 *
 */

#ifndef _CODER_CODERRAND_API_H
#define _CODER_CODERRAND_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_coderRand_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern real_T coderRand(void);
extern void coderRand_api(int32_T nlhs, const mxArray *plhs[1]);
extern void coderRand_atexit(void);
extern void coderRand_initialize(void);
extern void coderRand_terminate(void);
extern void coderRand_xil_terminate(void);

#endif

/* End of code generation (_coder_coderRand_api.h) */
