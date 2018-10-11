#ifndef _SL2RMatrix_
#define _SL2RMatrix_

#include "kernel.h"
#include "o31_matrices.h"

typedef double SL2RMatrix[2][2];
typedef double O21Matrix[3][3];

extern void SL2R_to_o21(SL2RMatrix,O21Matrix);
extern void o21_to_o31(O21Matrix, O31_matrix&);
extern void SL2RMatrix_product(SL2RMatrix , SL2RMatrix, SL2RMatrix);
extern void SL2RMatrix_transpose(SL2RMatrix ,SL2RMatrix);
extern void SL2RMatrix_conjugate(SL2RMatrix ,SL2RMatrix, SL2RMatrix);
#endif
