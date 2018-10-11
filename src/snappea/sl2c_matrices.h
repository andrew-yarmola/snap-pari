#ifndef _sl2c_matrices_
#define _sl2c_matrices_

#include "complex.h"

typedef Complex SL2CMatrix[2][2];

/************************************************************************/
/*                                                                      */
/*                              sl2c_matrices.c                         */
/*                                                                      */
/************************************************************************/

extern void     sl2c_copy(SL2CMatrix dest, const SL2CMatrix source);
extern void     sl2c_invert(const SL2CMatrix a, SL2CMatrix inverse);
extern void     sl2c_complex_conjugate(const SL2CMatrix a, SL2CMatrix conjugate);
void        sl2c_product(const SL2CMatrix a, const SL2CMatrix b, SL2CMatrix product);
extern void     sl2c_adjoint(const SL2CMatrix a, SL2CMatrix adjoint);
extern void     sl2c_normalize(SL2CMatrix a);
extern Complex  sl2c_determinant(const SL2CMatrix m);
extern Boolean  sl2c_matrix_is_real(const SL2CMatrix a);

Complex apply_matrix(const SL2CMatrix matrix, Complex z);

inline Complex operator * (const SL2CMatrix matrix, Complex z)
{
  return apply_matrix(matrix, z); 
}

#endif
