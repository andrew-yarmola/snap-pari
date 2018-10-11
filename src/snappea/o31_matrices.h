#ifndef _o31_matrices_
#define _o31_matrices_

#include "complex.h"

typedef int FuncResult;

/*
 *  Matrices in O(3,1) represent isometries in the Minkowski space
 *  model of hyperbolic 3-space.  The matrices are expressed relative
 *  to a coordinate system in which the metric is
 *
 *                          -1  0  0  0
 *                           0  1  0  0
 *                           0  0  1  0
 *                           0  0  0  1
 *
 *  That is, the first coordinate is timelike, and the remaining
 *  three are spacelike.  O(3,1) matrices represent both
 *  orientation_preserving and orientation_reversing isometries.
 */

typedef double O31Matrix[4][4];
typedef double GL4RMatrix[4][4];

/*
 *  An O31Vector is a vector in (3,1)-dimensional Minkowski space.
 *  The 0-th coordinate is the timelike one.
 */

typedef double O31Vector[4];

/*
 *  MatrixInt22 is a 2 x 2 integer matrix.  A MatrixInt22
 *  may, for example, describe how the peripheral curves of
 *  one Cusp map to those of another.
 */

/************************************************************************/
/*                                                                      */
/*                          o31_matrices.c                              */
/*                                                                      */
/************************************************************************/

extern double       gl4R_determinant(const GL4RMatrix m);
extern double       o31_trace(const O31Matrix m);

extern O31Matrix    O31_identity;

extern void         o31_copy(O31Matrix dest, const O31Matrix source);
extern void         o31_invert(const O31Matrix m, O31Matrix m_inverse);
extern FuncResult   gl4R_invert(const GL4RMatrix m, GL4RMatrix m_inverse);
extern void         o31_product(const O31Matrix a, const O31Matrix b, O31Matrix product);
extern Boolean      o31_equal(const O31Matrix a, const O31Matrix b, double epsilon);
extern double       o31_deviation(const O31Matrix m);
extern void         o31_GramSchmidt(O31Matrix m);
extern void         o31_conjugate(const O31Matrix m, const O31Matrix t, O31Matrix Tmt);
extern double       o31_inner_product(const O31Vector u, const O31Vector v);
extern void         o31_matrix_times_vector(const O31Matrix m, const O31Vector v, O31Vector product);
extern void         o31_constant_times_vector(double r, const O31Vector v, O31Vector product);
extern void         o31_copy_vector(O31Vector dest, const O31Vector source);
extern void         o31_vector_sum(const O31Vector a, const O31Vector b, O31Vector sum);
extern void         o31_vector_diff(const O31Vector a, const O31Vector b, O31Vector diff);
/*
 *  These functions all do what you would expect.
 *  o31_conjugate() replaces m with (t^-1) m t.
 */

#endif
