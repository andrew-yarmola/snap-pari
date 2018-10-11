#ifndef _matrix_conversion_
#define _matrix_conversion_

#include "o31_matrices.h"
#include "Moebius_transformations.h"


/************************************************************************/
/*                                                                      */
/*                          matrix_conversion.c                         */
/*                                                                      */
/************************************************************************/

extern void Moebius_to_O31(const MoebiusTransformation *A, O31Matrix B);
extern void O31_to_Moebius(const O31Matrix B, MoebiusTransformation *A);
/*
 *  Convert matrices back and forth between SL(2,C) and O(3,1).
 */

extern void Moebius_array_to_O31_array( const MoebiusTransformation arrayA[],
					O31Matrix               arrayB[],
					int                     num_matrices);
extern void O31_array_to_Moebius_array( const O31Matrix         arrayB[],
					MoebiusTransformation   arrayA[],
					int                     num_matrices);
/*
 *  Convert arrays of matrices back and forth between SL(2,C) and O(3,1).
 */

extern Boolean O31_determinants_OK( const O31Matrix   arrayB[],
				    int         num_matrices,
				    double      epsilon);
/*
 *  Returns TRUE if all the O31Matrices in the array have determinants
 *  within epsilon of plus or minus one, and FALSE otherwise.
 */


#endif
