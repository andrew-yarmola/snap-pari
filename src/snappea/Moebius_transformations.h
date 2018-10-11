#ifndef _Moebius_transformations_
#define _Moebius_transformations_

#include "complex.h"
#include "sl2c_matrices.h"

#ifndef NOPARI
#include "pariclass.hh"
#endif

/*
 *  SnapPea represents a Moebius transformation as a matrix
 *  in SL(2,C) plus a specification of whether the Moebius
 *  transformation is orientation_preserving or orientation_reversing.
 *
 *  If mt->parity is orientation_preserving, then mt->matrix is
 *  interpreted in the usual way as the Moebius transformation
 *
 *                      az + b
 *              f(z) = --------
 *                      cz + d
 *
 *
 *  If mt->parity is orientation_reversing, then mt->matrix is
 *  interpreted as a function of the complex conjugate z' ("z-bar")
 *
 *                      az' + b
 *              f(z) = ---------
 *                      cz' + d
 */

/*
 *  The values of MatrixParity should not be changed.
 *  (They must correspond to the values in the parity[] table in tables.c.)
 */

typedef int MatrixParity;
enum
{
    orientation_reversing = 0,
    orientation_preserving = 1
};


struct MoebiusTransformation
{
    SL2CMatrix      matrix;
    MatrixParity    parity;

#ifndef NOPARI  
    pari            acc_matrix;
#endif

    MoebiusTransformation() {}
    MoebiusTransformation(const Complex& a, const Complex& b, 
			  const Complex& c, const Complex& d, 
			  MatrixParity par = orientation_preserving); 
    MoebiusTransformation(const Complex& diag, MatrixParity par = orientation_preserving);
    MoebiusTransformation(const Complex a[3], const Complex b[3], 
			  MatrixParity par = orientation_preserving);
    MoebiusTransformation(const Complex a[3]); // Map to 0,1,Infinity. 
    MoebiusTransformation(const Complex& a, const Complex& b, 
			  const Complex& c);


    void copy(const MoebiusTransformation& mt);

    MoebiusTransformation(const MoebiusTransformation& mt)
      { copy(mt); }
    MoebiusTransformation& operator = (const MoebiusTransformation& mt)
      { copy(mt); return *this; }

    ~MoebiusTransformation() {}

    Complex z_transform_length() const; 
};

void Moebius_copy(MoebiusTransformation *dest, const MoebiusTransformation *source);

void Moebius_invert(const MoebiusTransformation *mt,
		    MoebiusTransformation *mt_inverse);
void Moebius_product(const MoebiusTransformation *a,
		     const MoebiusTransformation *b,
		     MoebiusTransformation *product);

Boolean read(MoebiusTransformation& m, FILE *in = stdin);

MoebiusTransformation get_transform(const Complex a[4], const Complex b[4]);

MoebiusTransformation Moebius_z_trans(const Complex& d);
MoebiusTransformation Moebius_x_trans(const Complex& d);

void print(const MoebiusTransformation& mt); 

Boolean is_identity(const MoebiusTransformation& mt, double eps=1e-7); 
double id_distance(const MoebiusTransformation& mt);

/*
 *  These functions do what you would expect.
 */


inline Complex operator * (const MoebiusTransformation& mt, Complex z)
{
  return apply_matrix(mt.matrix, mt.parity ? z : complex_conjugate(z)); 
}

#ifndef NOPARI
pari operator * (const MoebiusTransformation& mt, pari const& z);
#endif

inline MoebiusTransformation operator * (const MoebiusTransformation& a, 
					 const MoebiusTransformation& b)
{
  MoebiusTransformation product;
  Moebius_product(&a,&b,&product);
  return product; 
}

inline MoebiusTransformation& operator *= (MoebiusTransformation& a, 
					   const MoebiusTransformation& b)
{
  MoebiusTransformation product;
  Moebius_product(&a,&b,&product);
  a = product; 
  return a; 
}

inline MoebiusTransformation inverse(const MoebiusTransformation& mt)
{
  MoebiusTransformation inverse;
  Moebius_invert(&mt,&inverse);
  return inverse; 
}

std::ostream& operator << (std::ostream& out, MoebiusTransformation const& mt);

#endif
