/*
 *  Moebius_transformations.c
 */

#include "Moebius_transformations.h"
#ifndef NOPARI  
#include "pariwrap.hh"
#endif
#include <algorithm>

using std::ostream;
using std::max;

// const trace_init ti("Moebius_transformations.c"); 

MoebiusTransformation::MoebiusTransformation(const Complex& diag, MatrixParity par)
{ 
  matrix[0][0] = diag;           matrix[0][1] = Complex(0.,0.); 
  matrix[1][0] = Complex(0.,0.); matrix[1][1] = Complex(1.,0.)/diag; 
  parity = par;
#ifndef NOPARI  
  acc_matrix = idmat(2);
#endif
}

MoebiusTransformation::MoebiusTransformation(const Complex& a, const Complex& b, 
					     const Complex& c, const Complex& d, 
					     MatrixParity par)
{
  matrix[0][0] = a; matrix[0][1] = b; matrix[1][0] = c; matrix[1][1] = d; 
  parity = par; 
#ifndef NOPARI  
  acc_matrix = idmat(2);
#endif
}

void MoebiusTransformation::copy(const MoebiusTransformation &mt)
{
    sl2c_copy(matrix, mt.matrix);
    parity = mt.parity;

    /* copy accurate part if set */
#ifndef NOPARI  
    if (mt.acc_matrix.type()!=1) 
      acc_matrix = gcopy(mt.acc_matrix);
#endif  
}


void Moebius_invert(
    const MoebiusTransformation *mt,
    MoebiusTransformation   *mt_inverse)
{
    sl2c_invert(mt->matrix, mt_inverse->matrix);

    if (mt->parity == orientation_reversing)
	sl2c_complex_conjugate(mt_inverse->matrix, mt_inverse->matrix);

    /* invert accurate part if set */
#ifndef NOPARI
    if (mt->acc_matrix.type()!=1) {
      mt_inverse->acc_matrix = (mt->parity == orientation_preserving) ?
	adj(mt->acc_matrix) : gconj(adj(mt->acc_matrix));
    }
#endif  

    mt_inverse->parity = mt->parity;
}

void Moebius_copy(MoebiusTransformation *dest, const MoebiusTransformation *source)
{ dest->copy(*source); }

void Moebius_product(
    const MoebiusTransformation *a,
    const MoebiusTransformation *b,
    MoebiusTransformation   *product)
{
    SL2CMatrix  factor1,
		factor2;

    sl2c_copy(factor1, a->matrix);
    sl2c_copy(factor2, b->matrix);

    if (a->parity == orientation_reversing) 
	sl2c_complex_conjugate(factor2, factor2);

    sl2c_product(factor1, factor2, product->matrix);
    
    /* multiply accurate part if set */
#ifndef NOPARI
    if (a->acc_matrix.type()!=1 && b->acc_matrix.type()!=1) {
      pari acc_factor2 = (a->parity == orientation_preserving) ?
			 b->acc_matrix : gconj(b->acc_matrix);
      product->acc_matrix = a->acc_matrix * acc_factor2; 
    }
#endif

    product->parity = (a->parity == b->parity) ?
		      orientation_preserving:
		      orientation_reversing;
}

MoebiusTransformation::MoebiusTransformation(const Complex u[3])
: parity(orientation_preserving)
{ 
  Complex u10 = u[1] - u[0]; 
  Complex u02 = u[0] - u[2]; 
  Complex u12 = u[1] - u[2]; 
  Complex a = complex_sqrt(u10/(u12 * u02));
  Complex c = a * u10/u12; 

  matrix[0][0] = a;
  matrix[0][1] = -a * u[0]; 
  matrix[1][0] = c;
  matrix[1][1] = -c * u[2]; 
}

MoebiusTransformation::MoebiusTransformation(const Complex& u, const Complex& v, const Complex& w)
: parity(orientation_preserving)
{
  Complex u10 = v - u; 
  Complex u02 = u - w; 
  Complex u12 = v - w; 
  Complex a = complex_sqrt(u10/(u12 * u02));
  Complex c = a * u10/u12; 

  matrix[0][0] = a;
  matrix[0][1] = -a * u; 
  matrix[1][0] = c;
  matrix[1][1] = -c * w; 
}


MoebiusTransformation::MoebiusTransformation(const Complex a[3], const Complex b[3], MatrixParity par)
: parity(par)
{
  Complex an[3], ad[3], bn[3], bd[3]; 

  int i;
  for (i = 0; i < 3; i++) {
    if (complex_big(a[i])) {
      an[i] = Complex(1.,0.); 
      ad[i] = Complex(0.,0.); 
    } else {
      an[i] = a[i];
      ad[i] = Complex(1.,0.); 
    }
    if (complex_big(b[i])) {
      bn[i] = Complex(1.,0.); 
      bd[i] = Complex(0.,0.); 
    } else {
      bn[i] = b[i];
      bd[i] = Complex(1.,0.); 
    }
  }

    /* No, I didn't get these by hand, I used Mathematica. */

    matrix[0][0] =
      -(ad[1]*ad[2]*an[0]*bd[2]*bn[0]*bn[1]) + 
	ad[0]*ad[2]*an[1]*bd[2]*bn[0]*bn[1] + 
	  ad[1]*ad[2]*an[0]*bd[1]*bn[0]*bn[2] - 
	    ad[0]*ad[1]*an[2]*bd[1]*bn[0]*bn[2] - 
	      ad[0]*ad[2]*an[1]*bd[0]*bn[1]*bn[2] + 
		ad[0]*ad[1]*an[2]*bd[0]*bn[1]*bn[2];

    matrix[0][1] = 
      ad[1]*an[0]*an[2]*bd[2]*bn[0]*bn[1] - 
	ad[0]*an[1]*an[2]*bd[2]*bn[0]*bn[1] - 
	  ad[2]*an[0]*an[1]*bd[1]*bn[0]*bn[2] + 
	    ad[0]*an[1]*an[2]*bd[1]*bn[0]*bn[2] + 
	      ad[2]*an[0]*an[1]*bd[0]*bn[1]*bn[2] - 
		ad[1]*an[0]*an[2]*bd[0]*bn[1]*bn[2]; 
 
    matrix[1][0] = 
      ad[0]*ad[2]*an[1]*bd[1]*bd[2]*bn[0] - 
	ad[0]*ad[1]*an[2]*bd[1]*bd[2]*bn[0] - 
	  ad[1]*ad[2]*an[0]*bd[0]*bd[2]*bn[1] + 
	    ad[0]*ad[1]*an[2]*bd[0]*bd[2]*bn[1] + 
	      ad[1]*ad[2]*an[0]*bd[0]*bd[1]*bn[2] - 
		ad[0]*ad[2]*an[1]*bd[0]*bd[1]*bn[2]; 

    matrix[1][1] =  
      -(ad[2]*an[0]*an[1]*bd[1]*bd[2]*bn[0]) + 
	ad[1]*an[0]*an[2]*bd[1]*bd[2]*bn[0] + 
	  ad[2]*an[0]*an[1]*bd[0]*bd[2]*bn[1] - 
	    ad[0]*an[1]*an[2]*bd[0]*bd[2]*bn[1] - 
	      ad[1]*an[0]*an[2]*bd[0]*bd[1]*bn[2] + 
		ad[0]*an[1]*an[2]*bd[0]*bd[1]*bn[2];

  Complex root_det = complex_sqrt(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);

  if (complex_small(root_det)) {
    fprintf(stderr, "attempt to construct Moebius transformation mapping between degenerate triples\n");
    return;
  }

  matrix[0][0] /= root_det;
  matrix[0][1] /= root_det;
  matrix[1][0] /= root_det;
  matrix[1][1] /= root_det;
}

Boolean is_identity(const MoebiusTransformation& mt, double eps)
{
  if (!mt.parity) return FALSE;
  double esq = eps*eps; 
  if (complex_modulus_squared(mt.matrix[0][1]) > eps) return FALSE; 
  if (complex_modulus_squared(mt.matrix[1][0]) > eps) return FALSE; 
  if (complex_modulus_squared(mt.matrix[0][0]-Complex(1.,0.)) > eps) return FALSE; 
  if (complex_modulus_squared(mt.matrix[1][1]-Complex(1.,0.)) > eps) return FALSE; 
  return TRUE; 
}

#ifndef NOPARI
pari operator * (const MoebiusTransformation& mt, pari const& z)
{
  static pari big = real(1)/pEPSILON;
  // printf("Mob:"); z.print(); 
  if (gabs(z) > big) {
    // mt.acc_matrix[0][1].print(); 
    if (gabs(mt.acc_matrix[0][1]) < pEPSILON) {
      // printf("Mob:infinity\n");
      return pHUGE; 
    }
    return mt.acc_matrix[0][0]/mt.acc_matrix[0][1]; /* a/c */ 
  }
  pari denom = mt.acc_matrix[0][1] * z + mt.acc_matrix[1][1]; 

  // denom.print();
  if (gabs(denom) < pEPSILON) {
    // printf("Mob:infinity\n");
    return pHUGE; 
  }

  /* result is (az + b)/(cz + d). */ 

  pari res = (mt.acc_matrix[0][0] * z + mt.acc_matrix[1][0])/denom; 
  return (mt.parity==orientation_preserving) ? res : gconj(res); 
}
#endif

MoebiusTransformation get_transform(const Complex a[4], const Complex b[4])
{
  MoebiusTransformation mt; 
  Complex A[4];
  int j; 

  if (!same_point(cross_ratio(a), cross_ratio(b))) {
    if (!same_point(cross_ratio(a), complex_conjugate(cross_ratio(b)))) {
      fprintf(stderr, "attempt to find moebius transformation between non congruent tetrahedra\n"); 
      return MoebiusTransformation(Complex(1.,0.));
    }
    mt.parity = orientation_reversing; 
    for (j=0; j<4; j++) A[j] = complex_conjugate(a[j]); 
  } else {
    mt.parity = orientation_preserving;
    for (j=0; j<4; j++) A[j] = a[j]; 
  }

  Complex k = (b[2] - b[0]) * (A[2] - A[1]) / ((b[2] - b[1]) * (A[2] - A[0])); 

  Complex b1k = b[1] * k; 

  Complex normalization = complex_sqrt(Complex(1.,0.) / (k * (A[1] - A[0]) * (b[1] - b[0]))); 

  mt.matrix[0][0] = normalization * (b1k - b[0]); 
  mt.matrix[0][1] = normalization * (b[0] * A[1] - b1k * A[0]); 
  mt.matrix[1][0] = normalization * (k - Complex(1.,0.)); 
  mt.matrix[1][1] = normalization * (A[1] - k * A[0]); 

  return mt; 
}

// Checks if the MoebiusTransformation is complex 
// multiplication then returns its complex length. 

Complex MoebiusTransformation::z_transform_length() const
{
  if (!complex_small(matrix[0][1]) || !complex_small(matrix[1][0])) {
    printf("z_transform_length() applied to non_diagonal MoebiusTransformation\n"); 
    return Complex(0.,0.); 
  }

  return (complex_log(matrix[0][0],0.0) - complex_log(matrix[1][1],0.0))/Complex(2.,0.); 
}

MoebiusTransformation Moebius_z_trans(const Complex& d)
{
  return MoebiusTransformation(complex_exp(d/2.0));
}

MoebiusTransformation Moebius_x_trans(const Complex& d)
{
  Complex a = complex_cosh(d/2.0); 
  Complex b = complex_sinh(d/2.0); 
  return MoebiusTransformation(a,b,b,a); 
}

void print(const MoebiusTransformation& mt)
{
  printf("["); 
  fwprint(stdout, mt.matrix[0][0], 8); 
  printf(","); 
  fwprint(stdout, mt.matrix[0][1], 8); 
  printf(";"); 
  fwprint(stdout, mt.matrix[1][0], 8); 
  printf(","); 
  fwprint(stdout, mt.matrix[1][1], 8); 
  printf("]");
  if (!mt.parity) printf("*"); 
}

ostream& operator << (ostream& out, MoebiusTransformation const& mt)
{
  out << '[' << mt.matrix[0][0] << ',' << mt.matrix[0][1] << ';' << 
		mt.matrix[1][0] << ',' << mt.matrix[1][1] << ']';
  if (!mt.parity) out << '*';
  return out; 
}

Boolean read(MoebiusTransformation& m, FILE *in)
{
  char s[2]; 

  fscanf(in, " %1s ", s); if (s[0] != '[') return 0;
  if (!read(m.matrix[0][0], in)) return 0;
  fscanf(in, " %1s ", s); if (s[0] != ',') return 0;
  if (!read(m.matrix[0][1], in)) return 0;
  fscanf(in, " %1s ", s); if (s[0] != ';') return 0;
  if (!read(m.matrix[1][0], in)) return 0;
  fscanf(in, " %1s ", s); if (s[0] != ',') return 0;
  if (!read(m.matrix[1][1], in)) return 0;
  fscanf(in, " ; %d %1s", &m.parity, s);

  return (s[0]==']'); 
}

double id_distance(MoebiusTransformation const& m)
{
  double d = 0., zd;
  if (m.parity != orientation_preserving)
    d = 1.; // Don't want to return zero for z->conj(z).

  Complex z; 
  int i, j; 
  for (i=0; i<2; ++i) {
    for (j=0; j<2; ++j) {
      z = m.matrix[i][j];

      // Adjust diagonal terms so that identity becomes zero.
      // [-1,0;0,-1] also represents the identity.
      if (i==j) z.real += ((z.real > 0.) ? -1.:1.);

      zd = max(fabs(z.real),fabs(z.imag));
      if (zd > d) d = zd; 
    }
  }
  return d; 
}
