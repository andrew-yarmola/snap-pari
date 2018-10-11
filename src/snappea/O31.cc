#include "O31.h"
#include <iomanip>
#include <cstdio>

using std::ostream;
using std::ios;

void GL4R_matrix::copy(const GL4R_matrix& m)
{
  int i,j;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      mx[i][j] = m.mx[i][j];
}

GL4R_matrix::GL4R_matrix(int) // constructs identity matrix
{
  int i,j;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      mx[i][j] = (i == j) ? 1.0 : 0.0;
}

O31_vector O31_matrix::row(int r) const
{
  O31_vector v;
  int i; 
  for (i=0; i<4; i++) v[i] = mx[r][i];
  return v; 
}

O31_vector O31_matrix::col(int c) const
{
  O31_vector v;
  int i; 
  for (i=0; i<4; i++) v[i] = mx[i][c];
  return v; 
}

double GL4R_matrix::distance_from_identity() const
{
  int i,j;
  double d, maxd = 0.0; 
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      d = fabs(mx[i][j] - ((i == j) ? 1.0 : 0.0));
      if (d > maxd) maxd = d; 
    }
  }
  return maxd;
}

#if 0
R4_vector operator * (const GL4R_matrix& m, const R4_vector& v)
{ 
  O31_vector mv; 
  o31_matrix_times_vector(m.mx, v, mv); 
  return mv; 
}

O31_vector operator * (const O31_matrix& m, const O31_vector& v) 
{ 
  O31_vector mv; 
  o31_matrix_times_vector(m.mx, v, mv); 
  return mv; 
}
#endif

/*
 *  gl4R_invert will consider a matrix to be singular iff one of the
 *  absolute value of one of the pivots is less than SINGULAR_MATRIX_EPSILON.
 */
#define SINGULAR_MATRIX_EPSILON     double(1e-2)

#define COLUMN_PRODUCT(m, i, j)     \
    (-m[0][i]*m[0][j] + m[1][i]*m[1][j] + m[2][i]*m[2][j] + m[3][i]*m[3][j])

GL4R_matrix inverse(const GL4R_matrix& m)
{
  GL4R_matrix m_inverse; 
  if (!compute_inverse(m,m_inverse)) {
    fprintf(stderr, "inversion of GL4R_matrix failed\n"); 
    return m;
  }
  return m_inverse; 
}

GL4R_matrix transpose(const GL4R_matrix& m)
{
  GL4R_matrix t;
  int i, j; 
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      t.mx[i][j] = m.mx[j][i];
    }
  }
  return t; 
}

void GL4R_matrix::print(FILE* fp) const
{
  int i, j; 

  fprintf(fp, "[");
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      fwprint(fp, mx[i][j], 8);
      if (j<3) fprintf(fp, ","); 
      else if (i<3) fprintf(fp, ";\n "); 
    }
  }
  fprintf(fp,"]\n"); 
}

ostream& operator << (ostream& out, const GL4R_matrix& m)
{
  int i, j, w; 
  ios::fmtflags old = out.setf(ios::fixed); 
  w = out.precision() + 2; 

  out << "[";
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      fwprint(out, m.mx[i][j], w);
      if (j<3) out << ", "; 
      else if (i<3) out << ";\n "; 
    }
  }
  out << "]\n"; 

  out.flags(old);
  return out; 
}


int GL4R_matrix::read(FILE* fp)
{
  char c[2]; 
  int i, j; 

  fscanf(fp, " %1s ", c); if (c[0] != '[') return 0; 
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      if (j<3) {
	if (fscanf(fp, " %lf , ", &mx[i][j]) != 1) return 0; 
      } else {
	fscanf(fp, " %lf %1s \n", &mx[i][j], c); 
	if ((i<3 && c[0] != ';') || (i==3 && c[0] != ']')) return 0; 
      }
    }
  }
  return 1; 
}


void R4_vector::copy(const R4_vector& v)
{
  register int  i;

  for (i = 0; i < 4; i++)
    vec[i] = v.vec[i]; 
}

R4_vector::R4_vector(double x, double y, double z, double w)
{
  vec[0] = x; vec[1] = y; vec[2] = z; vec[3] = w; 
}

void R4_vector::operator += (const R4_vector& rhs)
{
  int i;
  for (i=0; i<4; i++)
    vec[i] += rhs.vec[i]; 
}

void R4_vector::operator -= (const R4_vector& rhs)
{
  int i;
  for (i=0; i<4; i++)
    vec[i] -= rhs.vec[i]; 
}

void R4_vector::operator *= (double r)
{
  int i;
  for (i=0; i<4; i++)
    vec[i] *= r; 
}

void R4_vector::operator /= (double r)
{
  int i;
  for (i=0; i<4; i++)
    vec[i] /= r; 
}

void R4_vector::negate()
{
  int i;
  for (i=0; i<4; i++)
    vec[i] = -vec[i]; 
}

double operator * (const R4_vector& a, const R4_vector& b) // Euclidean dot
{
  int i; 
  double dot = 0.0; 

  for (i=0; i<4; i++)
    dot += a[i] * b[i];
  return dot; 
}

void R4_vector::print(FILE* fp) const
{
  fprintf(fp, "[%lf, %lf, %lf, %lf]", vec[0], vec[1], vec[2], vec[3]);
}

ostream& operator << (ostream& out, const R4_vector& v)
{
  return out << "[" << v.vec[0] << ", " << v.vec[1] << ", " << v.vec[2] << ", " << 
    v.vec[3] << "]";
}

int R4_vector::read(FILE* fp)
{
  if (fscanf(fp, "[%lf, %lf, %lf, %lf]", vec[0], vec[1], vec[2], vec[3]) != 4)
    return 0; 

  return 1; 
}

bool close(const R4_vector& a, const R4_vector& b, double eps)
{
  int i; 
  // Assumes both normalized st. first cpt. == 1.0. 
  for (i=1; i<4; i++) {
    if (fabs(a[i] - b[i]) > eps) return false; 
  }
  return true; 
}


// Gives a point on the light cone. 
O31_vector::O31_vector(const Complex& z)
{
  double d = (z.real * z.real + z.imag * z.imag + 1.0)/2.0; 
  vec[0] = 1.0;
  vec[1] = 1.0 - 1.0/d;
  vec[2] = z.real/d;
  vec[3] = z.imag/d;
}

O31_vector::operator Complex () const
{
  return Complex(vec[2],vec[3])/(vec[0] - vec[1]);
}

O31_matrix O31_x_trans(const Complex& tr)
{
  O31_matrix trans = O31_matrix(1); // initialize to identity

  trans.mx[0][0] = cosh(tr.real); 
    trans.mx[0][1] = sinh(tr.real); 
  trans.mx[1][0] = sinh(tr.real); 
    trans.mx[1][1] = cosh(tr.real); 

  trans.mx[2][2] = cos(tr.imag); 
    trans.mx[2][3] = -sin(tr.imag); 
  trans.mx[3][2] = sin(tr.imag); 
    trans.mx[3][3] = cos(tr.imag); 

  return trans;
}

O31_matrix O31_y_trans(const Complex& tr)
{
  O31_matrix trans = O31_matrix(1); // initialize to identity

  trans.mx[0][0] = cosh(tr.real); 
    trans.mx[0][2] = sinh(tr.real); 
  trans.mx[2][0] = sinh(tr.real); 
    trans.mx[2][2] = cosh(tr.real); 

  trans.mx[1][1] = cos(tr.imag); 
    trans.mx[1][3] = sin(tr.imag); 
  trans.mx[3][1] = -sin(tr.imag); 
    trans.mx[3][3] = cos(tr.imag); 

  return trans;
}

O31_matrix O31_z_trans(const Complex& tr)
{
  O31_matrix trans = O31_matrix(1); // initialize to identity

  trans.mx[0][0] = cosh(tr.real); 
    trans.mx[0][3] = sinh(tr.real); 
  trans.mx[3][0] = sinh(tr.real); 
    trans.mx[3][3] = cosh(tr.real); 

  trans.mx[1][1] = cos(tr.imag); 
    trans.mx[1][2] = -sin(tr.imag); 
  trans.mx[2][1] = sin(tr.imag); 
    trans.mx[2][2] = cos(tr.imag); 

  return trans;
}

double operator * (const O31_vector& a, const O31_vector& b) // minkowski dot
{
  int i; 
  double dot; 

  dot = -(a[0] * b[0]);
  for (i=1; i<4; i++)
    dot += a[i] * b[i];
  return dot; 
}

bool hyp_close(O31_vector const& a, O31_vector const& b, double eps)
{
  return (1.0 - a*b) < eps*eps/2.0; 
}

void GL4R_matrix::invert()
{
  if (gl4R_invert(mx, mx)!=func_OK)
    fprintf(stderr, "inversion of GL4R_matrix failed\n"); 
}


O31_vector O31_matrix::image_of_origin() const
{
  O31_vector result;
  int i; 
  for (i=0; i<4; i++) result[i] = mx[i][0]; 
  return result; 
}

bool GL4R_matrix::is_identity(double eps) const
{
  int i, j; 
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      if (fabs(mx[i][j] - (i==j ? 1.0:0.0)) > eps) return false; 
    }
  }
  return true; 
}

bool close(GL4R_matrix const& a, GL4R_matrix const& b, double eps)
{
  int i, j; 
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      if (fabs(a(i,j) - b(i,j)) > eps) 
	return false; 
    }
  }
  return true; 
}

const O31_vector O31_origin = O31_vector(1.,0.,0.,0.); 

void Moebius_array_to_O31_array(const MoebiusTransformation arrayA[], O31_matrix arrayB[], int num_matrices)
{
  int i; 
  for (i=0; i<num_matrices; i++) arrayB[i] = O31_matrix(arrayA[i]); 
}

void O31_array_to_Moebius_array(const O31_matrix arrayB[], MoebiusTransformation arrayA[], int num_matrices)
{
  int i; 
  for (i=0; i<num_matrices; i++) arrayA[i] = O31_matrix(arrayB[i]); 
}

void snappea_print_o31(const O31_matrix& mx, FILE* fp)
{
  int i, j; 

  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      fprintf(fp, "%.16f ", mx(i,j)); 
    }
    fprintf(fp, "\n"); 
  }
  fprintf(fp, "\n"); 
}

