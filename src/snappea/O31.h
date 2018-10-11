// -*-C++-*-
#ifndef _O31_
#define _O31_

#include "matrix_conversion.h"
// This includes o31_matrices and Moebius_transformations as well as 
// the conversion functions for going between them. 

// Wrapper for Jeff`s o31_matrices.[ch].

class R4_vector;

class GL4R_matrix {
protected:
  double mx[4][4];

public:
  void copy(const GL4R_matrix& m); 

  GL4R_matrix() {}
  GL4R_matrix(int); // constructs identity matrix
  GL4R_matrix(const GL4R_matrix& m) { copy(m); }

  double operator () (int r, int c) const { return mx[r][c]; }
  double& operator () (int r, int c) { return mx[r][c]; }

  O31Vector const* M() const { return mx; }
  O31Vector* M() { return mx; }

  double distance_from_identity() const;
  bool is_identity(double eps) const; 

  void O31_GramSchmidt() { ::o31_GramSchmidt(mx); }
  void invert();

  GL4R_matrix& operator = (const GL4R_matrix& rhs)
    { copy(rhs); return *this; }
  void operator *= (const GL4R_matrix& m) 
    { o31_product(mx, m.mx, mx); }
  void left_mul(const GL4R_matrix& m) // left multiplication: m * *this
    { o31_product(m.mx, mx, mx); }

  void print(FILE* fp = stdout) const; 
  int read(FILE* fp = stdin); 

  friend GL4R_matrix operator * (const GL4R_matrix& a, const GL4R_matrix& b)
  { GL4R_matrix p; o31_product(a.mx, b.mx, p.mx); return p; }
  friend R4_vector operator * (const GL4R_matrix& m, const R4_vector& v);

  friend GL4R_matrix transpose(const GL4R_matrix& m); 
  friend double determinant(GL4R_matrix const& m)
    { return gl4R_determinant(m.mx); }
  friend double trace(GL4R_matrix const& m)
    { return o31_trace(m.mx); }
  friend GL4R_matrix inverse(const GL4R_matrix& m); 
  friend bool compute_inverse(const GL4R_matrix& m, GL4R_matrix& m_inverse)
    { return gl4R_invert(m.mx, m_inverse.mx) == func_OK; }

  friend std::ostream& operator << (std::ostream& out, const GL4R_matrix& m);
  friend bool close(GL4R_matrix const& a, GL4R_matrix const& b, double eps);
};

class R4_vector {
protected:
  double vec[4]; 

public:
  void copy(const R4_vector& v); 

  R4_vector() {}
  R4_vector(double x, double y, double z, double w);
  R4_vector(const R4_vector& v) { copy(v); }

  R4_vector& operator = (const R4_vector& rhs)
    { copy(rhs); return *this; }
  void operator += (const R4_vector& rhs); 
  void operator -= (const R4_vector& rhs); 
  void operator *= (double r); 
  void operator /= (double r); 
  void negate(); 
  void left_mul(const GL4R_matrix& m) // left multiplication: m * *this
    { o31_matrix_times_vector(m.M(), vec, vec); }

  // Unsafe conversion operator for compatibility with old code.
  operator double* () const
    { return (double*)vec; }

  friend double operator * (const R4_vector& a, const R4_vector& b); // dot

  double norm() const 
    { return sqrt(*this * *this); }
  double norm_squared() const
    { return *this * *this; }

  void normalize() // result has norm == 1
    { *this /= norm(); } 

  double operator [] (int i) const { return vec[i]; }
  double& operator [] (int i) { return vec[i]; }

  friend std::ostream& operator << (std::ostream& out, const R4_vector& v);

  void print(FILE* fp = stdout) const; 
  int read(FILE* fp = stdin); 

  friend R4_vector operator * (const GL4R_matrix& m, const R4_vector& v)
  { R4_vector mv; o31_matrix_times_vector(m.mx, v.vec, mv.vec); return mv; }

  friend bool close(const R4_vector& a, const R4_vector& b, double eps);
};

inline R4_vector operator - (const R4_vector& a)
{ R4_vector neg = a; neg.negate(); return neg; }

inline R4_vector operator + (const R4_vector& a, const R4_vector& b)
{ R4_vector sum = a; sum += b; return sum; }

inline R4_vector operator - (const R4_vector& a, const R4_vector& b)
{ R4_vector diff = a; diff -= b; return diff; }

inline R4_vector operator * (double r, const R4_vector& v)
{ R4_vector prod = v; prod *= r; return prod; }

inline R4_vector operator / (const R4_vector& v, double r)
{ R4_vector quot = v; quot /= r; return quot; }

class O31_vector; 

class O31_matrix : public GL4R_matrix {
public:
  O31_matrix() {}
  O31_matrix(int) : GL4R_matrix(1) {} // constructs identity matrix
  O31_matrix(const O31Matrix m) { o31_copy(mx, m); }
  O31_matrix(const O31_matrix& m) : GL4R_matrix(m) {}
  O31_matrix(const MoebiusTransformation& m)
    { Moebius_to_O31(&(MoebiusTransformation&)m, mx); }

  O31_vector row(int i) const;
  O31_vector col(int i) const;

  friend O31_matrix operator * (const O31_matrix& a, const O31_matrix& b)
  { O31_matrix p; o31_product(a.mx, b.mx, p.mx); return p; }
  friend O31_vector operator * (const O31_matrix& m, const O31_vector& v);

  O31_matrix& operator = (const O31_matrix& rhs)
    { copy(rhs); return *this; }
  void operator *= (const O31_matrix& m) 
    { o31_product(mx, m.mx, mx); }
  void left_mul(const O31_matrix& m) // left multiplication: m * *this
    { o31_product(m.mx, mx, mx); }

  friend O31_matrix inverse(const O31_matrix& m)
    { O31_matrix m_inverse; o31_invert(m.mx, m_inverse.mx); return m_inverse; }
  void invert()
    { o31_invert(mx, mx); }
  double deviation() const // checks that matrix is actually in O(3,1)
    { return o31_deviation(mx); }

  void conjugate_by(const O31_matrix& m) 
    { o31_conjugate(mx, m.mx, mx); }

  O31_vector image_of_origin() const; 
  double height() const { return mx[0][0]; }

  operator MoebiusTransformation() const
    { MoebiusTransformation m; O31_to_Moebius(mx, &m); return m; }

  friend O31_matrix O31_x_trans(const Complex& tr); 
  friend O31_matrix O31_y_trans(const Complex& tr); 
  friend O31_matrix O31_z_trans(const Complex& tr); 
};

O31_matrix O31_x_trans(const Complex& tr); 
O31_matrix O31_y_trans(const Complex& tr); 
O31_matrix O31_z_trans(const Complex& tr); 

class O31_vector : public R4_vector {
public:
  O31_vector() {}
  O31_vector(const O31Vector v) { o31_copy_vector(vec, v); }
  O31_vector(double x, double y, double z, double w) : R4_vector(x,y,z,w) {}
  O31_vector(R4_vector const& v) : R4_vector(v) {}
  O31_vector(Complex const& z);

  O31_vector& operator = (const O31_vector& rhs)
    { copy(rhs); return *this; }

  operator Complex () const; 

  friend double operator * (const O31_vector& a, const O31_vector& b); // minkowski dot

  double norm() const 
    { return sqrt(fabs(*this * *this)); }
  double norm_squared() const
    { return *this * *this; }

  void normalize() // result has norm_squared == +/-1
    { *this /= norm(); } 

  friend O31_vector operator * (const O31_matrix& m, const O31_vector& v)
  { O31_vector mv; o31_matrix_times_vector(m.mx, v.vec, mv.vec); return mv; }
};

inline O31_vector operator - (const O31_vector& a)
{ O31_vector neg = a; neg.negate(); return neg; }

inline O31_vector operator + (const O31_vector& a, const O31_vector& b)
{ O31_vector sum = a; sum += b; return sum; }

inline O31_vector operator - (const O31_vector& a, const O31_vector& b)
{ O31_vector diff = a; diff -= b; return diff; }

inline O31_vector operator * (double r, const O31_vector& v)
{ O31_vector prod = v; prod *= r; return prod; }

inline O31_vector operator / (const O31_vector& v, double r)
{ O31_vector quot = v; quot /= r; return quot; }

// both must be normalized. 
bool hyp_close(O31_vector const& a, O31_vector const& b, double eps); 

extern const O31_vector O31_origin; 

/*
 *  Convert arrays of matrices back and forth between SL(2,C) and O(3,1).
 */
void Moebius_array_to_O31_array(const MoebiusTransformation arrayA[], O31_matrix arrayB[], int num_matrices);
void O31_array_to_Moebius_array(const O31_matrix arrayB[], MoebiusTransformation arrayA[], int num_matrices);

void snappea_print_o31(const O31_matrix& mx, FILE* fp);

#endif
