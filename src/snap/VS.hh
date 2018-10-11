#ifndef _VS_hh_
#define _VS_hh_

#include "rc_matrix.hh"

#undef zero

struct VS {
  typedef n_vector  V;
  typedef rc_matrix M;
  typedef double   SC;

  static const SC zero; 
  static const SC half; 
  static const SC one; 
  static const double eq_EPS; 
  static const double ne_EPS;

  static double scalar(int n)
  { return n; }
  static n_vector vector(int dim)
  { return n_vector(dim); }
  static n_vector vector(const char* s)
  { return n_vector(s); }
};

inline void normalize(n_vector& V)
{ V.normalize(); }
inline double norm(n_vector const& V)
{ return V.norm(); }
inline void dehomogenized_copy(n_vector& c, n_vector const& v)
{ c.dehomogenized_copy(v); }
inline double size(double const& x)
{ return fabs(x); }
inline rc_matrix id_mat(int d)
{ return rc_matrix(d, 1.0); }
inline int dim(n_vector const& v)
{ return v.dim; }
inline int rows(rc_matrix const& m)
{ return m.r; }
inline int cols(rc_matrix const& m)
{ return m.c; }
inline double dotprod(n_vector const& a, n_vector const& b)
{ return a*b; }
inline n_vector row(rc_matrix const& m, int r)
{ return m.row(r); }
inline n_vector col(rc_matrix const& m, int c)
{ return m.col(c); }

#endif
