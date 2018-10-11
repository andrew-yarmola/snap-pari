#ifndef _p_VS_hh_
#define _p_VS_hh_

#include "pariwrap.hh"
#include <iostream>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef PARI_2_2_OR_LATER
#include "pari_oldnames.hh"
#endif

#undef zero
// Avoid clash with pari normalize func. 
#define normalize my_normalize

struct p_VS {
  typedef pari V;
  typedef pari M;
  typedef pari SC;

  static const pari& zero; 
  static const pari& half; 
  static const pari& one;
 
  static double eq_EPS; 
  static double ne_EPS;

  static pari scalar(int n)
  { return integer(n); }
  static pari vector(const char* s)
  { return p_lisexpr(s); }
  static pari vector(int dim)
  { return rvector(dim); }
};

inline double size(pari const& p)
{ return gabs(p).double_value(); }

inline pari norm(pari const& p)
{ return sqrt(gnorml2(p)); }
inline int dim(pari const& v)
{ return length(v); }
inline void my_normalize(pari& p)
{ p /= norm(p); }

inline pari id_mat(int d)
{ return p_idmat(d); }
inline int rows(pari const& m)
{ return length(m[0]); }
inline int cols(pari const& m)
{ return length(m); }
inline pari row(pari const& m, int r)
{ return gtrans(m)[r]; }
inline pari col(pari const& m, int c)
{ return m[c]; }
inline pari inverse(pari const& m, bool* inv)
{ *inv = true; return ginv(m); }

pari dotprod(pari const& a, pari const& b);
void dehomogenized_copy(pari& c, pari const& v);

std::ostream& operator << (std::ostream& out, pari const& p);

int find_pivot(pari const& c, pari const& eps);
pari li_col_mask(pari const& m0, pari const& eps);

#endif
