#include "p_VS.hh"
#include <iomanip>

using std::ostream;
using std::setw; 

const pari& p_VS::zero(pZERO);
const pari& p_VS::half(pHALF);
const pari& p_VS::one(pONE);

double p_VS::eq_EPS=1e-11;
double p_VS::ne_EPS=1e-9;

void dehomogenized_copy(pari& c, pari const& v)
{
  int dim = length(c); 
  if (dim != length(v)-1) {
    std::cout << "Dimensions wrong for dehomogenized_copy\n";
    return; 
  }
  pari div = v[0]; 
  int i; 
  for (i=0; i<dim; i++) c[i]=v[i+1]/div;
}

ostream& operator << (ostream& out, pari const& p)
{
  char buf[50]; 
  int i, n, j, m; 
  switch (p.type()) {
  case t_INT:
    out << p.int_value(); 
    break; 
  case t_REAL:
    if (abs(p) < pEPSILON) {
      out << "  0.000000";
      break; 
    }
    if (glog(abs(p)) > integer(230)) {
      if (p>pZERO) out << "       BIG";
      else out << "      -BIG";
      break;
    }
    sprintf(buf, "% 10f", p.double_value());
    out << buf;
    break; 
  case t_COMPLEX:
    if (abs(p) < pEPSILON) {
      out << "  0.000000 +0.000000*i";
      break;
    }
    if (glog(abs(p)) > integer(230)) {
      out << "              Infinity";
      break;
    }
    sprintf(buf, "% 10f%+10f*i", p[0].double_value(), p[1].double_value());
    out << buf;
    break;
  case t_VEC:
  case t_COL:
    n = length(p);
    out << '[';
    for (i=0;i<n;++i) {
      if (i) out << ',';
      out << setw(10) << p[i];
    }
    out << ']';
    break;
  case t_MAT:
    n = length(p); // cols
    if (!n) { out << "Mat([])"; break; }
    m = length(p[0]); // rows
    for (j=0; j<m; ++j) {
      out << (j ? ' ' : '[');
      for (i=0; i<n; ++i) {
	if (i) out << ',';
	out << setw(10) << p[i][j];
      }
      out << ((j<m-1) ? ';':']') << std::endl;
    }
    break;
  default:
    out << "pari(t=" << p.type() << ")";
    break;
  }
  return out; 
}

pari dotprod(pari const& a, pari const& b)
{ 
  if (a.type()==t_VEC) {
    if (b.type()==t_VEC)
      return a*gtrans(b); 
    return a*b;
  } 
  if (b.type()==t_VEC)
    return gtrans(a)*gtrans(b); 
  return gtrans(a)*b;
}

int find_pivot(pari const& c, pari const& eps)
{
  int i, p, r = length(c); 
  pari max_abs; // initialized to 0. 
  for (i=0; i<r; ++i) {
    if (abs(c[i]) > max_abs) {
      max_abs = abs(c[i]);
      p = i; 
    }
  }
  if (max_abs < eps) p = -1; 
  return p; 
}

pari li_col_mask(pari const& m0, pari const& eps)
{
  pari m = gcopy(m0); 
  int i, j, p, n = length(m); // num cols
  pari mask; // initialized to integer(0)
 
  for (i=0; i<n; ++i) {
    p = find_pivot(m[i],eps);
    if (p<0) continue; 

    // record this row in mask. 
    mask += pow(pTWO, integer(i)); 

    for (j=i+1; j<n; ++j)
      m[j] -= (m[j][p]/m[i][p])*m[i];
  }
  return mask;
}

