/*
** Copyright (C) 2003 Oliver A. Goodman <oag@ms.unimelb.edu.au>
**  
** This file is part of Snap.
** 
** Snap is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software 
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cfloat>
#include <iomanip>
#include "eqs_interval.hh"
#include "adaptive_curve.hh"
#include "printable.hh"
#include <cstring>

using std::strcmp;
using std::swap; 
using std::reverse;
using std::ios;

// put this in bq_curve
double eqs_interval::epsilon = 1.e-6;
int eqs_interval::talk = 0; 

#define EPS eqs_interval::epsilon
#define TALK eqs_interval::talk

const Mark off_end(-10); 

double watch_for_x = 2.0; 

static bool small(double a)
{
  return fabs(a) < EPS;
}

static bool close(double a, double b)
{
  return fabs(b-a) < EPS; 
}

static int signof(double x)
{
  return (x > EPS ? 1:(x < -EPS ? -1:0));
}

static bool same(Complex const& a, Complex const& b)
{ 
  return complex_close(a,b,EPS); 
}

static bool close_mats(const R_matrix<3>& a, const R_matrix<3>& b)
{
  int i, j; 
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      if (!close(a[i][j], b[i][j])) break;
    }
    if (j<3) break; 
  }
  return i==3; 
}

static double evaluate(const R_matrix<3>& R, double x, double y)
{
  R_vector<3> x_pow, y_pow; 

#ifdef __DECCXX
  x_pow[0] = 1.0; 
  if (x == DBL_MAX) {
    x_pow[1] = 0.;
    x_pow[2] = 1e50;
  } else {
    x_pow[1] = x; 
    x_pow[2] = x*x; 
  }
  y_pow[0] = 1.0; 
  if (y == DBL_MAX) {
    y_pow[1] = 0.;
    y_pow[2] = 1e50;
  } else {
    y_pow[1] = y; 
    y_pow[2] = y*y; 
  }
#else
  x_pow[0] = 1.0; 
  x_pow[1] = x; 
  x_pow[2] = x*x; 
  y_pow[0] = 1.0; 
  y_pow[1] = y; 
  y_pow[2] = y*y; 
#endif 
  return y_pow * (R * x_pow); 
}

double evaluate(const R_matrix<3>& R, Complex const& z)
{
  return evaluate(R, z.real, z.imag); 
}

Complex evaluateD(const R_matrix<3>& R, Complex const& z)
{
  Complex D;

  R_vector<3> x_pow, y_pow; 

  x_pow[0] = 0.;
  x_pow[1] = 1.;
  x_pow[2] = 2*z[0];
  y_pow[0] = 1.;
  y_pow[1] = z[1];
  y_pow[2] = z[1]*z[1];

  D[0] = y_pow * (R * x_pow); 

  x_pow[0] = 1.;
  x_pow[1] = z[0];
  x_pow[2] = z[0]*z[0];
  y_pow[0] = 0.;
  y_pow[1] = 1.;
  y_pow[2] = 2*z[1];

  D[1] = y_pow * (R * x_pow); 

  return D; 
}

static void get_derivative(int xy, R_matrix<3> const& M, R_matrix<3>& DM)
{
  int i, j;
  for (i=0; i<3; i++) {
    for (j=1; j<3; j++) {
      if (xy==0) 
	DM[i][j-1] = j*M[i][j];
      else
	DM[j-1][i] = j*M[j][i];
    }
    if (xy==0) DM[i][2] = 0.;
    else DM[2][i] = 0.;
  }
}

static double bound(const R_matrix<3>& M)
{
  double val, bound = 0.0;
  int i, j;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      val = fabs(M[i][j]);
      if (val > bound) bound = val; 
    }
  }
  return bound; 
}

inline static Complex pnt(int xy, double x, double y)
{
  return xy ? Complex(y,x) : Complex(x,y);
}

static vector<double> critical_values(const vector<R_poly>& f, R_poly& f_disc,
				      double lo, double hi)
{
  // Get the roots of the discriminant of f as poly in y. 
  f_disc = f[1]*f[1] - 4.0*f[2]*f[0];

  vector<double> rts = roots(f_disc, lo+EPS, hi-EPS);
  rts.insert(rts.begin(), lo);
  rts.push_back(hi); 

  return rts; 
}

static R_matrix<3> R_mat(vector<R_poly> const& f)
{
  R_matrix<3> R(0.); 
  int i, j;
  for (i=0; i<f.size(); ++i)
    for (j=0; j<f[i].size(); ++j)
      R(i,j) = f[i][j];
  return R;
}

bool bq_curve::factorize(vector<bq_curve>& factors) const
{
  factors.resize(0); 

  vector<R_poly> g(f(0)); // take a copy, we'll want to modify

  if (g.size()==0) {
    factors.push_back(*this);
    return false; // zero poly is irreducible! 
  }

  // first check for vertical lines which are common factors 
  // of the polys in x. 

  int i; 
  R_poly xcommon = g[0];
  for (i=1; i<g.size() && degree(xcommon) > 0; ++i)
    xcommon = hcf(g[i],xcommon,EPS); 

  R_poly remainder; // to keep ediv happy

  if (degree(xcommon) > 0) {

    vector<double> rts = roots(xcommon);

    R_matrix<3> R(0.); 
    R(0,1) = 1.;
    for (i=0; i<rts.size(); ++i) {
      R(0,0) = -rts[i];
      factors.push_back(bq_curve(R,mark(),arg()));
    }

    for (i=0; i<g.size(); ++i) {
      g[i] = ediv(g[i],xcommon,remainder,EPS);
    }
  }

  // deal with linear in y: irreducible unless we already
  // had a common factor in x. 

  if (g.size() < 3) {
    factors.push_back(bq_curve(R_mat(g),mark(),arg()));
    return factors.size() > 1;
  }

  R_poly disc = sqrt(g[1]*g[1] - 4.0*g[2]*g[0], EPS);

  if (!disc.size()) {
    factors.push_back(bq_curve(R_mat(g),mark(),arg()));
    return factors.size() > 1;
  }

  vector<R_poly> F(2);
  R_poly branch, cf;
  int k; 
  for (k=0; k<2; k++) {

    // Curves are y = (-g[1] +/- rD)/(2*g[2]), which are eqivalent to 
    // (2*g[2]) y + (g[1] -/+ rD) = 0, once we remove common factors from 
    // the coefficients. 

    branch = k ? (g[1] - disc) : (g[1] + disc);
    cf = hcf(branch, g[2], EPS);

    F[0] = ediv(branch/2., cf, remainder, EPS);
    F[1] = ediv(g[2], cf, remainder, EPS); 

    factors.push_back(bq_curve(R_mat(F),mark(),arg()));
  }

  return true;
}


static bool separate_curve_components(vector<R_poly> const& f, R_poly const& D, vector<R_poly> F[2])
{
  // D should already be chopped.
  R_poly rD = sqrt(D, EPS); 

  if (!rD.size()) 
    return false; // Curve is irreducible so nothing to be done. 

  R_poly cf, rem0, rem1, branch; 

  int i; 
  for (i=0; i<2; i++) {

    F[i].resize(2); 

    // Curves are y = (-f[1] +/- rD)/(2*f[2]), which are eqivalent to 
    // (2*f[2]) y + (f[1] -/+ rD) = 0, once we remove common factors from 
    // the coefficients. We return the coefficients in F[i][1] and F[i][0]. 

    branch = i ? (f[1] - rD) : (f[1] + rD);
    cf = hcf(branch, f[2], EPS);

    F[i][0] = ediv(branch, cf, rem0, EPS);
    F[i][1] = ediv(2.0*f[2], cf, rem1, EPS); 

    if (rem0.size() || rem1.size()) {
      cout << "Problem separating curve components.\n";
      return false; 
    }
  }
  return true; 
}

static bool has_line(vector<R_poly> const& f, vector<double> const& critical)
{
  // Check for vertical lines in zero set. 
  int i, j; 
  for (i=0; i<critical.size(); i++) {
    for (j=0; j<3; j++) {
      if (fabs(evaluate(f[j], critical[i])) > EPS) break; 
    }
    if (j==3) return true; 
  }
  return false; 
}

static void copy_to_polys(vector<R_poly>& f, const R_matrix<3>& R0, int xy)
{
  int i, j, deg = -1; 

  f.resize(3);

  // Copy R0 into polynomials. 
  for (i=0; i<3; i++) {
    f[i] = R_poly(3);
    for (j=0; j<3; j++) {
      f[i][j] = xy ? R0[j][i] : R0[i][j];
      if (deg < i && !small(f[i][j])) deg = i; 
    }
  }
  if (deg < 2) f.resize(deg+1);
}

R_poly get_resultant(vector<R_poly> const& f, vector<R_poly> const& g)
{
  if (f.size()==2) {
    if (g.size()==2)
      return f[0]*g[1] - f[1]*g[0];
    return sqr(f[1])*g[0] - f[0]*f[1]*g[1] + sqr(f[0])*g[2]; 
  }
  if (f.size()==3) {
    if (g.size()==2)
      return sqr(g[1])*f[0] - g[0]*g[1]*f[1] + sqr(g[0])*f[2]; 
    return sqr(f[0]*g[2] - f[2]*g[0]) + 
      (f[0]*g[1] - f[1]*g[0]) * (f[2]*g[1] - f[1]*g[2]); 
  }
  return R_poly(1,1.);
}

R_poly get_resultant(R_matrix<3> const& F, R_matrix<3> const& G, int xy=0)
{
  vector<R_poly> f(3), g(3);

  copy_to_polys(f,F,xy);
  copy_to_polys(g,G,xy);

  return get_resultant(f, g); 
}

// BQ_CURVE STUFF

void bq_curve::scale(double s)
{
  vector<double> pows(5);
  int i, j; 
  pows[0] = 1.0;
  for (i=1; i<5; i++) pows[i] = s*pows[i-1];
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      M[i][j] *= pows[i+j];
    }
  }
}

ostream& operator << (ostream& out, bq_curve const& m)
{
  out << m.M;
  if (fabs(m.arg()) > EPS) out << " arg " << m.arg(); 
  out << ' ' << m.mark_num; 
  return out; 
}

bool bq_curve::read(FILE* fp, bool get_marking)
{
  if (!M.read(fp)) return false;
  if (!get_marking) return true; 

  char buf[256];

  // check for optional argument
  if (strcmp("arg", buf)==0) {
    if (fscanf(fp, " %lf ", &_arg)!= 1) return false; 
    if (fscanf(fp, " %s ", buf)!= 1) return false;
  } else {
    _arg = 0.;
  }

  // get the marking
  if (fscanf(fp, " %s ", buf)!= 1) return false; 
  return mark_num.set_to(buf); 
}

void bq_curve::x_reflect()
{
  int i; 
  for (i=0; i<3; i++) M[i][1] = -M[i][1]; 
  F[0].resize(0);
  F[1].resize(0);
}

// BQ_CURVE STUFF

vector<R_poly> const& bq_curve::f(int xy) const
{
  if (!F[xy].size()) {

    copy_to_polys(((bq_curve*)this)->F[xy],Mx(),xy); 
  }
  return F[xy]; 
}

double bq_curve::get_yd(int xy, double x, double& disc, int* deg) const
{
  // Get values of coeffs of y at this x value. 
  double F[3]; 
  int i;
  for (i=0; i<f(xy).size(); i++) 
    F[i] = evaluate(f(xy)[i], x); 
  --i; 
  while (i >= 0 && small(F[i])) --i;

  if (deg) *deg = i; 
  disc = 0.;  

  if (i <1) return DBL_MAX;

  if (i==1) return -F[0]/F[1];

  // Compute y value from formula for roots of a quadratic. 
  F[0] /= F[2];
  F[1] /= F[2];
  double d = F[1]*F[1]/4.0 - F[0]; 
  if (d > 0.) disc = sqrt(d); 
  else if (d < -EPS) disc = -1.; 

  return -F[1]/2;
}

double bq_curve::Y(int xy, int sign, double x) const
{
  double disc;
  double y = get_yd(xy, x, disc);
  return y + sign * disc; 
}

Complex bq_curve::at(int xy, int sign, double x) const
{
  return pnt(xy, x, Y(xy,sign,x));  
}

bool bq_curve::get_branch(int xy, Complex const& p, int& sign, int dir) const
{
  int num_pts; 
  double disc;
  double px = p[xy];
  double y = get_yd(xy, px, disc, &num_pts);
  double ydiff = p[1-xy] - y;

  // vertical line. 
  if (num_pts < 0) {
    sign = 0; 
    return true; 
  }

  // both branches off at infinity
  if (num_pts == 0) { 
    return false; 
  }

  if (TALK > 2) {
    ios::fmtflags old_fmt = cout.flags();
    cout.unsetf(ios::fixed); // use scientific notation. 
    cout << "    pt " << p;
    cout << " abs(abs(ydiff)-disc) " << fabs(fabs(ydiff) - disc) << endl; 
    cout.flags(old_fmt); 
  }

  // check if we're on the curve
  if (!small(fabs(ydiff) - disc)) {
    return false; 
  }

  if (dir && num_pts == 1) {

    // deal with a branch swap. 
    R_poly Df2 = derivative(f(xy)[2]); 

    // for px increasing, branch will be sign of b/da.
    // where y = root of ay^2 + by + c. 
    double da = evaluate(Df2,px) * dir;
    double b  = evaluate(f(xy)[1],px);

    sign = signof(da * b); 
    return true; 
  }

  sign = signof(ydiff); 
  return true; 
}

#if 0
bool bq_curve::resultant_vanishes(int xy, const R_matrix<3>& M) const
{
  vector<R_poly> g(3);
  // should this be copy to polys?
  int i, j;
  for (i=0; i<3; i++) {
    g[i] = R_poly(3); 
    for (j=0; j<3; j++) {
      g[i][j] = xy ? M[j][i] : M[i][j];
    }
  }
  R_poly resultant = get_resultant(f(xy), g); 
  return small(upper_bound(resultant)/(bound(M)*bound(Mx())));
}
#endif

static inline double taxi(Complex const& z)
{
  return fabs(z.real) + fabs(z.imag); 
}

static bool always_polish = false;

bool polish_point(Complex& z, R_matrix<3> const& A, R_matrix<3> const& B)
{
  Complex AB(evaluate(A,z),evaluate(B,z));

  double abs_val = taxi(AB); 

  if (abs_val > 1e-2) return false; // not near a solution. 

  if (100. * abs_val < EPS && !always_polish) {
    if (TALK > 1) cout << "point " << z << " error " << abs_val << endl; 
    return true; // already very accurate
  }

  // if (TALK) cout << "polishing " << z << " abs_val " << abs_val << endl; 

  Complex DAB[2];
  Complex Dz, nz, NAB;
  double det; 

  int guard=10;
  while (--guard > 0) {
 
    DAB[0] = evaluateD(A,z);
    DAB[1] = evaluateD(B,z);
    det = DAB[0][0]*DAB[1][1] - DAB[1][0]*DAB[0][1];

    if (fabs(det) > 1e-3) {

      // do newton's method
      Dz = Complex(( DAB[1][1]*AB[0]-DAB[0][1]*AB[1])/det,
		   (-DAB[1][0]*AB[0]+DAB[0][0]*AB[1])/det);

    } else if (taxi(DAB[0]) > 1e-3) {
      
      // do steepest descent 
      Dz = AB[0] * DAB[0]/complex_modulus_squared(DAB[0]);
      if (TALK > 1) cout << "polish_point doing steepest descent\n";
      
    } else {
      if (TALK) cout << "small derivatives in polish_point\n";
      return true; // derivatives vanish.. don't know what to do.
    }

    nz = z - Dz;
    NAB = Complex(evaluate(A,nz),evaluate(B,nz));

    if (taxi(NAB) > abs_val/2.) {
      if (TALK && guard==9) { 
	cout << "polish_point didn't help\n";
	cout << "at " << z << " went from " << AB << " to " << NAB << endl; 
      }
      break; 
    } else {
      z = nz; 
      AB[0] = NAB[0]; 
      AB[1] = NAB[1];
      abs_val = taxi(AB);
    }
  }

  if (TALK > 1) cout << "point " << z << " error " << abs_val << endl; 
  return true; 
}

vector<Complex> intersection_points(bq_curve const& A, bq_curve const& B, int xy, bool irred)
{
  vector<Complex> points; 

  if (A.f(xy).size()==1 || B.f(xy).size()==1) {
    xy = 1-xy;
  }

  R_poly resultant = get_resultant(A.f(xy), B.f(xy)); 

  int i;

  if (TALK > 2) {
    for (i=0; i<A.f(xy).size(); i++)
      cout << "A " << A.f(xy)[i] << endl;
    for (i=0; i<B.f(xy).size(); i++)
      cout << "B " << B.f(xy)[i] << endl;
    cout << "    resultant A,B " << resultant << endl;
  }


  if (small(upper_bound(resultant))) {
    if (irred) return points; 


    if (TALK > 1) cout << "  factorizing to determine intersections\n";

    vector<bq_curve> Afac,Bfac;
    A.factorize(Afac);
    B.factorize(Bfac);
    if (Afac.size() + Bfac.size() < 3) {
      if (TALK > 1) cout << "  factorization failed, trying swapping x and y\n";
      return intersection_points(A,B,1-xy,true);
    }

    if (TALK > 1) {
      cout << "A has " << Afac.size() << " factors\n";
      cout << "B has " << Bfac.size() << " factors\n";
    }

    int i,j,k;
    vector<Complex> these_points;
    for (i=0; i<Afac.size(); ++i) {
      for (j=0; j<Bfac.size(); ++j) {
	these_points = intersection_points(Afac[i],Bfac[j],xy,true);
	for (k=0; k<these_points.size(); ++k)
	  points.push_back(these_points[k]);
      }
    }

    return points;
  }

  vector<double> crossing;
  vector<vector<double> > adr;
  adr = all_deriv_roots(resultant, -1.-EPS, 1.+EPS);
  if (adr.size()==1) {
    crossing = adr[0]; 
  } else {
    crossing.resize(adr[0].size()+adr[1].size()); 
    merge(adr[0].begin(), adr[0].end(), 
	  adr[1].begin(), adr[1].end(), crossing.begin());
  }

  if (TALK > 1) 
    cout << "    crossing " << ("xy"[xy]) << "-values " << crossing << endl;

  // Get bound on B so we can normalize when evaluating. 
  double bsize = bound(B.Mx());

  // Check for crossing points at each x value. 
  double x, y, disc, yd;
  Complex z;
  int j, num_pts;
  for (j=0; j<crossing.size(); j++) {
    x = crossing[j]; 

    y = A.get_yd(xy,x,disc,&num_pts);

    if (disc < 0.) continue; // A is empty
    if (num_pts < 0) {
      // A has a vertical line here
      // so find intersections with B.
      y = B.get_yd(xy,x,disc,&num_pts);
      if (disc < 0. || num_pts < 1) continue;

    } 

#if 0
    // We'd like to eliminate tangential intersections
    // but we can't use the order of the resultant
    // to do it because this can also be even when
    // there are two transverse intersections at the same x value.

    else {
      if (order(resultant,x)%2==0) {

	// this is a tangency unless B has a vertical line
	B.get_yd(xy,x,yd,&num_pts);
	if (num_pts >= 0) continue;
      }
    }
#endif

    // cout << "x " << x << " y " << y << " disc " << disc << endl;

    if (small(disc)) {
      z = pnt(xy,x,y);
      if (polish_point(z, A.Mx(), B.Mx())) 
	points.push_back(z);
      continue;
    }

    z = pnt(xy,x,y-disc);
    if (polish_point(z, A.Mx(), B.Mx())) points.push_back(z);
    z = pnt(xy,x,y+disc);
    if (polish_point(z, A.Mx(), B.Mx())) points.push_back(z);
  }

  // cout << points << endl; 

  return points; 
}

// BRANCH STUFF

ostream& operator << (ostream& out, const branch& b)
{
  out << ("XY"[b.xy]) << ("-0+"[b.sign+1]);
  return out; 
}

// HELPERS FOR EQS_INTERVAL

class pt_w_branch : public Complex {
  int pos[2]; // branch for x or y.
public:
  pt_w_branch() {}
  pt_w_branch(Complex const& q) 
    : Complex(q) { pos[0]=pos[1]=2; }

  bool get_branch(bq_curve const& curve, int xy, int& sign) const;
};

struct critical_value {
  double x;
  int branch[2]; // branch from below/above x. 

  critical_value(double _x, int blo, int bhi)
    : x(_x) { branch[0] = blo; branch[1] = bhi; }
  bool on_this_branch(int br, int dir) const
  { return branch[dir > 0 ? 0:1] * br >= 0; }
  friend ostream& operator << (ostream& out, critical_value const& c);
};

class curve_study {
  bq_curve const& model;
  bool have[2];
  R_poly f_disc[2]; 
  double box; 
  vector<critical_value> bswap[2];
  vector<critical_value> crit[2];
  vector<pt_w_branch> endpoint; 

  void get(int xy);

public:
  curve_study(bq_curve const& m, double b) : model(m), box(b)
  { have[0] = have[1] = false; }

  enum nx_typ { none, endpt, branch_swap, critical, dbl_critical };

  void add_endpoint(Complex const& p);
  int find_endpoint(Complex const& p) const;
  void remove_endpoint(Complex const& p); 
  bool get_an_endpoint(Complex& p) const;
  int num_endpoints() const { return endpoint.size(); }
  void add_endpoints(vector<Complex> const& pts); 

  bq_curve const& curve() const { return model; }

  void add_boundary_endpoints(); 
  bool on_box(Complex const& p) const;
  bool in_box(double x) const { return fabs(x) < box+EPS; }
  void add_local_maxima(int xy); 
  int slope_sign(Complex const& p) const; 
  Complex tangent_direction(Complex const& p) const;
  bool choose_direction(Complex const& p, Complex& t) const;
  nx_typ next_p(branch const& b, int dir, Complex& p); 
  bool get_branch(Complex const& p, Complex const& t, branch& b) const;
  bool get_branch(int xy, Complex const& p, int& sign) const
  { return model.get_branch(xy,p,sign); }
};

// using curve_study::nx_typ; 

class branch_interval {
  bq_curve const& iv; 
  branch const& b; 
  double x;
  int dir; 
  bool stop_short;

  Complex np; 
  curve_study::nx_typ nt; 

public:

  branch_interval(bq_curve const& m, branch const& br, double _x, int d)
    : iv(m), b(br), x(_x), dir(d), stop_short(false),
      np(pnt(br.xy,1e9*d,0.)), nt(curve_study::none) {}

  void test_critical(critical_value const& v);
  void test_swap(critical_value const& v);
  bool test_point(pt_w_branch const& p);
  bool in_range(double px) const;
  void go_past(double px); 
  bool stopping_short() const { return stop_short; }
  curve_study::nx_typ get_point(Complex& p) const; 
};

// PT_W_BRANCH

bool pt_w_branch::get_branch(bq_curve const& curve, int xy, int& sign) const
{
  // cache the branch if we don't already have it. 

  if (pos[xy]!=2) { 
    sign = pos[xy]; 
    return true; 
  }

  if (!curve.get_branch(xy,(*this),sign)) 
    return false; 

  ((pt_w_branch*)this)->pos[xy] = sign; 
  return true; 
}

ostream& operator << (ostream& out, vector<pt_w_branch> const& v)
{
  int i, n=v.size();
  out << '[';
  for (i=0; i<n; i++) {
    if (i) out << ", ";
    out << v[i];
  }
  out << ']';
  return out; 
}

// CRITICAL_VALUE

ostream& operator << (ostream& out, critical_value const& c)
{ 
  return out << '(' << c.x << ',' << c.branch[0] << ',' << c.branch[1] << ')';
}

ostream& operator << (ostream& out, vector<critical_value> const& v)
{
  int i, n=v.size();
  out << '[';
  for (i=0; i<n; i++) {
    if (i) out << ", ";
    out << v[i];
  }
  out << ']';
  return out; 
}

// BRANCH_INTERVAL STUFF

bool branch_interval::in_range(double px) const
{
  if ((px - x)*dir < EPS) return false; // p is behind us.
  if ((px - np[b.xy])*dir > -EPS) return false; // p is ahead of np
  return true; 
}

void branch_interval::test_critical(critical_value const& v)
{
  if (!in_range(v.x)) return; 
  if (!v.on_this_branch(b.sign, dir)) return;

  np = pnt(b.xy, v.x, DBL_MAX);
  nt = curve_study::critical; 
  stop_short = true; 
  // cout << "  ..found critical " << v << endl;
}

void branch_interval::test_swap(critical_value const& v)
{
  if (!in_range(v.x)) return; 
  if (!v.on_this_branch(b.sign, dir)) return;

  np = iv.at(b.xy, 0, v.x);
  nt = curve_study::branch_swap; 
  stop_short = false; 
  if (TALK > 1) cout << "  ..found swap " << v << endl;
}

bool branch_interval::test_point(pt_w_branch const& p)
{
  if (!in_range(p[b.xy])) return false;
  int sign; 
  if (!p.get_branch(iv, b.xy, sign)) return false; 
  if (sign * b.sign < 0) return false; 

  np = p;
  nt = curve_study::endpt; 
  stop_short = false; 
  if (TALK > 1) cout << "  ..found endpoint " << p << endl;
  return true; 
}

void branch_interval::go_past(double px)
{
  // cout << "go past " << px << endl; 
  if (in_range(px)) x = px; 
}

curve_study::nx_typ branch_interval::get_point(Complex& p) const
{
  if (!stop_short) {
    p = np; 
    return nt;
  } 

  // cout << "getting pt between " << x << " and " << np[b.xy] << endl; 
  
  double px = 0.5*(x + np[b.xy]);
  p = iv.at(b.xy, b.sign, px);
  return nt; 
}

// CURVE_STUDY STUFF

static const char* typ_name[] = {
  "none", "endpt", "branch_swap", "critical", "dbl_critical"
};

bool curve_study::on_box(Complex const& p) const
{
  if (fabs(fabs(p[0]) - box) < EPS) return true;
  if (fabs(fabs(p[1]) - box) < EPS) return true;
  return false;
}

void curve_study::get(int xy)
{
  if (have[xy]) return; 

  vector<R_poly> const& f(model.f(xy)); 
  vector<double> C; 
  int i, sign, num_roots;

  // Examine the roots of f[2].
  R_poly Df2;
  double x, y, b, disc, tmp; 
  if (f.size()==3) {
    C = roots(f[2], -box, box);
    if (C.size()) Df2 = derivative(f[2]);
    for (i=0; i<C.size(); i++) {
      x = C[i];

      // get the (at most) one point here
      y = model.get_yd(xy,x,disc,&num_roots);

      // cout << "xy"[xy] << " = " << x << " num roots " << num_roots;
      // cout << " yx = " << y << endl; 

      // vertical line? ignore it and hope for the best. 
      if (num_roots < 0) continue; 

      if (num_roots ==0) { 
	// Both branches go to infinity.
	crit[xy].push_back(critical_value(x,0,0));
	continue; 
      }

      // One finite point here so b = f[1](x) != 0.
      // Branch swap if Df[2](x) != 0.
      b = evaluate(f[1],x);
      sign = signof(b * evaluate(Df2,x));

      if (sign) {
	// Looks like y=1/x U y=0. 
	if (in_box(y))
	  bswap[xy].push_back(critical_value(x,-sign,sign));
	crit[xy].push_back(critical_value(x,sign,-sign));
      } else {
	// Looks like y=1/x^2 U y=0.
	sign = -signof(f[2][2] * b);
	crit[xy].push_back(critical_value(x,sign,sign));
      }
    }
  } else if (f.size()==2) {

    // Roots of f[1].
    C = roots(f[1], -box, box);
    for (i=0; i<C.size(); i++) {
      // Some kind of vertical asymptote to be avoided. 
      crit[xy].push_back(critical_value(C[i],0,0));
    }
  }

  if (f.size() < 3) return; 
  // No other troublesome activity is possible in
  // this case. 

  // get the discriminant for this parameter direction
  f_disc[xy] = f[1]*f[1] - 4.0*f[2]*f[0];

  // critical values at zeros of discriminant
  C = roots(f_disc[xy], -box-EPS, box+EPS);

  // should really get roots to return multiplicities.
  R_poly ddisc  = derivative(f_disc[xy]); 
  int o; 

  for (i=0; i<C.size(); i++) {
    x = C[i];
    y = model.get_yd(xy,x,disc);

    if (!in_box(y)) continue;

    o = order(ddisc,x,EPS); 

    if (o==0) {
      // Curve becomes vertical here. 
      crit[xy].push_back(critical_value(x,0,0));

    } else if (o==1) {

      // Non-singular double point of the curve
      // will be a branch swap on both branches. 
      bswap[xy].push_back(critical_value(x,0,0));
    }
  }

  // cout << ("xy"[xy]) << "-critical points " << crit[xy] << endl; 
  // cout << ("xy"[xy]) << "-branch swaps    " << bswap[xy] << endl; 
  have[xy] = true; 
}

int curve_study::find_endpoint(Complex const& p) const
{
  int i, n=endpoint.size();
  for (i=0; i<n; i++) if (same(p,endpoint[i])) return i;
  return -1;
}

void curve_study::add_endpoint(Complex const& p)
{
  if (find_endpoint(p) >= 0) return; 
  endpoint.push_back(p); 
}

void curve_study::add_endpoints(vector<Complex> const& pts)
{
  int i, n=pts.size();
  for (i=0; i<n; ++i) endpoint.push_back(pts[i]);
}

void curve_study::remove_endpoint(Complex const& p)
{
  int i=find_endpoint(p);
  if (i<0) return; 
  endpoint.erase(endpoint.begin()+i);
} 

bool curve_study::get_an_endpoint(Complex& p) const
{
  if (!endpoint.size()) return false; 
  p = endpoint.front();
  return true; 
}

void curve_study::add_boundary_endpoints()
{
  double x, y, disc;
  bool one_branch; 
  int xy; 

  for (xy=1; xy>=0; --xy) {
    one_branch = model.f(xy).size() < 3; 
    
    for (x=box; x > -box-EPS; x -= 2*box) {
      y = model.get_yd(xy,x,disc);

      if (one_branch) { 
	if (in_box(y)) add_endpoint(pnt(xy,x,y));
	continue;
      }

      if (disc < EPS) continue;
      if (in_box(y-disc)) add_endpoint(pnt(xy,x,y-disc));
      if (in_box(y+disc)) add_endpoint(pnt(xy,x,y+disc));
    }
  }
}

void curve_study::add_local_maxima(int xy)
{
  get(xy); 

  if (!crit[xy].size()) return; // there are no critical points

  // maxima occur where derivative of discriminat is < 0.
  R_poly ddisc = derivative(f_disc[xy]); 

  double x, y; 
  int i, n = crit[xy].size(); 
  for (i=0; i<n; i++) {
    x = crit[xy][i].x;
    if (evaluate(ddisc,x) > -EPS) continue; 

    y = model.Y(xy, 0, x); 
    if (!in_box(y)) continue; 

    add_endpoint(pnt(xy,x,y));
  }

  // cout << "endpoints and maxima " << endl; 
  // cout << endpoint << endl; 
}

curve_study::nx_typ curve_study::next_p(branch const& b, int dir, Complex& np)
{
  get(b.xy); 

  double x=np[b.xy];

  branch_interval BI(model,b,x,dir);

  // Check for the first critical point..
  int i, n;
  n = crit[b.xy].size(); 
  for (i=0; i<n; i++)
    BI.test_critical(crit[b.xy][i]); 

  // or endpoint.. 
  n=endpoint.size();
  for (i=0; i<n; i++)
    BI.test_point(endpoint[i]); 

  // or branch swap.
  n = bswap[b.xy].size(); 
  for (i=0; i<n; i++)
    BI.test_swap(bswap[b.xy][i]); 

  // make sure we go past all y-critical pts.
  if (BI.stopping_short()) { 
    int yx = 1-b.xy;
    get(yx); 
    for (i=0; i<crit[yx].size(); i++)
      BI.go_past(model.Y(yx,0,crit[yx][i].x)); 
  }

  return BI.get_point(np); 
}

int curve_study::slope_sign(Complex const& p) const
{
  Complex t = tangent_direction(p);
  return signof(t.real*t.imag);
}

Complex curve_study::tangent_direction(Complex const& p) const
{
  // dM is an inward pointing normal
  Complex dM = evaluateD(model.Mx(),p);

  // rotate through -pi/2 to get tangent with positive side to the left.
  return Complex(dM.imag,-dM.real);
}

bool curve_study::choose_direction(Complex const& p, Complex& t) const
{
  int d[2] = {0,0};

  // get direction constraints for boundary points.
  if (fabs(fabs(p.real)-box) < EPS) 
    d[0] = -signof(p.real);
  if (fabs(fabs(p.imag)-box) < EPS) 
    d[1] = -signof(p.imag); 
  
  // for interior points (should have xy==1) go up. 
  if (d[0]==0 && d[1]==0) d[1]=1; 

  // check if we're clipping a corner
  if (d[0] * d[1] * t[0] * t[1] < 0.) return false; 

  // make sign of t[i] agree with any nonzero signs of d[i].
  if (d[0] * t[0] + d[1] * t[1] < 0.) t = -t; 

  return true; 
}  

bool curve_study::get_branch(Complex const& p, Complex const& t, branch& b) const
{
  // choose parameter
  b.xy = (fabs(t.imag) > fabs(t.real)) ? 1:0;
  return model.get_branch(b.xy, p, b.sign, signof(t[b.xy]));
}

// EQS_INTERVAL STUFF

ostream& operator << (ostream& out, const eqs_interval& i)
{
  out << "[";
  int j, n = i.BR.size();
  for (j=0; j<n; j++) {
    if (j) { 
      out << ','; 
      if (j%3==0) out << "\n ";
    }
    out << i.P[j] << ", " << i.BR[j];
  }
  if (i.P.size() > n) {
    out << ','; 
    if (j%3==0) out << "\n ";
    out << i.P[n];
  }
  out << "] " << i.mark();
  return out; 
}

double eqs_interval::xp(int e) const
{
  int xy;
  double _x=0.;
  if (e==0) return _x;
  int i, n=P.size();
  for (i=0; i<e; i++) {
    xy = BR[i].xy; 
    _x += fabs(P[(i+1)%n][xy] - P[i%n][xy]);
  }
  return _x;
}

double eqs_interval::x(int e) const
{
  return e==0 ? 0. : xp(BR.size());
}

void eqs_interval::reverse()
{
  std::reverse(P.begin(), P.end()); 
  std::reverse(BR.begin(), BR.end()); 
}

bool eqs_interval::read(FILE* fp)
{
  // get the matrix
  if (!C.read(fp)) return false; 

  char buf[256]; 

  // FIX THIS

  BR.resize(1);
  P.resize(2);

  // read whether parameter is x or y. 
  if (strcmp("xrange", buf)==0)
    BR[0].xy = 0; 
  else if (strcmp("yrange", buf)==0)
    BR[0].xy = 1; 
  else return false; 

  // get the parameter interval
  double x[2];
  if (fscanf(fp, " %lf %lf ", &x[0], &x[1]) != 2) return false; 

  // get the branch
  if (fscanf(fp, " %s ", buf)!= 1) return false;
  if (strcmp("upper", buf)==0)
    BR[0].sign = 1; 
  else if (strcmp("lower", buf)==0)
    BR[0].sign = -1; 
  else return false; 

  P[0] = at(0,x[0]); 
  P[1] = at(0,x[1]); 

  return true; 
}

void eqs_interval::write(ostream& out)
{
  out << C << endl; 
  out << *this; 
}

int eqs_interval::branch_num(double t, double& bx) const
{
  // check_me();

  if (t < 0.) {
    bx = x(0,0);
    return 0;
  }
  int i,n=BR.size();
  double bstep,bsize;
  for (i=0; i<n; ++i) {
    bstep = x(i,1) - x(i,0);
    bsize = fabs(bstep);
    if (t < bsize) break;
    t -= bsize;
  }
  if (i==n) {
    bx = x(n-1,1);
    return n-1;
  }
  bx = (bstep > 0) ? x(i,0) + t : x(i,0) - t;
  return i; 
}

int eqs_interval::branch_num(Complex const& p) const
{
  int sign[2] = {2,2};

  int xy;
  int b, n = BR.size();
  for (b=0; b<n; b++) {

    // parameter for this branch
    xy = BR[b].xy; 

    if (!in_range(b,p[xy])) continue; 

    // work out which branch for this parameter
    if (sign[xy]==2) 
      if (!C.get_branch(xy,p,sign[xy])) return -1; 

    // right branch? 
    if (BR[b].sign * sign[xy] >= 0) 
      return b; 
  }

  return -1; // never found p.
}

Complex eqs_interval::at(int b, double x) const
{
  return C.at(BR[b].xy, BR[b].sign, x);  
}


Complex eqs_interval::operator () (double t) const
{
  double x;
  int b = branch_num(t,x);
  return at(b,x);
}

void eqs_interval::close()
{
  // assume P.front()==P.back()

  P.pop_back(); 
  if (BR.front()==BR.back()) {
    P.erase(P.begin());
    BR.erase(BR.begin());
  }
}

bool eqs_interval::grow_interval(curve_study& CS, Complex& p, Complex& t)
{
  if (CS.find_endpoint(p) < 0) return false; 

  branch br;
  if (!CS.get_branch(p,t,br)) return false; 

  if (!P.size()) {

    // start a new interval. 
    P.push_back(p);
  } else if (!same(p,P.back())) {

    // continuing an interval from the wrong point!
    return false; 
  } else if (br==BR.back()) {

    // extend an existing interval
    P.pop_back();
    BR.pop_back();
  }
  
  int guard=10; 
  int dir = signof(t[br.xy]); 

  if (!dir) return false; 

  curve_study::nx_typ typ; 
  while (--guard > 0) {

    if (TALK > 1) 
      cout << "  direction " << dir << " branch " << br << endl; 

    // look along the curve for the next point. 
    typ = CS.next_p(br, dir, p);  

    if (TALK > 1) { 
      cout << "  next " << br << " " << p << endl; 
      cout << "  type " << typ_name[typ] << endl; 
    }

    if (typ==curve_study::none) return false; 
    // This can happen when we start at a point
    // on the boundary of the box where the curve
    // is tangent to the box. 

    // add a segment to this interval. 
    BR.push_back(br);
    P.push_back(p); 

    if (typ==curve_study::endpt) break;

    if (typ==curve_study::branch_swap) {

      br.sign = -br.sign; 

    } else { // typ == critical

      br.xy = 1 - br.xy; 
      if (CS.slope_sign(p) < 0) dir = -dir; 
  
      if (!CS.get_branch(br.xy, p, br.sign)) return false;
    }

  }

  // get new tangent and point it the way we were going.
  t = CS.tangent_direction(p);
  if (dir * t[br.xy] < 0.) t = -t; 

  check_me();

  return true;
}

bool eqs_interval::make_interval(curve_study& CS)
{
  Complex p0, p, t; 
  bool ok = false; 

  // get a point and starting direction
  while (CS.get_an_endpoint(p)) {

    t = CS.tangent_direction(p);
    if (!CS.choose_direction(p, t)) {

      // must have been on a corner of the box.
      CS.remove_endpoint(p); 
      continue; 
    }

    // cout << "pt, tgt = " << p << ' ' << t << endl;

    P.resize(0);
    BR.resize(0);

    p0 = p; 
    while ((ok = grow_interval(CS, p, t))) {
      
      if (same(p,p0)) {
	close();
	break;
      }
      if (CS.on_box(p)) break; 
      
      // we're heading through one of our local maxima
      CS.remove_endpoint(p);
    }

    CS.remove_endpoint(p0);
    CS.remove_endpoint(p);

    if (ok) break; 
  }

  if (ok) C = CS.curve();

  return ok;
}

void eqs_interval::x_reflect()
{
  C.x_reflect(); 
  int i; 
  for (i=0; i<P.size(); i++) P[i].real = -P[i].real; 
  for (i=0; i<BR.size(); i++)
    if (BR[i].xy==1) BR[i].sign = -BR[i].sign;
}

list<Complex> eqs_interval::get_polyline() const
{
  return adaptive_curve(*this, x(0), x(1), 0.2, 0.05);
}

bool eqs_interval::subdivide(Complex const& p, eqs_interval& a, eqs_interval& b, double eps) const
{
  int bn = branch_num(p);
  if (bn < 0) return false; 
  int nb = BR.size();

  if (complex_close(p, end(0), eps) || 
      complex_close(p, end(1), eps)) return false; 

  a.C = C;
  a.BR.resize(bn+1);
  a.P.resize(bn+2);
  int i;
  for (i=0; i<bn+1; i++) { 
    a.BR[i] = BR[i];
    a.P[i] = P[i];
  }
  a.P[bn+1] = p;

  b.C = C;
  b.BR.resize(nb-bn);
  b.P.resize(nb-bn+1);
  b.P[0] = p; 
  for (i=bn; i<nb; i++) { 
    b.BR[i-bn] = BR[i];
    b.P[i-bn+1] = P[i+1];
  }
  return true; 
}

eqs_interval::eqs_interval(eqs_interval const& i, int a, int b)
 : C(i.C)
{
  int j, n=i.P.size();

  if (a < b) {
    P.resize(b-a+1);
    BR.resize(b-a);
    for (j=a; j<b; ++j) {
      P[j-a] = i.P[j];
      BR[j-a] = i.BR[j];
    }
    P[b-a] = i.P[b];
  } else {
    // this interval is supposed to be closed
    P.resize(n-(a-b)+1);
    BR.resize(n-(a-b));
    for (j=a; j<n; ++j) {
      P[j-a] = i.P[j];
      BR[j-a] = i.BR[j];
    }
    for (j=0; j<b; ++j) {
      P[j+n-a] = i.P[j];
      BR[j+n-a] = i.BR[j];
    }
    P[n-(a-b)] = i.P[b];
  }

  if (!check_me()) {
    cout << "failed in eqs_interval(i," << a << ',' << b << ")" << endl;
    cout << i << endl; 
  }
}

void eqs_interval::insert_point(Complex const& p, int bn, vector<int>& indices)
{
  int pn=bn+1;
  bool adjust=false;

  // check if it's an endpoint of branch bn. 
  int bne=pn%P.size(); 
  if (same(p,P[bn])) {
    pn = bn; 
  } else if (same(p,P[bne])) {
    pn = bne;
  } else {
    vector<Complex>::iterator ci=P.begin()+bn+1;
    P.insert(ci,p);
    vector<branch>::iterator bi=BR.begin()+bn;
    BR.insert(bi,*bi);
    adjust = true; 
  }

  // cout << "indices " << indices << endl;
  // cout << "pn " << pn << " adjust=" << adjust << endl;

  // insert pn into indices keeping it sorted.
  int i;
  for (i=0; i<indices.size(); i++) if (indices[i] >= pn) break;
  if (i < indices.size() && indices[i]==pn && !adjust) return;

  indices.insert(indices.begin()+i,pn);

  // add 1 to indices above pn if a point was inserted.
  if (adjust) for (++i; i<indices.size(); ++i) ++indices[i];

  // cout << "new indices " << indices << endl;
}

int eqs_interval::subdivide(vector<Complex> const& S, vector<int>& indices, vector<Complex>* hits)
{
  int i, n=S.size(), bn, count=0;
  Complex p;
  indices.resize(0);

  if (!closed()) indices.push_back(0);

  for (i=0; i<n; i++) {
    p = S[i];
    if (!closed() && (same(p, end(0)) || 
		      same(p, end(1)))) {
      if (hits) hits->push_back(p);
      continue;
    }
    bn = branch_num(p);

    if (TALK > 1) 
      cout << "  point " << p << " branch num " << bn << endl; 

    if (bn < 0) continue;
    insert_point(p, bn, indices);
    if (hits) hits->push_back(p);
    ++count;
  }

  if (!closed()) 
    indices.push_back(P.size()-1);
  else if (indices.size()) 
    indices.push_back(indices.front()); 

  // cout << "possible intersections\n" << S << endl;
  // cout << "subdivided\n";
  // cout << *this << endl;
  // if (hits) cout << "hits " << *hits << endl; 
  // cout << "subdivision indices " << indices << endl;

  return count;
}

// This function might be better replaced with
// a slightly slower one which gives a midpoint, 
// guaranteed not to be close to either endpoint.

Complex eqs_interval::a_point(int a, int b) const
{
#if 0
  if (b <= a) b+= P.size(); 
  return (*this)((xp(a)+xp(b))/2);
#endif

  if (a < b) {
    if (a+1 < b) return P[(a+b)/2];
  } else {
    int n = P.size(); 
    if (n > a-b+1) return P[((a+b+n)/2)%n];
    // a == n-1
  }

  // return midpoint of branch a. 
  return at(a, (x(a,0)+x(a,1))/2.);
}

bool operator == (eqs_interval const& a, eqs_interval const& b)
{
  return same(a,b,EPS);
}

bool same(eqs_interval const& a, eqs_interval const& b, double eps)
{
  int i;
  for (i=0; i<2; i++) {
    if (complex_close(a.end(0), b.end(i),eps)) break; 
  }
  if (i==2 || !complex_close(a.end(1),b.end(1-i),eps)) return false; 
  
  double val = upper_bound(get_resultant(a.C.Mx(), b.C.Mx(), a.BR[0].xy)); 
  return val < eps * bound(a.C.Mx()) * bound(b.C.Mx()); 
}

bool eqs_interval::check_me() const
{
  if (P.size() > BR.size() + 1 || P.size() < BR.size()) {
    cout << "sizes of P and BR got out of sync\n";
    cout << "P.size() = " << P.size() << endl;
    cout << "BR.size() = " << BR.size() << endl;
    return false; 
  }

  int i; 
  branch br; 
  for (i=0; i<BR.size(); i++) {
    br = BR[i]; 
    if (br.xy != 0 && br.xy != 1) {
      cout << "invalid xy in BR[" << i << "] " << br.xy;
      cout << " " << mark() << endl; 
      return false; 
    }
    if (br.sign < -1 || br.sign > 1) {
      cout << "invalid sign in BR[" << i << "] " << br.sign << endl;
      return false; 
    }
  }
  return true; 
}

// INTERVAL_LIST STUFF

void interval_list::print(ostream& out) const
{
  out << LPSeq(L); 
}

void interval_list::print_sorted(ostream& out, double eps) const
{
  Complex prev; 
  list<eqs_interval>::const_iterator it; 
  for (it = L.begin(); it!=L.end(); it++) {
    if (it!=L.begin() && !complex_close(it->end(0), prev, eps)) {
      out << endl; 
    }
    out << *it << endl; 
    prev = it->end(1); 
  }
}

void interval_list::make_intervals(bq_curve const& curve)
{
  curve_study C(curve, 1.);

  C.add_boundary_endpoints(); 
  C.add_local_maxima(0); // 0 means for x. 

  eqs_interval i; 

  int ne;
  while ((ne = C.num_endpoints())) {
    L.push_back(i);
    if (!L.back().make_interval(C)) 
      L.pop_back();
    if (ne == C.num_endpoints()) {
      cerr << "make_intervals failed to use all endpoints\n";
      break;
    }
  }
}

inline static double dot(Complex const& a, Complex const& b)
{
  return a.real*b.real + a.imag*b.imag; 
}

bool interval_list::chop(bq_curve const& R)
{
  list<eqs_interval>::iterator ii=L.begin(), k;
  vector<int> indices;
  interval_list subs;
  vector<Complex> intersections;
  bool modified = false; 
  int i, n, a, b, s;
  bool keep, prev_kept;

  if (TALK > 1) {
    for (k=L.begin(); k!=L.end(); ++k) {
      Complex p = k->end(1); 
      cout << "  " << p << ' ' << evaluate(R.Mx(),p) << endl; 
    }
  }

  while (ii!=L.end()) {
    intersections = intersection_points(R, ii->C, ii->BR[0].xy);

    if (TALK > 1) {
      cout << "  intersections with " << ii->mark() << endl;
      cout << PSeq(intersections) << endl; 
    }

    if (ii->subdivide(intersections, indices)) {

      prev_kept = false; 
      n = indices.size(); 
      for (i=0; i<n; ++i) {

	if (i<n-1) 
	  keep = evaluate(R.Mx(), ii->a_point(indices[i],indices[i+1])) > 0.;
	else 
	  keep = false; 

	if (keep && !prev_kept) 
	  a = indices[i]; 
	if (!keep && prev_kept)
	  subs.L.push_back(eqs_interval(*ii,a,indices[i]));

	prev_kept = keep;
      }

      // replace *ii with subs.
      ii = L.erase(ii); 
      L.splice(ii,subs.L);

      modified = true; 
    } else {

      Complex p = ii->a_point(0,ii->P.size()-1);
      s = signof(evaluate(R.Mx(), p));

      if (s==0) // clipping along ii. 
	s = signof(dot(evaluateD(R.Mx(),p),evaluateD(ii->C.Mx(),p)));

      if (s > 0) {
	++ii;
      } else {
	ii=L.erase(ii);
	modified = true;
      } 
    }
  }

  for (ii=L.begin(); ii!=L.end(); ++ii) 
    if (!ii->check_me()) {
      cout << *ii << endl; 
      return -1; 
    }


  return modified; 
}

bool interval_list::end_connects(list<eqs_interval>::const_iterator it0) const
{
  list<eqs_interval>::const_iterator it(it0);
  Complex p=it->end(1);
  for (it=it0; it!=L.end(); ++it)
    if (same(p,it->end(0))) return true;
  for (it=L.begin(); it!=it0; ++it)
    if (same(p,it->end(0))) return true;
  return false;
}

void interval_list::tidy()
{
  list<eqs_interval>::iterator it, path, gap;
  if (!L.size()) return;

  // record start of a path
  Complex p, p0=L.front().end(0);

  for (path=L.begin(); path!=L.end();) {

    // find the first break
    p = path->end(1);
    for (it=path, ++it; it!=L.end(); ++it) {
      if (!same(it->end(0),p)) break;
      p = it->end(1);
    }

    // finished if this path continues to the end.
    if (it==L.end()) return;

    // finished this path if end matches start
    if (same(p0,p)) {
      path = it;
      p0 = path->end(0);
      continue;
    }

    // ok so we've found a gap. 
    gap = it;

    // look for matching end
    for (++it; it!=L.end(); ++it) {

      if (same(it->end(0),p)) {

        // put the connecting stuff into the gap
        L.splice(gap, L, it, L.end());
        path = it; 
        continue;
      }
    }

    // path didn't close.. oh well, keep going with rest of it
    path = gap; 
    p0 = path->end(0);
  }
}

static void toggle_membership(vector<Complex>& b, Complex const& p)
{
  int i, n = b.size();
  for (i=n-1; i>=0; --i) {
    if (same(p, b[i])) {
      b.erase(b.begin()+i);
      return;
    }
  }
  if (i<0) b.push_back(p);
}

void interval_list::get_boundary(vector<Complex>& b) const
{
  list<eqs_interval>::const_iterator it;

  for (it=L.begin(); it!=L.end(); ++it) {
    toggle_membership(b, it->end(0));
    toggle_membership(b, it->end(1));    
  }
}

// return values:
//  0 nothing changed
// -1 loop broken
//  1 whole loops discarded
//  2 new intervals created

int interval_list::clip(bq_curve const& R)
{
  if (!chop(R)) return 0;

  if (!L.size()) return 1; // nothing left. 

  vector<Complex> intersections;
  get_boundary(intersections);

  if (TALK) {
    int oldprec = cout.precision(16);
    cout << "clipping curve " << R << endl; 
    cout.precision(oldprec);

    cout << "chopped\n" << LPSeq(L);
    cout << "intersections\n" << PSeq(intersections) << endl;  
  }

  curve_study C(R, 1.);
  C.add_endpoints(intersections);

  bool loop_broken=false; 
  bool intervals_created=false; 
  eqs_interval I;
  I.set_mark(R.mark()); 
  I.C = R;
  Complex p,t;
  list<eqs_interval>::iterator it, i2;
  for (it=L.begin(); it!=L.end(); ++it) {
    if (it->closed()) continue;
    if (end_connects(it)) continue;

    p = it->end(1);
    t = C.tangent_direction(p);

    if (TALK) {
      cout << "point " << p << " tangent " << t << endl;
    }

    // insert new interval after it and point to it
    ++it; L.insert(it, I); --it;

    if (!it->grow_interval(C,p,t)) {
      loop_broken = true;
      cout << "grow interval failed\n";
      it = L.erase(it); 
      --it; continue;
    }
    if (!end_connects(it)) {
      if (TALK) {
	Complex p=it->end(1), p2; 
	cout << "broken at " << p << endl; 
	i2 = it; ++i2; p2 = i2->end(0); 
	cout << "distance " << complex_modulus(p2-p) 
	     << " same = " << same(p,p2) << endl; 
      }
      loop_broken=true;
    }
    intervals_created = true; 
  }

  tidy();

  // cout << endl; 
  // print();

  if (loop_broken) return -1; 
  return intervals_created ? 2:1;
}

bool interval_list::sort_intervals(double eps)
{
  bool all_loops = true; 

  Complex p;
  list<eqs_interval>::iterator loop, k, j=L.begin(), best_k;
  Mark mark_num; 

  loop = j;
  while (j!=L.end()) {
    p = j->end(1);
    mark_num = j->mark();
    j++;

    // Look for a matching endpoint for the current curve. 
    best_k = L.end();
    for (k=j; k!=L.end(); k++) {
      if (complex_close(p, k->end(0), eps) || complex_close(p, k->end(1), eps)) {
	best_k = k;
	if (mark_num == k->mark()) break; // Use this one rather than any other. 
      }
    }
    k = best_k; 

    // Move it to the appropriate place and set its direction. 
    if (k!=L.end()) {
      if (!complex_close(p, k->end(0), eps)) k->reverse();
      if (j!=k) {
	L.splice(j, L, k);
	j--; 
      }
      continue; 
    }

    if (!complex_close(p, loop->end(0), eps)) {
      all_loops = false; 

      if (j==L.end()) break; 
      // Not a loop. Shift any remaining segments that attach to the 
      // start of the present segment (still called a loop). 

      while (true) {
	p = loop->end(0);
	for (k=j, k++; k!=L.end(); k++)
	  if (complex_close(p, k->end(1), eps) || complex_close(p, k->end(0), eps)) break;
	if (k==L.end()) break; 
	if (!complex_close(p, k->end(1), eps)) k->reverse();
	L.splice(loop, L, k);
	loop--;
      }
    } 
    loop = j;

  }
  return all_loops; 
}

void interval_list::subdivide(Complex const& p, double eps)
{
  eqs_interval a, b; 
  list<eqs_interval>::iterator i; 
  for (i=L.begin(); i!=L.end(); i++)
    if (i->subdivide(p, a, b, eps)) break; 
  if (i==L.end()) return; 

  L.erase(i++); 
  L.insert(i, a);
  L.insert(i, b); 
}

// Modifies lists of intervals A and B. If they each correspond to the 
// same subset of the plane they should both be empty when this function returns. 

void interval_difference(interval_list& A, interval_list& B, double eps, int report)
{
  list<eqs_interval>::iterator ii; 
  int j; 
  for (ii = A.L.begin(); ii != A.L.end(); ii++) {
    for (j=0; j<2; j++)
      B.subdivide(ii->end(j), eps); 
  }
  for (ii = B.L.begin(); ii != B.L.end(); ii++) {
    for (j=0; j<2; j++)
      A.subdivide(ii->end(j), eps); 
  }

  if (report) {
    cout << endl; 
    for (ii=A.L.begin(); ii!=A.L.end(); ii++)
      cout << 'A' << *ii << endl; 
    for (ii=B.L.begin(); ii!=B.L.end(); ii++)
      cout << 'B' << *ii << endl; 
  }

  list<eqs_interval>::iterator ib; 
  ii = A.L.begin();
  while (ii != A.L.end()) {
    for (ib = B.L.begin(); ib != B.L.end(); ib++)
      if (same(*ii, *ib, eps)) break; 
    if (ib == B.L.end()) {
      // No matching interval for *ii; 
      ii++; 
      continue;
    }
    A.L.erase(ii++);
    B.L.erase(ib); 
  }
}

void interval_list::x_reflect()
{
  list<eqs_interval>::iterator i; 
  for (i=L.begin(); i!=L.end(); i++) i->x_reflect(); 
}

bool interval_list::has_mark(Mark m) const
{
  list<eqs_interval>::const_iterator it;
  for (it = L.begin(); it!=L.end(); it++)
    if (it->mark() == m) return true; 
  return false; 
}

void interval_list::print_markings() const 
{
  list<eqs_interval>::const_iterator it; 
  cout << "(";
  for (it = L.begin(); it!=L.end(); it++) {
    if (it != L.begin()) cout << ' '; 
    cout << it->mark(); 
  }
  cout << ")\n"; 
}

void interval_list::compute_bounds(bbox& b, eqsfunc const& surface, double& rad) const
{
  rad = 0.;
  b = bbox(Zero,Zero); 

  if (!L.size()) return;

  list<Complex> segment, curve; 
  list<eqs_interval>::const_iterator k;
  double r; 

  for (k=L.begin(); k!=L.end(); k++) {

    segment = adaptive_curve(linkpic_func(*k, surface), 
			     k->x(0), k->x(1), .2, 0.05);
    curve.splice(curve.end(), segment); 

    r = surface(k->end(0).real, k->end(0).imag).distance();
    if (r > rad) rad = r; 
  }

  b = bbox(curve); 
  b.expand(.02); 
}



// For debugging since function name overloading confuses gdb. 
double eval_poly(R_poly const& p, double x)
{
  return evaluate(p, x); 
}

