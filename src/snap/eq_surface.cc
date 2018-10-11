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
#include "eq_surface.hh"
#include <float.h>
#include "complex_lll.hh"  // for get_dep. 
#include "eqs_interval.hh" // for compute_bounds();

static double atan_2(double y, double x, double approx_arg)
{
    if (x == 0.0  &&  y == 0.0)
    {
	printf("atan2(0,0) encountered");
	return approx_arg;
    }

    double arg = atan2(y, x);
    while (arg - approx_arg > PI)
	arg -= TWO_PI;
    while (approx_arg - arg > PI)
	arg += TWO_PI;

    return arg;
}

// class uc_point

uc_point crossing_point(uc_point a, double a_val, uc_point b, double b_val)
{
  return uc_point((1.0/(a_val - b_val)) * (a_val * b - b_val * a), 
		    ((a_val * b.arg - b_val * a.arg)/(a_val - b_val)));
}

double uc_point::distance() const
{
  double d = 1.0 - norm_squared();
  if (d <= 0.0) return DBL_MAX;
  return acosh(sqrt((1.0 - z()*z())/d));
}

void uc_point::uc_print(ostream& out) const
{
  out << '[' << Complex(*this) << ", ";
  double d = distance();
  if (d > 1e4) out << "Infinity";
  else out << d; 
  out << ']'; 
}


void uc_point::uc_print(ostream& out, vector<Complex> const& H) const
{
  Complex a = Complex(*this); 
  double C[2]; 
  int i;
  get_dep(H, a, C);
  for (i=0; i<2; i++) C[i] -= floor(C[i]); 
  out << '[' << Complex(C[0],C[1]) << ']';
}

uc_point& uc_point::transform_by(const O31_matrix& m, double appr_arg)
{
  *this = uc_point(m * *this, appr_arg); 
  set_correct_arg(); 
  return *this;
}

uc_point& uc_point::set_correct_arg()
{
  double ar = atan_2(y(), x(), arg);
  arg = ar;
  return *this;
}

uc_point uc_transformation::operator () (const uc_point& p) const
{ 
  if (z_axis) // Screw motion along z axis. 
    return uc_point(m * p, ar + p.arg); 

  uc_point pc(p); 
  pc.transform_by(m,ar); 
  return pc; 
}

// class eqsfunc

// Construct the equidistant surface between the geodesic 0,infinity and
// its image a,b under the given Moebius transformation. 

eqsfunc::eqsfunc(const MoebiusTransformation& m)
{
  Complex a, b; 
  a = m.matrix[0][1]/m.matrix[1][1];
  b = m.matrix[0][0]/m.matrix[1][0];

  complex_pair_to_ol_od(a, b, gz, gx); 
  set_implicit_desc(); 
  compute_bounds();
}

void eqsfunc::set_implicit_desc()
{
  O31_matrix tr = O31_x_trans(gz) * O31_y_trans(gx);
  O31_vector pt, dir;

  pt = tr * O31_vector(1.,0.,0.,0.);
  dir = tr * O31_vector(0.,1.,0.,0.); 

  a = pt[2] * pt[2] - dir[2] * dir[2];
  b = pt[3] * pt[3] - dir[3] * dir[3];

  c = 1.0 + pt[1] * pt[1] - dir[1] * dir[1];

  d = 2.0 * pt[2] * pt[3] - 2.0 * dir[2] * dir[3];
  e = 2.0 * pt[2] * pt[1] - 2.0 * dir[2] * dir[1];
  f = 2.0 * pt[3] * pt[1] - 2.0 * dir[3] * dir[1];

  g = -2.0 * pt[2] * pt[0] + 2.0 * dir[2] * dir[0];
  h = -2.0 * pt[3] * pt[0] + 2.0 * dir[3] * dir[0];
  i = -2.0 * pt[1] * pt[0] + 2.0 * dir[1] * dir[0];

  j = -1.0 + pt[0] * pt[0] - dir[0] * dir[0];

  K = sin(gx.imag)/sinh(gx.real);
  T = O31_x_trans(gz) * O31_y_trans(0.5 * gx);

  // cout << '[' << gx << ',' << gz << ']' << endl; 
}

void eqsfunc::compute_bounds()
{
  R_matrix<3> Sph(0.0);
  Sph[0][0] = 1.0; 
  Sph[2][0] =-1.0; 
  Sph[0][2] =-1.0; 
  Sph[2][2] =-K_value()*K_value(); 
  bq_curve S(Sph, Mark(), position().imag); 

  double rad;
  interval_list boundary;
  boundary.make_intervals(S);
  boundary.compute_bounds(box,*this,rad);
}

double eqsfunc::operator() (const uc_point& p) const
{
  if (fabs(p.arg - (gz.imag)) > PI) return 1.0; 

  return a * p.x() * p.x() + b * p.y() * p.y() + c * p.z() * p.z() + d * p.x() * p.y() + 
    e * p.x() * p.z() + f * p.y() * p.z() + g * p.x() + h * p.y() + i * p.z() + j; 
}

uc_point quadratic::operator() (const uc_point& p1, const uc_point& p2) const
{
  double v1 = (*this)(p1); 
  double v2 = (*this)(p2);
  double vm = (*this)(crossing_point(p1,-1,p2,1)); 

  double A =  2.0*v1 + 2.0*v2 - 4.0*vm; 
  double B = -3.0*v1 -     v2 + 4.0*vm; 
  double C =      v1; 

#if 0
  if (euc_distance(p1,p2) > 0.5) {
    cout << "** points too far apart in crossing point func. **\n"; 
  }
#endif

  if (fabs(A) < 1e-6) {
    // cout << "** quadratic too flat, using linear appr. **\n"; 
    return crossing_point(p1,v1,p2,v2); 
  }
  double sgn = (v1 < v2) ? 1.0 : -1.0;

  double t = (-B + sgn * sqrt(B*B - 4.0*A*C))/(2.0*A); 
  if (t < -.001 || t > 1.001) {
    cout << "** t out of range in crossing point func. trying other sign. **\n"; 
    t = (-B - sgn * sqrt(B*B - 4.0*A*C))/(2.0*A); 
  }

  uc_point res(t*p2 + (1.0-t)*p1, t*p2.arg + (1.0-t)*p1.arg);
  if (fabs((*this)(res)) > 1e-5)
    cout << "** crossing point func. produced inaccurate result. **\n"; 
  return res; 
}

uc_point eqsfunc::operator() (const uc_point& p1, const uc_point& p2) const
{
  return quadratic::operator() (p1, p2); 
}

uc_point eqsfunc::operator() (const double& x, const double& y) const
{ 
  O31_vector pt = T * O31_vector(1.0, y, K*x*y, x);

  double arg = atan_2(pt[3], pt[2], gz.imag); 

  return uc_point(pt, arg);
}

R_matrix<3> eqsfunc::restriction_matrix(const eqsfunc& f) const
{
  // Set up matrix of Kxz - yw. 
  GL4R_matrix M0;
  int i, j;
  for (i=0; i<4; i++) 
    for (j=0; j<4; j++) 
      M0(i,j)=0.0;
  M0(0,2)=-1.0;
  M0(1,3)=f.K;
  M0(3,1)=f.K;
  M0(2,0)=-1.0;

  // Compute matrix of restriction 
  O31_matrix Tdiff = inverse(f.T) * T; 
  GL4R_matrix Mn = transpose(Tdiff) * M0 * Tdiff; 
  R_matrix<3> R; 
  R[0][0] = Mn(0,0);
  R[0][1] = 2.0 * Mn(0,3);
  R[0][2] = Mn(3,3);
  R[1][0] = 2.0 * Mn(0,1);
  R[1][1] = 2.0 * (Mn(0,2) * K + Mn(1,3));
  R[1][2] = 2.0 * Mn(2,3) * K;
  R[2][0] = Mn(1,1);
  R[2][1] = 2.0 * Mn(1,2) * K;
  R[2][2] = Mn(2,2) * K * K;

  return R;
}

// Because an eqsfunc represents a quadratic function, its derivative
// is a linear function of the point and can thus conveniently be
// represented by a matrix and vector combination A, V. The derivative
// of eqsfunc at a point P is then A*P + V. 

void eqsfunc::get_derivative(R_matrix<3>& A, R_vector<3>& V) const
{
  A[0][0] = 2.0*a; 
  A[1][1] = 2.0*b; 
  A[2][2] = 2.0*c; 
  A[0][1] = d; 
  A[1][0] = d; 
  A[0][2] = e; 
  A[2][0] = e; 
  A[1][2] = f; 
  A[2][1] = f; 
  V[0] = g;
  V[1] = h;
  V[2] = i; 
}

void spherefunc::get_derivative(R_matrix<3>& A, R_vector<3>& V) const
{
  A[0][0] = -2.0; 
  A[1][1] = -2.0; 
  A[2][2] = -2.0; 
  A[0][1] = 0.0; 
  A[1][0] = 0.0; 
  A[0][2] = 0.0; 
  A[2][0] = 0.0; 
  A[1][2] = 0.0; 
  A[2][1] = 0.0; 
  V[0] = 0.0;
  V[1] = 0.0;
  V[2] = 0.0; 
}

// p is an approximation to a point on the common intersection of the
// surfaces s0, s1 and s2. The function uses Newton's method to find 
// a more accurate value for the common intersection point (closest to p). 

bool polish_point(quadratic const& s0, quadratic const& s1, quadratic const& s2, uc_point& p)
{
  R_vector<3> err;
  R_matrix<3> D; 
  int i, max_steps = 20; 

  R_matrix<3> DsA[3];
  R_vector<3> DsV[3];

  s0.get_derivative(DsA[0],DsV[0]); 
  s1.get_derivative(DsA[1],DsV[1]); 
  s2.get_derivative(DsA[2],DsV[2]); 

  for (i=0; i<3; i++)
    D.set_row(i, DsA[i] * p + DsV[i]);
  // cout << D << endl;
  if (fabs(D.determinant()) < .001) 
    return false; 

  while (true) {

    err[0] = s0(p);
    err[1] = s1(p);
    err[2] = s2(p);

    for (i=0; i<3; i++)
      D.set_row(i, DsA[i] * p + DsV[i]);

    // err:= D^-1 * err. 
    D.gauss(err); 

    p -= err; 
    if (err.is_zero() || --max_steps == 0) break;
  }
  return (max_steps != 0); 
}


Complex inverse_eqsfunc::operator() (const uc_point& p) const
{
  point npt = Tinv * p; 
  return Complex(npt.y(), npt.z()); 
}

#if 0
polylist<uc_point> eqsfunc::get_polylist(double step) const
{
  polylist<uc_point> m = polylist<uc_point>(*this, -1., 1., step, -1., 1., step);
  m.clip(spherefunc()); 

  return m;
}
#endif

void complex_pair_to_ol_od(const Complex& a, const Complex& b, Complex& ol, Complex& od)
{
  Complex g_mean = complex_sqrt(a * b); 

  if (fabs((a/g_mean).arg()) > PI/2) 
    g_mean = -g_mean; 
  ol = complex_log(g_mean, 0.0); 

  Complex nz_a = a/g_mean; 

  od = complex_log((One + nz_a)/(One - nz_a), 0.0);
}

eqspt_gen::eqspt_gen(const quadratic& orig, const quadratic& cl, const vector<eqsfunc>& su)
: surfaces(&su), original(&orig), clipping(&cl), sphere(1.0) 
{}

eqspt_gen::eqspt_gen()
: surfaces(0), original(0), clipping(0), sphere(1.0) 
{}

uc_point eqspt_gen::operator() (const uc_point& a, const uc_point& b, const int& mark) const
{
  uc_point c0 = (*clipping)(a,b); 
  if (mark==0) return c0; 

  // these didn't have side effects did they?
  // (*clipping)(a); 
  // (*clipping)(b); 

  uc_point c; 

  c = c0; 
  const quadratic *s; 
  if (mark>0) s = &((*surfaces)[mark-1]); 
  else s = &sphere; 
  if (!polish_point(*original, *clipping, *s, c)) {
    // cout << "** polishing failed in eq_surface.cc. **\n";  
    return c0; 
  }

#if 0
  cout << '[' << (*original)(a) << ',' << (*original)(b) << ',';
  cout << (*s)(a) << ',' << (*s)(b) << ';' << mark << ']' << endl;

  if (alt_distance(c, c0) > 0.05) {
    cout << "** inaccurate result in eq_surface.cc. **\n"; 
    return c0;
  }
#endif

  return c; 

}
