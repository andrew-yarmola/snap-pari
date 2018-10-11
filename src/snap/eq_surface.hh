#ifndef _tube_
#define _tube_
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


#include "point.hh"
#include "vfunction.hh"
#include <vector>
#include "bbox.hh"

using std::vector;
using std::ostream;

// uc is universal cover of H^3, (Klein/Beltrami/projective model) minus the z-axis. 
// uc_point adds an argument to point to say which branch of 
// uc the point lies on. (x,y) = k (cos(arg),sin(arg)) for suitable k>0. 

class uc_point : public point {
public:
  double arg;

  uc_point() {}
  uc_point(const point& pt, double ar) : point(pt), arg(ar) {}
  uc_point(double a, double b, double c, double ar)
    : point(a,b,c), arg(ar) {}

  uc_point& transform_by(const O31_matrix& m, double appr_arg); 
  uc_point& set_correct_arg();

  double distance() const; 
  operator Complex () const 
    { return Complex(atanh((*this)[2]), arg); }
  void uc_print(ostream& out) const; 
  void uc_print(ostream& out, vector<Complex> const& H) const; 
};

uc_point crossing_point(uc_point a, double a_val, uc_point b, double b_val);

inline int operator == (uc_point const& a, uc_point const& b)
{
  return euc_distance(a,b) < 1e-6 && 
    fabs(a.arg-b.arg) < 1e-6; 
}

inline int operator != (uc_point const& a, uc_point const& b)
{ return !(a==b); }

class uc_transformation 
: public unary_vfunction<uc_point, uc_point> 
{
  O31_matrix m; 
  double ar; 
  bool z_axis; 
public:
  uc_transformation() {}

  uc_transformation(const O31_matrix& mx, double a) 
    : m(mx), ar(a), z_axis(false) {}

  uc_transformation(const Complex& t) 
    : m(O31_x_trans(t)), ar(t.imag), z_axis(true) {}

  // Memberwise copy and assignment OK. 

  virtual uc_point operator() (const uc_point& p) const;
};

class quadratic 
: public unary_vfunction<uc_point, double>,
  public binary_vfunction<uc_point, uc_point, uc_point>
{
public:
  virtual void get_derivative(R_matrix<3>& A, R_vector<3>& V) const = 0;
  virtual double operator() (const uc_point& pt) const = 0;
  // crossing pt func.
  virtual uc_point operator() (const uc_point& p1, const uc_point& p2) const; 
};


class spherefunc 
: public quadratic
{
  double rad_sq; 
public:
  spherefunc(double r = 1.0) { rad_sq = r*r; }

  double operator() (const uc_point& p) const
    { return rad_sq - p.norm_squared(); }
  void get_derivative(R_matrix<3>& A, R_vector<3>& V) const;
}; 

class eqsfunc 
: public quadratic,
  public binary_vfunction<double, double, uc_point>
{
  Complex gx; // coords of geodesic, z-axis translated along x-axis, then z-axis. 
  Complex gz; 

  double a,b,c,d,e,f,g,h,i,j;

  double K; // the "curvature" constant of the surface, K = sin(gx.imag)/sinh(gx.real); 
  O31_matrix T; // gives the transformation from normalized location to location on tube.
  // T = O31_x_trans(gz) * O31_y_trans(0.5 * gx);

  bbox box;

  void set_implicit_desc(); 
  void compute_bounds();
public:

  eqsfunc() {}
  eqsfunc(const MoebiusTransformation& m); 
  eqsfunc(const Complex& d, const Complex& a) : gx(d), gz(a) 
  { set_implicit_desc(); compute_bounds(); } 
  eqsfunc(eqsfunc const& e, Complex const& z) : gx(e.gx), gz(e.gz+z), box(e.box+z)
  { set_implicit_desc(); }

  // memberwise copy and assignment are OK. 

  uc_point operator() (const double& x, const double& y) const; // parametrization of surface
  double operator() (const uc_point& pt) const; // implicit function for surface
  uc_point operator() (const uc_point& p1, const uc_point& p2) const; // crossing pt func.

  R_matrix<3> restriction_matrix(const eqsfunc& f) const;
  void get_derivative(R_matrix<3>& A, R_vector<3>& V) const;
  double K_value() const { return K; }
  Complex const& geodesic_distance() const { return gx; }
  Complex const& position() const { return gz; }
  bbox const& bounds() const { return box; }

  friend eqsfunc operator + (const eqsfunc& F, const Complex& z) 
    { return eqsfunc(F, z); }
  friend class inverse_eqsfunc;
};

// inverse_eqsfunc is a minimal function-object which we construct
// from an eqsfunc. its purpose is to recover parameter values from
// points on the equidistant surface. although this task 
// could, in principle, just as easily be handled by eqsfunc itself, 
// the required "Complex operator() (const uc_point& pt) const;"
// would conflict with the existing "double operator() (const uc_point& pt) 
// const;". 
class inverse_eqsfunc
: public unary_vfunction<uc_point, Complex>
{
  O31_matrix Tinv;
public:
  inverse_eqsfunc() {}
  inverse_eqsfunc(const eqsfunc& ef) : Tinv(inverse(ef.T)) {}

  Complex operator() (const uc_point& pt) const;
  // Returns the parameter value of point pt in this equidistant surface. 
};

class eqspt_gen
: public ternary_vfunction<uc_point, uc_point, int, uc_point>
{
  const vector<eqsfunc>* surfaces;
  const quadratic* original;
  const quadratic* clipping; 
  spherefunc sphere; 
public:

  eqspt_gen(); 
  eqspt_gen(const quadratic& orig, const quadratic& cl, const vector<eqsfunc>& su); 

  uc_point operator() (const uc_point& a, const uc_point& b, const int& mark) const; 
  // Returns the point on original, clipping and surfaces[mark-1] which should
  // be somewhere between a and b, for mark >0. If mark==0, returns point
  // on geodesic a,b and clipping. If mark <0, returns point on original, clipping
  // and sphere. 
};

bool polish_point(quadratic const& s0, quadratic const& s1, quadratic const& s2, uc_point& p);

void complex_pair_to_ol_od(const Complex& a, const Complex& b, Complex& ol, Complex& od);

#endif
