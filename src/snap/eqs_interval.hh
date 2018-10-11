#ifndef _eqs_interval_
#define _eqs_interval_
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


#include <list>
#include <iostream>
#include <algorithm>
#include "R_poly.hh"
#include "R_matrix.hh"
#include "eq_surface.hh"
#include "mark.hh"
#include "bbox.hh"

using std::list;
using std::ostream;

double evaluate(const R_matrix<3>& R, Complex const& z);
Complex evaluateD(const R_matrix<3>& R, Complex const& z);

// Each bq_curve comes from an equidistant surface. The mark num
// keeps track of which equidistant surface, the arg keeps track of 
// which branch of the cover of H^3-axis we are in. 

class bq_curve
{
  R_matrix<3> M;       // R[r][c] is coeff of y^r x^c. Eval as Y * R * X.
  vector<R_poly> F[2]; // F[0][r] is coeff of y^r as poly in x. 
                       // F[1][r] is coeff of x^r as poly in y.
  Mark mark_num;
  double _arg; 
public:
  bq_curve() {}
  bq_curve(const R_matrix<3>& m, Mark mk=Mark(), double a=0.) : 
    M(m), mark_num(mk), _arg(a) {}

  vector<R_poly> const& f(int xy) const;
  // void set_f(int xy, vector<R_poly> const& f) { F[xy] = f; }

  double get_yd(int xy, double x, double& disc, int* num_pts=0) const;
  double Y(int xy, int sign, double x) const;
  Complex at(int xy, int sign, double x) const;
  bool get_branch(int xy, Complex const& p, int& sign, int dir=0) const;
  // bool resultant_vanishes(int xy, const R_matrix<3>& M) const;

  void scale(double s); 
  void x_reflect();

  Mark mark() const { return mark_num; }
  void set_mark(Mark mk) { mark_num = mk; }
  R_matrix<3> const& Mx() const { return M; }
  double arg() const { return _arg; }

  bool factorize(vector<bq_curve>& factors) const;

  friend ostream& operator << (ostream& out, bq_curve const& m); 
  bool read(FILE* fp=stdin, bool get_marking=true);
};

vector<Complex> intersection_points(bq_curve const& A, bq_curve const& B, int xy, bool irred=false);

struct branch {
  int xy;          // 0=x, 1=y
  int sign;        // -1 for lower, 1 for upper

  branch() {}
  branch(int _xy, int s) : xy(_xy), sign(s) {}

  friend ostream& operator << (ostream& out, const branch& b);
  friend bool operator == (branch const& a, branch const& b)
  { return a.xy==b.xy && a.sign==b.sign; }
};

class curve_study;

class eqs_interval {
  bq_curve C;

  vector<Complex> P;  // points along the curve
  vector<branch> BR;  // branch specifier between points

public:

  static double epsilon;
  static int talk;

  eqs_interval() {}

  Mark mark() const { return C.mark(); }
  void set_mark(Mark m) { C.set_mark(m); }

  Complex operator () (double param) const; 
  double x(int e) const; 
  Complex end(int e) const { return (e==0) ? P.front():P.back(); }
  bool closed() { return P.size()==BR.size(); }

  list<Complex> get_polyline() const; 
  void x_reflect();

  bool read(FILE* fp);
  void write(ostream& out); 

  friend ostream& operator << (ostream& out, const eqs_interval& i);
  friend bool operator == (eqs_interval const& a, eqs_interval const& b); 
  friend bool same(eqs_interval const& a, eqs_interval const& b, double eps);

  bool check_me() const; 

private:

  eqs_interval(eqs_interval const& i, int a, int b); // subinterval

  Complex at(int b, double x) const;
  double x(int b, int e) const 
  { int j=(b+e)%P.size(); return BR[b].xy ? P[j].imag : P[j].real; }

  bool in_range(int b, double t) const { return (t-x(b,0))*(x(b,1)-t) >= 0.; }
  void reverse();

  int branch_num(double t, double& xb) const;
  int branch_num(Complex const& p) const;

  bool subdivide(Complex const& p, eqs_interval& a, eqs_interval& b, double eps) const;
  int subdivide(vector<Complex> const& S, vector<int>& indices, vector<Complex>* used=0);
  void close(); 
  void insert_point(Complex const& p, int bn, vector<int>& indices);
  double xp(int e) const; 
  Complex a_point(int a, int b) const;

  friend class interval_list; 

  bool grow_interval(curve_study& CS, Complex& p, Complex& t);
  bool make_interval(curve_study& CS);
};

class interval_list {
public:
  list<eqs_interval> L;

  void make_intervals(bq_curve const& R);
  bool chop(bq_curve const& R); 
  int clip(bq_curve const& R);

  bool sort_intervals(double eps); 

  void print(ostream& out) const;
  void print_sorted(ostream& out, double eps) const;
  void print_markings() const;

  void x_reflect();
  bool has_mark(Mark m) const;
  void compute_bounds(bbox& b, eqsfunc const& surface, double& rad) const;

  friend void interval_difference(interval_list& A, interval_list& B, double eps, int report=0);

private:

  void push_back(eqs_interval const& i) { L.push_back(i); }
  void subdivide(Complex const& p, double eps);
  bool end_connects(list<eqs_interval>::const_iterator it0) const;
  void tidy();
  void get_boundary(vector<Complex>& b) const;
};

typedef list<eqs_interval>::iterator eqs_iter;
typedef list<eqs_interval>::const_iterator eqs_const_iter;

class close_test {
  double epsilon;
public:
  close_test(double e) : epsilon(e) {}
  bool operator () (double a, double b) const 
    { return fabs(b-a) < epsilon; }
};

class linkpic_func {
  const eqs_interval& ivl;
  const eqsfunc& surface; 
public:
  linkpic_func(const eqs_interval& i, const eqsfunc& s) : ivl(i), surface(s) {}

  Complex operator() (double x) const
    { Complex z(ivl(x)); return surface(z.real, z.imag); }
}; 

#endif 
