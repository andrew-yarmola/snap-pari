#ifndef _point_
#define _point_
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

#include "snappea/O31.h"
#include "R_matrix.hh"

using std::cout;

class point : public R_vector<3> {
public:
  point() {}

  double& x() { return (*this)[0]; }
  double& y() { return (*this)[1]; }
  double& z() { return (*this)[2]; }  

  point(double _x, double _y, double _z)
    { x()=_x; y()=_y; z()=_z; }
  point(const O31_vector& v) 
    { x()=v[2]/v[0]; y()=v[3]/v[0]; z()=v[1]/v[0]; }
  point(R_vector<3> const& v) : R_vector<3>(v) {}

  double x() const { return (*this)[0]; }
  double y() const { return (*this)[1]; }
  double z() const { return (*this)[2]; }

  operator O31_vector() const
    { return O31_vector(1.0,z(),x(),y()); }

};

inline point operator * (double r, const point& v)
{ return r*((R_vector<3>)v); }

point operator * (const O31_matrix& m, const point& p); 

double hyp_distance(const point& a, const point& b); 
double euc_distance(const point& a, const point& b);
double alt_distance(const point& a, const point& b);

#define SMALL 1e-6

inline bool num_zero(double d)
{
  return (d>0.0) ? d < SMALL : -d < SMALL;
}

inline void mma_print(const point& p, ostream& out = cout)
{ 
  out << '{' << p[0] << ',' << p[1] << ',' << p[2] << '}'; 
}

inline void gv_print(const point& p, ostream& out = cout)
{ 
  out << p[0] << ' ' << p[1] << ' ' << p[2]; 
}

inline double sort_value(const point& p) 
{ 
  return p[0]; 
}

inline point crossing_point(point a, double a_val, point b, double b_val)
{
  return (1.0/(a_val - b_val)) * 
    (a_val * (R_vector<3>)b - b_val * (R_vector<3>)a);
}

inline double dot(const point& a, const point& b)
{
  return a*b;
}

#if 0
// Now these all return R_vector<3> instead of point
// but that should be OK since a point can be implicitly 
// constructed from an R_vector<3>. 

inline point operator * (double d, const point& p)
{
  return point(d * p[0], d * p[1], d * p[2]); 
}

inline point operator - (const point& a, const point& b)
{
  return point(a[0]-b[0], a[1]-b[1], a[2]-b[2]); 
}

inline point operator + (const point& a, const point& b)
{
  return point(a[0]+b[0], a[1]+b[1], a[2]+b[2]); 
}

inline int operator == (const point& a, const point& b)
{
  return num_zero(a[0]-b[0]) && num_zero(a[1]-b[1]) && num_zero(a[2]-b[2]); 
}

inline int operator != (const point& a, const point& b)
{ 
  return !(a==b);
}

inline int operator < (const point& a, const point& b)
{
  return a[0] < b[0];
}
#endif

#endif
