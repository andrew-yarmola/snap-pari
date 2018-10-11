#ifndef _O31_line_
#define _O31_line_
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


// A based geodesic class using O31 vectors to specify a point and direction. 

#include "snappea/O31.h"
#include "lines.hh"

class O31_line {
  O31_matrix trans; // line is image of x-axis.

  void set_line(O31_vector const& pt, O31_vector const& dir);
public:

  O31_line() {}
  O31_line(const line& l);

  O31_line(const O31_matrix& m) : trans(m) {}
  O31_line(O31_vector p, O31_vector d);
  O31_line(const O31_vector& a, const O31_vector& b, int); // a,b on light cone.

  operator line() const; 
  operator O31_matrix() const { return trans; } 

  O31_matrix& T() { return trans; }
  O31_matrix const& T() const { return trans; }

  O31_vector point() const
  { return trans.col(0); }
  O31_vector point(double s, double c) const;
  O31_vector point(double t) const
  { return point(sinh(t),cosh(t)); }

  O31_vector direction() const
  { return trans.col(1); }
  O31_vector direction(double t) const
  { return direction(sinh(t),cosh(t)); }
  O31_vector direction(double s, double c) const;

  O31_vector end(int i) const; 


  // void normalize(); 
  void transform_by(const O31_matrix& mx) { trans.left_mul(mx); }
  void adjust_basepoint(O31_vector p); 

  int crossing_parameter(const O31_vector& normal, double& s, double& c) const;
  bool is_in_plane(const O31_vector& normal) const; 

  O31_vector projection(const O31_vector& p) const;

  double point_parameter(const O31_vector& p) const
  { return asinh(p * trans.col(1)); }

  friend O31_line operator * (const O31_matrix& mx, const O31_line& l)
  { return mx * l.trans; }

  friend double distance_to_org(const O31_line& l); // distance to origin. 
  friend double distance_to_org(const O31_line& l, const O31_vector& p);

  friend ostream& operator << (ostream& out, const O31_line& l);

  static double epsilon;
};

int compare(O31_line const& a, O31_line const& b, double eps=O31_line::epsilon);

#endif 
