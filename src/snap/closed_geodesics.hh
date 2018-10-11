#ifndef _closed_geodesics_
#define _closed_geodesics_
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


#include "snappea/SnapPea.h"
#include "O31_line.hh"
#include <vector>
#include <set>

using std::list;

class DirichletDomain;

// Often it is useful to pass around information about a 
// geodesic of the manifold, but when printing it out, all 
// we want to see is a number. Since this number is not 
// really intrinsic to the geodesic we keep it separate. 

struct GSpec {
  int index;
  const GeodesicWord* g;

  GSpec(int gn=-1, const GeodesicWord* gwp=0) 
    : index(gn), g(gwp) {}
};

// An interval represents a directed, parametrized geodesic segment in 
// hyperbolic space. 
// An interval represents the image under trans of the  
// segment from lo to hi of the x-axis (w,x,0,0). 

// Seen from the point of view of a Dirichlet domain for a 3-manifold, a
// closed geodesic looks like a bunch of geodesic segments. The function
// all_crossing_lifts computes such a set of segments from a closed 
// geodesic given by point, direction and length. 

class interval {
  O31_line L;
  double lo, hi; 
  FGWord word; 
  GSpec g;
public:

  static double ideal_end_eps;

  interval() {}

  interval(O31_line const& ln, double l, double h, FGWord const& wd, GSpec const& gs) 
    : L(ln), lo(l), hi(h), word(wd), g(gs) {}

  void transform(O31_matrix const& m, FGWord const& w)
  { L.transform_by(m); word.left_multiply_by(w); }

  // Default copy and assignment OK. 

  O31_vector point(double t) const
  { return L.point(t); }
  O31_vector end(int i) const
  { return L.point((i==0)?lo:hi); }
  O31_vector direction(int i) const
  { return L.direction((i==0)?lo:hi); }

  O31_vector midpoint() const
  { return L.point((lo + hi)/2.0); }

  FGWord const& the_word() const { return word; }
  O31_line const& the_line() const { return L; }
  GSpec const& geodesic() const { return g; } 
  int gnum() const { return g.index; } 

  bool in_range(double t) const 
  { return t > lo - O31_line::epsilon && t < hi + O31_line::epsilon; }
  void set_range(double l, double h) 
  { lo = l; hi = h; }

  friend class DirichletDomain;

  friend bool operator == (interval const&, interval const&); 
  friend ostream& operator << (ostream& out, interval const& ivl);
};

bool all_crossing_lifts(const WEPolyhedron *poly, O31_line const& L, GSpec const& g, list<interval>& l, bool full_set);
bool all_crossing_lifts(const WEPolyhedron *poly, interval const& I, list<interval>& l, bool full_set);

bool find_geodesic(const WEPolyhedron* domain,
		   std::vector<GeodesicWord> const& geodesics, 
		   O31_line const& g, 
		   int start_g_num, int& gn, int& ori, FGWord& conj);

double outradius(list<interval> const& l);

void dirichlet_normalize(const WEPolyhedron* polyhedron, O31_vector& point, O31_matrix& tr, FGWord& nw, int report);
bool is_symmetry_of_geodesic_lifts(const WEPolyhedron* poly, list<interval> const& li, O31_matrix T, int report);

// line_wd stands for line with distance (from the origin). 

class line_wd {
  O31_line l;
  double d; // distance from origin
public:
  int gnum; 

  line_wd() {}
  line_wd(O31_line const& l0, int gn=0) : l(l0), d(distance_to_org(l)), gnum(gn) {}

  O31_line const& L() const { return l; }
  double distance() const { return d; }

  friend bool operator < (line_wd const& a, line_wd const& b);
  friend line_wd operator * (O31_matrix const& M, line_wd const& L);
};


int get_lifts(const WEPolyhedron* poly, list<interval> const& il, double radius, std::set<line_wd>& S);

#endif
