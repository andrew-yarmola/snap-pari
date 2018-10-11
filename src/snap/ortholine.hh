#ifndef _ortholengths_
#define _ortholengths_
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
#include "snappea/length_spectrum.h"
#include "closed_geodesics.hh"

#define black rb_tree_black // Avoid name collision with black from color.hh. 
#include <set>
#undef black

using std::set;

const double ort_end_eps = 1e-5;

// A CGRay is a ray orthogonal to a closed geodesic. 

struct CGRay {
  GSpec g;      // The closed geodesic.
  Complex pos;  // Position, mod (g.g->length, 2*Pi*I).

  CGRay() {}
  CGRay(GSpec const& gs, Complex const& p=Zero) : g(gs), pos(p) {}
  CGRay(Complex const& p, int gn) : g(gn), pos(p) {}
};

int cmp(CGRay const& a, CGRay const& b, double eps=ort_end_eps);
ostream& operator << (ostream& out, CGRay const& r);

/* An Ortholine keeps track of the position of a geodesic with respect
   to some given base geodesic. The geodesic is given by giving the
   common orthogonal and the complex distance of the geodesic from the
   base geodesic along the common orthogonal. The common orthogonal is
   given by giving the log of its positive endpoint after normalizing
   so that the base geodesic goes to [0,infinity].
*/


class Ortholine {
  Complex _distance; // The orthodistance.
  CGRay end[2];      // Ends of the ortholine segment.
  int low_end;       // Which end is lower wrt. cmp.

public:
  FGWord word;       // Conjugacy carrying geodesic at one end to the other. 

  Ortholine() : low_end(0) {}
  Ortholine(Complex const& d, CGRay const& a, CGRay const& b);
  Ortholine(MoebiusTransformation const& m, GSpec const& a, GSpec const& b); 

  // Ortholine(MoebiusTransformation const& m, int na, int nb, FGWord const& w); 

  const Complex& position(int i=1) const { return end[i].pos; }
  const Complex& distance() const { return _distance; }
  int geodesic_num(int i) const { return end[i].g.index; }

  // Complex& position(int i=1) { return end[i].pos; }
  // Complex& distance() { return _distance; }

  line get_line(int i) const; 

  void sort_ends();

  void set_to(Complex const& d, Complex const& a, Complex const& b);

  void set_from_Moebius(MoebiusTransformation const& M); 
  operator MoebiusTransformation() const; 

  friend bool operator < (const Ortholine& a, const Ortholine& b); 
  friend ostream& operator << (ostream& out, const Ortholine& o); 
  
private:
  void set_low_end() { low_end = (cmp(end[0],end[1])==1) ? 1:0; }
};

typedef set<Ortholine> Ortholine_set; 

void get_ortholines(Tile* tiling, list<interval> const& ivls, double tile_rad_lo, double tile_rad_hi, Ortholine_set& ort_set, double cutoff=1000., bool check_range=true);

void new_ortholines(const WEPolyhedron* poly, list<interval> const& li, std::vector<GeodesicWord> const& GL, double radius, Ortholine_set& OS);

inline bool operator == (const Ortholine& a, const Ortholine& b)
{
  return !(a < b || b < a);
}

bool lxless(Complex const& a, Complex const& b, double eps);

#endif
