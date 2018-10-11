#ifndef _tilt_polytope_
#define _tilt_polytope_
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


#include "VS.hh"
#include "polytope_T.hh"
#include "snappea/SnapPea.h"
#include "picfile.hh"

typedef face<VS> cface;
typedef polytope<VS>::face_cp face_cp;
typedef nface<VS> n_face;

struct tilt_polytope; 

struct tp_face {
  const tilt_polytope* T; 
  n_vector tied_cusp_sizes; // set of (ntc) cusp sizes for tied cusps
  n_vector cusp_size; // set of (nc) cusp sizes yeilding triangulation for this face.

  face_cp F; // provides access to tilt sum hyperplane and live/dead.
  list<n_vector> tilt; // vectors of dim = ntc. 

  tp_face(const tilt_polytope* tp, n_vector const& tcs);

  void insert_tilt(n_vector const& new_tilt);
  void initialize(Triangulation* manifold, n_vector& tilt_sum);
  bool locate_outside_vertex(n_vector& V) const; 
  void save_picture(picfile& pic) const;
};

struct tilt_polytope {
  Triangulation* manifold;
  vector<int> ties; // size = nc. 
  n_vector tie_ratios; // dim = nc. 
  int ntc; // num_tie_classes
  polytope<VS> P; // dimension ntc+1. 
  list<tp_face> FL;
  list<cface> dead_F;

  tilt_polytope(Triangulation* m);
  tilt_polytope(Triangulation* m, vector<int> const&, n_vector const&);

  bool compute(int report, int limit=1000);
  void print() const;
  void save_picture(picfile& pic) const; 
  void get_cusp_sizes(n_vector& cusp_size, n_face const& f) const;
  void expand_with_ties(n_vector& EV, n_vector const& CV) const;
  void contract_with_ties(n_vector& CV, const double EV[]) const;

private:
  void add_face(n_vector const& cusp_size); 
  void set_default_ties();
};

bool on_boundary(n_vector const& v);

#endif
