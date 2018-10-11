#ifndef _CONNECTED_FACES_H_
#define _CONNECTED_FACES_H_
/*
** Copyright (C) 2004 Oliver A. Goodman <oag@ms.unimelb.edu.au>
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

// For inclusion in tube.cc only, so we'll assume
// types used there are already defined. 

#include "uc_polylist.hh"
#include "triangulation_builder.hh"

class connected_faces {
  const tube* T;

  vector<O31_matrix> pairing_matrices; // One per ortholine. 

  vector<uc_polygon> FL; // Corners of the faces.
  vector<int> face_num; // Relation with previous numbering scheme. 

  triangulation_builder* TB;

  O31_matrix pairing_matrix(int face) const; 
  O31_matrix cf_pairing_matrix(int face) const; 
  face_corner adjacent_corner(point const& p, Mark const& nbr_face, uc_vertex& v) const;
  face_corner paired_corner(const point& p, int face, uc_vertex& v) const;
  int subface_mark_number(Mark const& mark, point const& p, vector<int> const& subface_start, vector<int> const& subface_count) const;
  uc_vertex vertex(face_corner const& fc) const;
  Complex edge_holonomy(face_corner fc) const;
  bool initialize_triangulation_bldr(int report=0);
  bool get_triangulation_data(TriangulationData& TD, int report=0);

public:

  connected_faces(const tube* t) : TB(0), T(t) {}

  O31_matrix trans(Mark const& mark) const
  { return O31_x_trans(T->holonomy(mark)); }
  void set_pairings(vector<Ortholine> const& ol);
  bool compute_connected_faces(vector<tube_face>& faces);
  void orient_connected_faces();

  void print_connected_faces() const; 
  void print_face_corners() const;

  Triangulation* drill(); 
  void print_triangulation_data();
  void recover_holonomies(Triangulation* m, Complex H[2]) const;
};

#endif
