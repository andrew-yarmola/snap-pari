#ifndef _triangulation_builder_
#define _triangulation_builder_
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
#include <vector>
#include <map>

using std::map;
using std::vector;
using std::ostream;

enum fc_type { fc_back, fc_front, fc_top, fc_ideal };

class face_corner {
public:
  face_corner() : face(-1), vertex(0), type(fc_back) {}
  face_corner(int f, int v, fc_type t=fc_back) : face(f), vertex(v), type(t) {}

  face_corner operator () (fc_type t)
    { return face_corner(face,vertex,t); }

  int face;
  int vertex; 
  fc_type type; 
}; 

bool operator == (face_corner const& a, face_corner const& b);
bool operator != (face_corner const& a, face_corner const& b);
bool operator < (face_corner const& a, face_corner const& b);
ostream& operator << (ostream& out, face_corner const& fc);


typedef map<face_corner, face_corner> fc_map; // We use this a lot. 
typedef face_corner face_vertex; // Use the name face_vertex when referring to cells. 
typedef vector<face_vertex> cell_face; 
typedef vector<cell_face> cell; 

class triangulation_builder {

  fc_map paired, adjacent, rev_adjacent, next; 
  bool ok; 

  vector<cell> CX;
  fc_map index;
  fc_map tfc_dual_edge; // Gives edge of combinatorial tube dual to a given tet face corner. 

public:

  triangulation_builder() : ok(true) {}

  void add_face_corner(face_corner const& f, face_corner const& p, face_corner const& a, face_corner const& n);

  void print_faces() const; 
  void glue_face(int n);
  bool get_triangulation_data(TriangulationData& TD, int report=0); 
  bool compute_cell_complex();
  void print_cell_complex() const;
  face_corner dual_edge(face_vertex const& fv) const;

  void split_all_cell_faces(); 
  void split_cell_face(face_vertex const& X, int offset);
  bool all_cell_faces_triangular() const; 
  void eliminate(face_corner const& e);
  void eliminate_edge(face_corner const& e);
  int vertex_order(face_corner const& fc);
  void remove_redundant_edges(int report);
  void remove_redundant_vertices();

private:

  void error(const char* msg);
  void find_tet_vertex(face_corner const& fc, int& tet, int& vertex) const; 

  face_vertex fv_paired(face_vertex const& fv) const;
  face_vertex fv_adjacent(face_vertex const& fv) const;
  face_vertex fv_next(face_vertex const& fv) const;
};

#endif
