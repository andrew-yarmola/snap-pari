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
#include "drill_Dirichlet.hh"

using std::cout;
using std::endl;

void lookup_index(const WEPolyhedron* P, const WEFace* f, const WEEdge* e, 
		  int& f_index, int& e_index)
{
  const WEFace *face;
  f_index = 0;
  const WEEdge *edge;
  bool fwd; 

  for (face = P->face_list_begin.next;
       face != &P->face_list_end;
       face = face->next, f_index++) {

    if (f!=face) continue; 

    edge = face->some_edge;
    for (e_index = 0; e_index < face->num_sides; e_index++) {

      if (e==edge) return; 

      // Find the next edge around this face. 

      fwd = (edge->f[left] == face);
      edge = edge->e[fwd ? tip:tail][fwd ? left:right]; 

    }
  }
  cout << "Warning: lookup_index failed.\n";
}


void lookup_index(const WEPolyhedron* P, const WEFace* f, const O31_vector& v, 
		  int& f_index, int& e_index)
{
  const WEFace *face;
  f_index = 0;
  const WEEdge *edge;
  bool fwd; 

  for (face = P->face_list_begin.next;
       face != &P->face_list_end;
       face = face->next, f_index++) {

    if (f!=face) continue; 

    edge = face->some_edge;
    for (e_index = 0; e_index < face->num_sides; e_index++) {

      fwd = (edge->f[left] == face);

      if (close(edge->v[fwd ? tip:tail]->x, v, 1e-6)) return; 

      // Find the next edge around this face. 

      edge = edge->e[fwd ? tip:tail][fwd ? left:right]; 

    }
  }
  cout << "Warning: lookup_index failed.\n";
}


void print_dirichlet_info(const WEPolyhedron* P)
{
  triangulation_builder T;
  get_dirichlet_info(P,T);
  T.print_faces(); 
}

void get_dirichlet_info(const WEPolyhedron* P, triangulation_builder& T)
{
  const WEFace *face, *fa;
  int f_index = 0, fa_index, fp_index; 
  const WEEdge *edge, *ea;
  int e_index, ea_index, ep_index, en_index; 
  bool fwd; 
  O31_matrix *group_elt; 
  O31_vector paired_vertex;

  for (face = P->face_list_begin.next;
       face != &P->face_list_end;
       face = face->next, f_index++) {

    edge = face->some_edge;
    for (e_index = 0; e_index < face->num_sides; e_index++) {

      fwd = (edge->f[left] == face);

      // Find the adjacent (face, edge) pair. 

      fa = edge->f[fwd ? right:left];
      ea = edge->e[fwd ? tip:tail][fwd ? right:left];
      lookup_index(P, fa, ea, fa_index, ea_index); 

      // Find the paired corner. 

      group_elt = face->mate->group_element;
      paired_vertex = *group_elt * edge->v[fwd ? tip:tail]->x;
      paired_vertex /= paired_vertex[0]; 

      lookup_index(P, face->mate, paired_vertex, fp_index, ep_index);
      
      // Find the next edge around this face. 

      edge = edge->e[fwd ? tip:tail][fwd ? left:right]; 
      en_index = e_index+1;
      if (en_index==face->num_sides) en_index = 0; 

      // Insert into the triangulation builder. 

      T.add_face_corner(face_corner(f_index, e_index), face_corner(fp_index, ep_index), 
			face_corner(fa_index,ea_index),face_corner(f_index,  en_index));

    }
  }
}

