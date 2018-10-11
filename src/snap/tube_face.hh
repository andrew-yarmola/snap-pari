#ifndef _tube_face_
#define _tube_face_
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

#include "eqs_interval.hh"
#include "eq_surface.hh"
#include "uc_polylist.hh"
#include "bbox.hh"

class color; 
class picfile;

class tube_face {
  eqsfunc surface; 
  interval_list edges; 
  Mark mark_num; 

  // da_spec's of all clippings applied to this face. 
  vector<da_spec> clipping_surfaces;
  
  // Bounds.
  double rad; 
  bbox box;
public:
  static double epsilon;
  static double boundary_epsilon;

  tube_face() {}
  tube_face(eqsfunc const& G, Mark mark_num, bool transpose=false); 

  int clip(eqsfunc const& G, Mark mark_num, bool ucover=true); 
  // int clip(tube_face const& f, bool transpose, bool ucover); 
  
  int make_face(da_spec const& face, vector<da_spec> const& nbrs, bool report=false);
  bool validate(vector<da_spec> const& nbrs, bool report=false); 

  bool is_empty() const { return edges.L.size()==0; }
  bool has_edge_at_infinity() const { return edges.has_mark(Mark()); }
  void clean_up();
  void translate(Complex const& z); 

  void set_to_partner(tube_face const& f, Complex const& oa); 
  void sort_edges(); 
  uc_polylist get_boundary() const; 

  void add_clipping_intervals(eqsfunc const& G, bool ucover);
  // void save_last_clipping_matrix(ostream& out) const; 

  void picfile_print(picfile& pic, Complex const& center, double size) const; 
  void picfile_print(picfile& pic, color const& col) const; 
  void picfile_print(picfile& pic, color const& col, const list<Complex>& trans) const;

  void uc_print(ostream& out, vector<Complex> const& H) const; 
  void print_clippings(ostream& out) const; 
  void print_edge_markings() const; 

  bool has_matching(eqs_const_iter const& e, Mark const& mk, vector<Complex> const& H) const;
  eqs_const_iter begin() const { return edges.L.begin(); } 
  eqs_const_iter end()   const { return edges.L.end(); } 
  eqs_iter begin() { return edges.L.begin(); } 
  eqs_iter end()   { return edges.L.end(); } 
  void erase(eqs_iter& e) { edges.L.erase(e++); }
  Mark mark() const { return mark_num; }

  friend ostream& operator << (ostream& out, tube_face const& f); 
  friend bool faces_match(tube_face const& a, tube_face const& b, double& urad, double eps, int report);
  friend bool verify_tube_face(Complex const& d0, vector<Complex> const& da_spec, tube_face*& tf); 
};

#endif
