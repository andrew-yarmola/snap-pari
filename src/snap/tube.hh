#ifndef TUBE_H
#define TUBE_H
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

#include "ortholine.hh"
#include "tube_face.hh"
#include "picfile.hh"

class connected_faces;
class clipping_cache; 

class tube {

  vector<Complex> H; // Holonomies (lll reduced)

  vector<tube_face> faces; // Two per ortholine. 
  vector<Ortholine> ortholines; 

  connected_faces* CF;
  clipping_cache* CC;

public:
  static double ideal_cutoff;

public:

  tube();
  ~tube();

  void clear(); 

  void set_holonomies(Complex const& me, Complex const& lo);

  int compute_tube(Ortholine_set& os, Ortholine_set& used, double rigorous_to_rad, bool report=false);
  int compute_tube(vector<Ortholine> const& ort, bool report=false);
  int compute_tube(Triangulation* T, GroupPresentation* G, bool report=false);

  bool check_face_pairings(double& urad, int report = 0) const; 
  int num_faces() const; 
  bool tidy_edges(bool report);

  void get_holonomies(Complex& m, Complex& l) const
    { m = H[0]; l = H[1]; }
  vector<Ortholine> const& the_ortholines() const 
    { return ortholines; }

  void print_ortholines() const; 

  // Debugging.
  void print_clip_list(vector<Ortholine> const& ort) const;

  // Tweaking. 
  void transpose_face(int n);
  void transpose_no_faces();

  void save_picture(const string& print_options, picfile& pic) const;
  void save_natural_face_picture(picfile& pic) const;

  void print_faces(bool natural) const;
  void print_face_matrices(int n) const; 

  // connected faces & drilling
  bool compute_connected_faces();

  void print_connected_faces() const; 
  void print_face_corners() const;

  Triangulation* drill(); 
  void print_triangulation_data();

  // These functions temporarily disabled. 
  void save_facepairing_picture(picfile& pic) const {}
  void covered_curve_picture(picfile& pic, Complex const& curve) const {}
  void save_gv_picture(string const& opts) const {}

  Complex holonomy(Mark const& mark) const; 

  // Not used anywhere
  void add_unclipped_face(Ortholine const& ol); 
  void sort_edges();

private:
  int add_face(Ortholine const& ol, bool report=false);

};

int triangulation_is_tubal(Triangulation* T, GroupPresentation* G);

#endif
