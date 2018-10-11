#ifndef _env_U_
#define _env_U_
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

#include "env_D.hh"
#include "tube.hh"

class env_U : public env_D {
  static string gv_picture_opts;
  static string picture_options;

  // TUBE RELATED STUFF.
  tube U; 

  // For a tube domain around a geodesic
  int g_num; 
  double tube_tile_rad; 
  Ortholine_set tube_ort, used_ort;

  // For tubes in manifold with Dehn surgery sing.
  Complex holonomy[2];
  vector<Ortholine> core_ort;

public:  
  static void setup_menu();

protected:
  env_U() : g_num(-1), tube_tile_rad(0.) {}
  virtual ~env_U() {}

  void validate_event(int& what); 
  void process_event(int what); 

  virtual void print_settings() const; 

  virtual void T_shape_changed(); 
  virtual void T_changed();

  void print_peripheral_curve_info(int cusp, double cm, double cl);

  // MORE TUBE STUFF
  const tube& the_tube() const { return U; }
  int geodesic_number() const { return g_num; }
  bool tube_is_valid() const { return g_num!=-1; }
  void clear_tube();

  int compute_tube(int gn, bool report=false); 
  void drill_tube(int manifold_num);

  bool recheck_face_pairings();

  vector<Ortholine> get_core_ortholines() const; 
  bool compute_core_ortholines(); 
  int num_core_ortholines() const { return core_ort.size(); }
  void compute_core_tube(bool report=false);
  int verify_core_tube(); 
  void print_core_words() const;
  void add_core_word(FGWord const& w);
  void delete_core_word(int n);
  void delete_core_word(FGWord const& w);
  void print_clip_list() const; 
  double core_tube_radius() const; 
  void print_tube_ortholines() const;
  bool ortholine_to_word(Ortholine const& o, FGWord& wd) const;

  void save_picture(string const& file, string const& options) const;

  // Tweaking. 
  void transpose_face(int n) { U.transpose_face(n); }
  void transpose_no_faces() { U.transpose_no_faces(); }
};

#endif
