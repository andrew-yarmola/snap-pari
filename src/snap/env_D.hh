#ifndef _env_D_
#define _env_D_
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

#include "env_T.hh"
#include "dd_tiling.hh"

class env_D : public env_T {
  static double dirichlet_epsilon;
  static unsigned int Dirichlet_print_options; 

public:
  static void setup_menu();

protected:

  void validate_event(int& what); 
  void process_event(int what); 
  virtual void print_settings() const; 

protected:
  env_D() : domain(0) {}
  virtual ~env_D(); 

  virtual env_D* get_D(int i) =0;

  // Dirichlet domain. 
  bool get_Dirichlet_domain();

  // Tiling. 
  double current_tiling_radius() const { return tiling.radius; }
  Tile* the_tiling() const { return tiling.root; }

  // Geodesics. 
  bool compute_geodesics(double max_len);
  bool compute_geodesics(int min_n);
  int num_geodesics() const { return geodesics.size(); }
  GeodesicWord const& geodesic(int i) const { return geodesics[i]; }
  bool find_core_geodesic(int index, Complex& clen, 
			  int& gn, int& ori, FGWord& conj) const;

  void print_lengths(double max_len=0.) const; 

  // Ortholines. 
  vector<Ortholine> ortholines(vector<int> const& gnums, double cutoff); 

  // Symmetries. 
  int compute_symmetries(double lmx=3.0, double dmax=2.0, int report=0);
  bool symmetry_orbits(vector<geo_orbit>& orbits, int n_orbits, int n_geod, double max_length, int info, int report=0);

  void print_symmetry_action(int gnum, const vector<O31_matrix>* mats=0, bool all_commens=false);
  void print_geodesic_symmetries(vector<int> const& gnums) const;

  void hsymm_geodesic_action();

  // 3-d graphics (Geomview file format)
  void gv_print_Dirichlet(FILE* fp, vector<int> const& gnums) const;
  void gv_print_Dirichlet(FILE* fp, int gnum, int n_ort);

  virtual void T_shape_changed();
  virtual void T_changed();
  void clear(); 
  void conjugate(O31_matrix const& c); 

  WEPolyhedron* domain;
  dd_tiling tiling;
  std::vector<GeodesicWord> geodesics;
  std::vector<O31_matrix> symmetries; 
};

Complex nearby_complex(Complex z, Complex const& z0);
bool words_to_ortholines(const GroupPresentation* group, Triangulation* manifold, vector<Ortholine>& ort, bool ac);
bool get_word(FGWord& w, vector<GeodesicWord> const& geodesics);

#endif
