#ifndef _snap_
#define _snap_
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


#include <stdio.h>
#include <list>
#include "i_triangulation.hh"
#include "pariclass.hh"

class snap : public i_triangulation {
public:

  snap();
  snap(Triangulation *m);

  virtual ~snap();

  void compute_accurate_info(int polish=0, int report = 0);
  void update_precision(); 

  int accurate_info_set() const { return _acc_shapes.type()!=1; }
  const pari& acc_shapes() const { return _acc_shapes; }

  // SolutionType do_Dehn_surgery(const vector<double>& coeffs);
  // SolutionType find_complete_hyperbolic_structure(); 

  virtual GroupPresentation* group(int simplified) const;

  void reorient();
  bool canonize(double area_multipliers[]=0); 
  void simplify();
  void randomize();
  void fill(); 
  void test_bc_subdivision(); 

  virtual void change_peripheral_curves(MatrixInt22 change_matrices[]); 

  virtual void print() const; 

  void print_homology_matrix() const; 
  void print_npz2_homology() const; 
  pari cusp_homology_kernel() const;

  pari accurate_log_shapes() const; 
  pari acc_holonomies() const; 
  pari core_geodesics() const; 

  pari dilog_sum() const; 
  pari complex_volume() const; 
  pari complex_volume2() const; 
  pari accurate_volume() const; 
  pari accurate_chern_simons(int& correct) const; 
  pari raw_chern_simons() const; 

  pari eta_invariant() const;
  pari eta_fudge() const; 
  pari eta_difference() const; 
  void check_eta() const; 
  bool eta_available() const; 
  bool set_eta_info(); 
  void set_eta_invariant(pari x); 
  void set_eta_fudge(pari x); 
  void print_eta_info() const; 

  void find_a_flattening(); 
  void set_flattening(pari f); 
  pari the_flattening() const { return flattening; }

  bool bootstrap_CS_value(double& value, int report = 0); 

  std::string filling_name() const;

  int is_closed() const; 
  int is_complete() const;

protected:

  virtual void clear();
  virtual void triangulation_changing();
  virtual void clear_shape_dependencies(); 

  void set_acc_shapes(const pari& new_shapes) { _acc_shapes = new_shapes; }

private: 

  pari _acc_shapes; 

  bool _eta_invariant_known; 
  pari _eta_invariant;

  bool _eta_fudge_known; 
  pari _eta_fudge; 
  pari flattening; 
};

#include "int_matrix.hh"
int_matrix to_int_matrix(pari const& m);
void show_matprod(pari const& mat, pari const& v);
void complex_volume3(Triangulation* T, vector<int> const& e_or);

/* numeric.cc */ 

int dedekind(int p, int q); 
pari eta_rat_diff(pari cusp_hom_ker, pari surgery_coeffs, pari cusp_maps); 

pari factor_sl2Z_matrix(pari mx); 
pari cusp_term_adjustment(pari mx, pari pq); 

pari cusp_homology_kernel(pari chr, pari ehr);
pari cusp_homology_kernel(pari chr, pari ehr, pari sc);

pari signature_Y(pari k0, pari sc);

pari monomial(pari shapes, pari powers); 

#endif

