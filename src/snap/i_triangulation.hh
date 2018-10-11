#ifndef _i_triangulation_
#define _i_triangulation_
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
#include "n_vector.hh"

/* A simple wrapper class for Triangulations so we can conveniently
   use it as a base class for ideal triangulations with other
   stuff thrown in. This code depends purely on SnapPea and C++. 
*/ 

using std::vector;

class tilt_polytope;

class i_triangulation {
public:
  i_triangulation();
  i_triangulation(i_triangulation const& m);
  i_triangulation(Triangulation *m);
  virtual ~i_triangulation();

  // Basic inspectors. 
  Triangulation* M() const { return manifold; }
  vector<Complex> shapes() const;
  vector<double> surgery_coeffs() const;

  // Input and Output. 
  bool save_manifold_file(const char* name) const;
  void print_domain() const; 

  // Basic modifier.
  void set_manifold(Triangulation* m);

  // Dehn surgery. 
  virtual SolutionType do_Dehn_surgery(const vector<double>& coeffs);
  virtual SolutionType find_complete_hyperbolic_structure();
  virtual SolutionType do_Dehn_surgery(const vector<Complex>& holo, const vector<bool>& hism, vector<Complex>& tgt);

  // A few SnapPea functions. 
  void make_orientable(int report=0);
  virtual void reorient();
  virtual bool canonize(double area_multipliers[]=0); 
  virtual void simplify();
  virtual void randomize();
  virtual void split_edge(int e, int f[2]);
  virtual void fill(); 
  virtual void change_peripheral_curves(MatrixInt22 change_matrices[]);
  virtual void test_bc_subdivision();
  bool standardize_curves(Triangulation* ref);

  // Fundamental group. 
  virtual GroupPresentation* group(int simplified) const;

  // SnapPea 3-manifold invariants. 
  double volume(int* precision = 0) const; 
  double chern_simons(int *precision, int *known) const; 

  // The orthodistance invariant
  // vector<Complex> edge_orthodistances() const;

  // Tilt polytope
  bool compute_tilt_polytope(vector<int> const& ties, n_vector const& tie_ratios, int limit=1000);
  void clear_tilt_polytope(); 
  tilt_polytope* tp() const { return P; } 

protected:
  void position_vertices(); 
  void set_surgery_coeffs(const vector<double>& coeffs);
  void clear_groups(); 

  virtual void clear();
  virtual void triangulation_changing();
  virtual void clear_shape_dependencies(); 

protected:
  Triangulation* manifold;

private:
  // GroupPresentation* fg[2]; 
  tilt_polytope* P; 
};

#endif
