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
#include "i_triangulation.hh"
#include "snappea/kernel.h"
#include "snappea/unix_io.h"
#include "kernel_extras.hh"
#include "orthoangles.hh"
#include <iostream>
#include "tilt_polytope.hh"

using std::cerr;

// static say xx("i_triangulation.cc"); 

i_triangulation::i_triangulation()
  : manifold(0), P(0)
{
  // fg[0] = fg[1] = 0; 
}

i_triangulation::i_triangulation(Triangulation *m)
  : manifold(m), P(0)
{
  // fg[0] = fg[1] = 0; 
}

i_triangulation::i_triangulation(i_triangulation const& m)
  : P(0)
{
  // fg[0] = fg[1] = 0; 
  copy_triangulation(m.manifold, &manifold);
}

i_triangulation::~i_triangulation()
{
  if (P) delete P; 
  clear_groups(); 
  if (manifold) free_triangulation(manifold); 
}

GroupPresentation* i_triangulation::group(int s) const
{
  int sol = get_filled_solution_type(manifold);

  if (sol==not_attempted || sol==degenerate_solution || sol==no_solution) 
    return 0; 

  return fundamental_group(manifold,s,s,s,TRUE);
}

void i_triangulation::position_vertices()
{
  choose_generators(manifold, TRUE, FALSE); // Computes vertices in the triangulation. 
}

void i_triangulation::clear_groups()
{
  // if (fg[1]) { free_group_presentation(fg[1]); fg[1] = 0; }
  // if (fg[0]) { free_group_presentation(fg[0]); fg[0] = 0; }
}

void i_triangulation::clear_shape_dependencies()
{ clear_groups(); }

void i_triangulation::triangulation_changing()
{ clear_groups(); }

void i_triangulation::clear()
{ 
  if (P) { delete P; P=0; }
  clear_groups(); 
  if (manifold) { free_triangulation(manifold); manifold = 0; }
}

void i_triangulation::set_manifold(Triangulation* m)
{
  if (manifold) clear(); 
  manifold = m;
}

#if 0
bool i_triangulation::read_census_manifold(const char *path, short c, short n)
{  
  set_manifold(::read_census_manifold(path, c, n));

  return manifold != 0; 
}

bool i_triangulation::read_manifold_file(const char *path, const char *name)
{
  set_manifold(::read_manifold_file(path, name)); 

  return manifold != 0; 
}
#endif

bool i_triangulation::save_manifold_file(const char* name) const
{
  FILE* fp = fopen(name, "w");
  if (!fp) return false;
  write_manifold_file(fp, manifold);
  fclose(fp);
  return true;
}

void i_triangulation::print_domain() const 
{
  int i, n = get_num_generators(manifold);
  MoebiusTransformation *gens = new MoebiusTransformation[n];
  matrix_generators(manifold, gens, FALSE, TRUE);

  print_tetrahedra(manifold); 
  printf("\n"); 
  for (i = 0; i < n; i++) {
    printf("%c=", i+'a');
    print_tetrahedra(gens[i], manifold);
    printf("\n"); 
  }

  delete [] gens; 
}

vector<Complex> i_triangulation::shapes() const
{
  vector<Complex> shapes;
  shapes.reserve(get_num_tetrahedra(manifold));

  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    {
      shapes.push_back(t->shape[filled]->cwl[ultimate][0].rect);
    }
  return shapes;
}

vector<double> i_triangulation::surgery_coeffs() const
{
  Cusp *cusp; 

  int n_cusps = manifold->num_cusps; 
  vector<double> coeffs(2*n_cusps); 

  for (cusp = manifold->cusp_list_begin.next;
       cusp != &manifold->cusp_list_end;
       cusp = cusp->next)
    {
      coeffs[2*cusp->index] = cusp->m; 
      coeffs[2*cusp->index + 1] = cusp->l; 
    }
  return coeffs; 
}

void i_triangulation::set_surgery_coeffs(const vector<double>& coeffs)
{
  int i, n = get_num_cusps(manifold);
  double m, l; 

  for (i=0; i<n; i++) {
    m = coeffs[2*i];
    l = coeffs[2*i+1];
    set_cusp_info(manifold, i, (m==0.0 && l==0.0), m, l);
  } 
}

// SolutionType is described in SnapPea.h. 

SolutionType i_triangulation::do_Dehn_surgery(const vector<double>& coeffs)
{
  clear_shape_dependencies(); 
  set_surgery_coeffs(coeffs); 
  return do_Dehn_filling(manifold);
}

SolutionType i_triangulation::find_complete_hyperbolic_structure()
{
  if (num_filled_cusps(manifold) > 0) clear_shape_dependencies();
  set_surgery_coeffs(vector<double>(2*manifold->num_cusps, 0.0)); 
  return ::find_complete_hyperbolic_structure(manifold);
}

SolutionType i_triangulation::do_Dehn_surgery(const vector<Complex>& holo, const vector<bool>& hism, vector<Complex>& tgt)
{
  Boolean *h_is_mer = new Boolean[hism.size()];
  copy(hism.begin(), hism.end(), h_is_mer); 
  Complex *holo_vec = new Complex[holo.size()];
  copy(holo.begin(), holo.end(), holo_vec);
  Complex *tgt_vec  = new Complex[tgt.size()];

  clear_shape_dependencies();
  SolutionType res = 
    solve_for_holonomies(manifold, holo_vec, h_is_mer, tgt_vec);

  copy(tgt_vec, tgt_vec+tgt.size(), tgt.begin()); 
  delete[] tgt_vec; 
  delete[] holo_vec; 
  delete[] h_is_mer;

  return res; 
}

void i_triangulation::make_orientable(int report)
{
  if (get_orientability(manifold) != nonorientable_manifold) return; 

  if (report) printf("Passing to orientable double cover.\n"); 
  Triangulation* dc = double_cover(manifold); 

  /* Adjust the name. */ 
  char name[100]; 
  strcpy(name, get_triangulation_name(manifold)); 
  strcat(name, ".or"); 
  set_triangulation_name(dc, name);

  set_manifold(dc); 
}

void i_triangulation::reorient()
{
  if (get_orientability(manifold) == nonorientable_manifold) {
    printf("Manifold is nonorientable.\n");
    return; 
  }

  triangulation_changing(); 
  ::reorient(manifold); 

  /* Adjust the name. */ 
  char name[100]; 
  strcpy(name+1, get_triangulation_name(manifold)); 
  if (name[1] == '-') {
    set_triangulation_name(manifold, name+2);
  } else {
    name[0] = '-';
    set_triangulation_name(manifold, name);
  }
}

bool i_triangulation::canonize(double area_multipliers[])
{
  if (num_filled_cusps(manifold) > 0) clear_shape_dependencies();
  triangulation_changing(); 
  if (proto_canonize(manifold, area_multipliers)!=func_OK) {
    clear(); 
    return false; 
  }
  find_complete_hyperbolic_structure();
  return true;
}

void i_triangulation::randomize()
{
  triangulation_changing(); 
  randomize_triangulation(manifold);
}

void i_triangulation::split_edge(int e, int f[2])
{
  triangulation_changing(); 
  if (!check_peripheral_curves(manifold)) {
    cerr << "Problem with peripheral curves\n";
    return;
  }
  EdgeClass *edge = find_edge(manifold, e); 
  if (!edge) { 
    cerr << "Edge " << e << " not found!\n";
    return; 
  }
  direct_edges(manifold); 
  EdgeClass *new_e;
  if (uncancel_tetrahedra(edge,f,&new_e,&manifold->num_tetrahedra)!=func_OK) {
    cerr << "Edge splitting failed\n";
  }
  if (!check_gluings(manifold)) {
    cerr << "Gluings messed up!\n";
    return;
  }
  if (!check_peripheral_curves(manifold)) {
    cerr << "Peripheral curves were corrupted!\n";
    return;
  }
  tidy_peripheral_curves(manifold); 
  orient_edge_classes(manifold); 
  number_the_edge_classes(manifold); 
  number_the_tetrahedra(manifold);
}

void i_triangulation::simplify()
{
  triangulation_changing(); 
  basic_simplification(manifold);
}

void i_triangulation::fill()
{
  Triangulation* filled = fill_reasonable_cusps(manifold); 
  if (!filled) return; 
  set_triangulation_name(filled, (char*)get_filling_name(manifold,2).c_str());
  triangulation_changing();
  free_triangulation(manifold); 
  manifold = filled; 
}

void i_triangulation::test_bc_subdivision()
{
  Triangulation* bc = barycentric_subdivision(manifold); 
  free_triangulation(bc); 
}

void i_triangulation::change_peripheral_curves(MatrixInt22 change_matrices[])
{
  if (::change_peripheral_curves(manifold, change_matrices)!=func_OK)
    printf("Attempt to change peripheral curves failed.\n"); 
}

bool i_triangulation::standardize_curves(Triangulation* ref)
{
  // Seek an isometry. 
  
  IsometryList *isometry_list = 0;
  Boolean isometric, or_pres; 
  int isom_num, attempt=0;
  while (true) { 
    compute_isometries(ref, manifold, &isometric,&isometry_list,NULL); 
    if (!isometric) {
      printf("manifold in standardize_curves is not isometric!\n");
      return false; 
    }
    if (!isometry_list) {
      printf("no isometry list found in standardize_curves!\n");
      return false; 
    }

    // Choose an isometry which is orientation preserving..
    or_pres = orientation_preserving_isometry(isometry_list, isom_num);
    if (or_pres || attempt > 0) break; 

    free_isometry_list(isometry_list); 
    isometry_list = 0;

    // Or reverse the orientation of manifold.
    reorient(); 
    attempt++; 
  }

  if (!or_pres) {
    printf("standardize_curves was unable to find an or.pres. isometry!\n"); 
    return false; 
  }

  // In an isometry list the complete cusps are numbered from zero, 
  // while filled cusps are discounted entirely. Therefore we make a 
  // table so we can go from an IsometryList cusp index
  // to a Triangulation cusp index. 


  // Start initializing the list of curve change matrices while 
  // we're at it. 

  int i, j = 0; 
  int c = get_num_cusps(manifold); 
  MatrixInt22 *pcb = new MatrixInt22[c]; 
  vector<int> cusp_num(c);
  for (i=0; i<c; i++) {
    if (cusp_is_complete(manifold,i)) {
      cusp_num[j++] = i; 
    } else {
      pcb[i][0][0] = 1; 
      pcb[i][0][1] = 0; 
      pcb[i][1][0] = 0; 
      pcb[i][1][1] = 1; 
    }
  }
  int n_cusps = j; 

  int mx[2][2];
  int img_cusp; 
  for (j=0; j < n_cusps; j++) {
    isometry_list_cusp_action(isometry_list, isom_num, j, &img_cusp, mx); 
    i = cusp_num[img_cusp];
    pcb[i][0][0] = mx[0][0]; 
    pcb[i][0][1] = mx[1][0]; // change matrix is transpose
    pcb[i][1][0] = mx[0][1]; 
    pcb[i][1][1] = mx[1][1]; 
  }

  // Change curves on manifold to agree with those on ref. 

  change_peripheral_curves(pcb);

  delete [] pcb; 
  free_isometry_list(isometry_list); 

  return true; 
}



double i_triangulation::volume(int *precision) const
{
  return ::volume(manifold, precision);  
}

double i_triangulation::chern_simons(int *precision, int *known) const
{
  Boolean is_known, requires_init;
  double value; 
  get_CS_value(manifold,&is_known,&value,precision,&requires_init); 
  *known = (int)is_known; 
  return value;
}

void i_triangulation::clear_tilt_polytope()
{
  delete P; 
  P = 0; 
}

bool i_triangulation::compute_tilt_polytope(vector<int> const& ties, 
					    n_vector const& tie_ratios, int limit)
{
  triangulation_changing(); 
  P = new tilt_polytope(manifold, ties, tie_ratios); 
  return P->compute(0, limit); 
}

#if 0
vector<Complex> i_triangulation::edge_orthodistances() const
{
  int n = get_num_tetrahedra(manifold); 
  vector<Complex> orth(n, Infinity);
  get_edge_orthodistances(manifold, group(0), orth); 
  return orth; 
}
#endif
