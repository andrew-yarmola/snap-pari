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
#include "get_equations.hh"
#include "kernel_extras.hh"
#include "snappea/kernel.h" 
#include "pariwrap.hh"
#include <vector>
#include <math.h>

#ifdef PARI_2_2_OR_LATER
#include "pari_oldnames.hh"
#endif

using std::cout;

// static say xx("get_equations.cc"); 

/* Returns a pari column vector containing z1,...,zn. */ 
pari tet_shapes(Triangulation *manifold)
{
  int n;

  n = get_num_tetrahedra(manifold);

  pari shapes = cvector(n); 

  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    {
      shapes[t->index] = to_pari(t->shape[filled]->cwl[ultimate][0].rect);
    }

  return shapes; 
}

/* Returns log(z1),...,log(zn),log(1-z1),...,log(1-zn) */ 
pari log_shapes(Triangulation *manifold)
{
  int n;

  n = get_num_tetrahedra(manifold);

  pari shapes = cvector(2*n); 

  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    {
      shapes[t->index] = to_pari(t->shape[filled]->cwl[ultimate][0].log);
      shapes[t->index+n] = -to_pari(t->shape[filled]->cwl[ultimate][1].log);
    }

  return shapes; 
}

pari cusp_shapes(const Triangulation* manifold)
{
  int i, n_cusps = get_num_cusps((Triangulation*)manifold); 
  Cusp *cusp; 

  pari shapes = rvector(n_cusps);

  for (i = 0; i < n_cusps; i++) {
    cusp = find_cusp((Triangulation*)manifold, i);
    if (!cusp->is_complete) continue; // entry will be pZERO
    shapes[i] = to_pari(cusp->cusp_shape[current]); 
  }

  return shapes; 
}

pari holonomies(Triangulation *manifold)
{
  pari shapes = concat(log_shapes(manifold), pPI * pI); 
  return cusp_equations(manifold) * shapes; 
}

#if 0
pari log_with_history(pari z, ShapeInversion *history)
{
  int n=0, down=1;

  while (history) {
    switch (history->wide_angle) {
    case 0: /* crossing [-inf,0] */
      n += down; 
      break;
    default:
      break;
    }
    down = -down; 
    history = history->next; 
  }
  return log(z) + pTWO_PI * gi * integer(n); 
}

pari log_one_minus_z_with_history(pari z, ShapeInversion *history)
{
  int n=0, down=1;

  while (history) {
    switch (history->wide_angle) {
    case 1: /* crossing [1,+inf] */
      n -= down; 
      break;
    default:
      break;
    }
    down = -down; 
    history = history->next; 
  }
  return log(pONE-z) + pTWO_PI * gi * integer(n); 
}
#endif

void clear_accurate_shapes(Triangulation* manifold)
{
  Tetrahedron *t; 
  
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    {
      t->acc_shape   = pZERO; 
      t->acc_corners = pZERO;
    }
}

void set_accurate_shapes(Triangulation* manifold, const pari& shapes)
{
  Tetrahedron *t; 

  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    {
      t->acc_shape = shapes[t->index];
    }
}

#if 0
pari get_accurate_cusp_shapes(const Triangulation *manifold, const pari& tet_shapes)
{
  set_accurate_shapes((Triangulation*)manifold, tet_shapes); 

  compute_cusp_shapes((Triangulation*)manifold, filled);
  int i, n_cusps = get_num_cusps((Triangulation*)manifold); 

  pari cusp_shapes = cvector(n_cusps); 
  Cusp *cusp;
  for (i=0; i<n_cusps; i++) {
    cusp = find_cusp((Triangulation *)manifold, i);
    cusp_shapes[i] = cusp->acc_cusp_shape;
  }

  /* leaving polymod shapes in the tetrahedron structure will cause
     problems if get_generators is ever called */ 
  clear_accurate_shapes((Triangulation*)manifold);

  return cusp_shapes; 
}
#endif

// Must call set_accurate_shapes before get_accurate_cusp_shapes and 
// clear_accurate_shapes afterwards. 

pari get_accurate_cusp_shapes(const Triangulation *manifold)
{
  compute_cusp_shapes((Triangulation*)manifold, filled);
  int i, n_cusps = get_num_cusps((Triangulation*)manifold); 

  pari cusp_shapes = cvector(n_cusps); 
  Cusp *cusp;
  for (i=0; i<n_cusps; i++) {
    cusp = find_cusp((Triangulation *)manifold, i);
    if (!cusp->is_complete) continue; 
    cusp_shapes[i] = cusp->acc_cusp_shape;
  }

  return cusp_shapes; 
}

pari surgery_coeffs(Triangulation* manifold)
{
  Cusp *cusp; 

  int n_cusps = manifold->num_cusps; 
  pari coeffs = cvector(3 * n_cusps); 
  int m, l, d; 

  for (cusp = manifold->cusp_list_begin.next;
       cusp != &manifold->cusp_list_end;
       cusp = cusp->next)
    {
      if (cusp->is_complete) {
	// For complete cusp, (m, l, denominator) = (1,0,0).
	coeffs[cusp->index] = pONE; 
      } else {
	ratsc(cusp->m, cusp->l, m, l, d); 
	coeffs[cusp->index] = integer(m); 
	coeffs[cusp->index + n_cusps] = integer(l);
	coeffs[cusp->index + 2*n_cusps] = integer(d); 
      }
    }
  return coeffs; 
}

pari surgery_coeffs2(Triangulation* manifold)
{
  Cusp *cusp; 

  int n_cusps = manifold->num_cusps; 
  pari coeffs = cvector(2 * n_cusps); 
  int m, l, d; 
  bool valid = true; 

  for (cusp = manifold->cusp_list_begin.next;
       cusp != &manifold->cusp_list_end;
       cusp = cusp->next)
    {
      if (cusp->is_complete) {
	// For complete cusp, (m, l, denominator) = (1,0,0).
	coeffs[cusp->index] = pZERO; 
	coeffs[cusp->index + n_cusps] = pZERO; 
      } else {
	ratsc(cusp->m, cusp->l, m, l, d); 
	coeffs[cusp->index] = integer(m); 
	coeffs[cusp->index + n_cusps] = integer(l);
	if (d!=1) valid = false;
      }
    }
  if (!valid) warn("this computation is invalid for cone manifolds\n"); 
  return coeffs; 
}

pari get_filled_equations(Triangulation *manifold, bool with_signs, bool full)
{
  int i,n = manifold->num_cusps, n_tet = get_num_tetrahedra(manifold);
  pari coeff_mx = matrix(2*n + 1,n); 
  pari coeffs = surgery_coeffs(manifold); 

  for (i=0; i<n; i++) {
    coeff_mx[i][i] = coeffs[i]; 
    coeff_mx[i+n][i] = coeffs[i+n]; 
    coeff_mx[2*n][i] = -pTWO * coeffs[i+2*n]; 
  }

  // have to work with the equations as columns to use concat.
  
  int nr = (full) ? (3*n_tet) : (2*n_tet+1); 
  pari pi_i_eqn = matrix(1,nr); 
  if (!full) pi_i_eqn[0][2*n_tet] = pONE; 
  pari cusp_eqns = concat(gtrans(cusp_equations(manifold,full)),pi_i_eqn); 
  
  /* concatenate the two sets of equations */ 
  return concat(cusp_eqns * gtrans(coeff_mx), 
		gtrans(edge_equations(manifold, with_signs, full)));
}

pari get_complete_equations(Triangulation *manifold, bool with_signs, bool full)
{
  /* concatenate the two sets of equations */ 
  return concat(gtrans(cusp_equations(manifold, full)), 
		gtrans(edge_equations(manifold, with_signs, full))); 
}

pari edge_equations(Triangulation* manifold, bool with_signs, bool full)
{
  Tetrahedron *tet;
  EdgeIndex e;
  int i, r, etype; 

  /* Number tetrahedra and edge classes. */
  number_the_tetrahedra(manifold);
  number_the_edge_classes(manifold);

  int n_tet = get_num_tetrahedra(manifold); 
  int n_edges = manifold->edge_list_end.prev->index + 1;
  int n_cols = (full) ? (3*n_tet) : (2*n_tet+1); 

  pari equations = matrix(n_cols, n_edges); 

  /* pari matrices are defined and accessed as (columns, rows).
     pari matrices are automatically initialized to (integer) 0. 

     Equations expresses the gluing equations in the following form:

     Full: 
     Each row contains 3*n entries, [a1,b1,c1,...,an,bn,cn]
     representing the edge product z1^a1.(1/(1-z1))^b1.(1-1/z1)^c1....

     Normal: 
     Each row contains 2*n+1 entries, [a1,...,an,b1,...,bn,c]
     representing the equation z1^a1...zn^an.(1-z1)^b1.(1-zn)^bn = e^(pi i c).
  */

  /* Can't do anything if it's not an oriented manifold */ 
  if (manifold->orientability != oriented_manifold) return equations; 

  pari unit; 
    
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {

    i = tet->index; 

    unit = integer(((!with_signs)||tet->has_correct_orientation) ? 1 : -1); 

    for (e = 0; e < 6; e++) {   /* Look at each of the six edges. */

      r = tet->edge_class[e]->index;

      /*
       * edge3[e] = 0, 1, 2 according to whether the contribution
       * to the equation for this edge is z, 1/(1-z), (z-1)/z. 
       * (where z is the parameter for this tetrahedron). 
       */
      etype = edge3[e];

      if (full) {
	equations[3*i + etype][r] += unit;
      } else {
	if (etype == 0) 
	  equations[i][r] += unit;
	else if (etype == 1) 
	  equations[i+n_tet][r] -= unit;
	else {
	  equations[i][r] -= unit;
	  equations[i+n_tet][r] += unit;
	  equations[2*n_tet][r] += unit;
	}
      }	
    }
  } 

  if (!full) {
    /* Angle sum around edges is 2 PI */ 
    for (r = 0; r < n_edges; r++)
      equations[2*n_tet][r] -= pTWO; 
  }

  return equations; 
}


pari cusp_equations(Triangulation *manifold, bool full)
{
  Tetrahedron *tet;
  VertexIndex v;
  FaceIndex initial_side, terminal_side;
  int init[2], term[2];
  int i, r, contrib, ml; 
  int etype; 

  int n_tet = get_num_tetrahedra(manifold); 
  int n_cusps = manifold->num_cusps; 
  int n_cols = full ? (3*n_tet) : (2*n_tet+1); 

  pari equations = matrix(n_cols, 2*n_cusps); 
    

  /* pari matrices are automatically initialized to (integer) 0. 

     Equations expresses the cusp equations in the following form:

     Full: 
     Each row contains 3*n entries, [a1,b1,c1,...,an,bn,cn]
     representing the meridian or longitude product 
     z1^a1.(1/(1-z1))^b1.(1-1/z1)^c1....

     Normal:
     Each row contains 2*n+1 entries, [a1,...,an,b1,...,bn,c]
     representing the equation z1^a1...zn^an.(1-z1)^b1.(1-zn)^bn = e^(pi i c).
  */

  /* Can't do anything if it's not an oriented manifold */ 
  if (manifold->orientability != oriented_manifold) return equations; 
    
  /* Number tetrahedra and edge classes */
  number_the_tetrahedra(manifold);

  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {

    i = tet->index; 

    for (v = 0; v < 4; v++) {   /* Look at each ideal vertex. */
	
      if (!tet->cusp[v]) continue; // Ignore finite vertices. 

      /*
       *  Each ideal vertex contains two triangular cross sections,
       *  one right_handed and the other left_handed. (Since the manifold
       *  is oriented we consider only the right_handed cross section.)
       *  We want to compute the contribution of each angle of 
       *  each triangle to the holonomy. A directed
       *  angle is specified by its initial and terminal sides.
       */
      for (initial_side = 0; initial_side < 4; initial_side++) {
	if (initial_side == v) continue;
	terminal_side = remaining_face[v][initial_side];
	  
	for (ml = 0; ml < 2; ml++)  {       /* meridian,longitude */
	  /*
	   *  Note the intersection numbers of the meridian and
	   *  longitude with the initial and terminal sides (on 
	   *  the right-handed triangular cross section).
	   */
	  init[ml] = tet->curve[ml][0][v][initial_side];
	  term[ml] = tet->curve[ml][0][v][terminal_side];

	  contrib = (int)FLOW(init[ml],term[ml]);

	  etype = edge3_between_faces[initial_side][terminal_side];

	  if (full) {
	    r = 2 * tet->cusp[v]->index + ml;
	    equations[3*i + etype][r] += integer(contrib); 

	  } else {

	    r = tet->cusp[v]->index + ml*n_cusps;       

	    if (etype==0) {
	      equations[i][r] += contrib; 
	    } else if (etype==1) {
	      equations[i+n_tet][r] -= contrib;
	    } else {
	      equations[i][r] -= contrib;
	      equations[i+n_tet][r] += contrib;
	      equations[2*n_tet][r] += contrib;
	    }
	  }
	}
      } 
    } 
  }

  return equations; 
}

pari edge_homology_relations(Triangulation *manifold)
{
  pari mx = matrix(get_num_tetrahedra(manifold), manifold->num_generators); 

  /* get the relations coming from the edges */ 
  EdgeClass *edge;
  PositionedTet ptet, ptet0;
  int col = 0;

  for (edge = manifold->edge_list_begin.next;
       edge != &manifold->edge_list_end;
       edge = edge->next) {
    set_left_edge(edge, &ptet0);
    ptet = ptet0;
    do {

      switch (ptet.tet->generator_status[ptet.near_face]) {
      case outbound_generator:
	mx[col][ptet.tet->generator_index[ptet.near_face]] += pONE;
	break;

      case inbound_generator:
	mx[col][ptet.tet->generator_index[ptet.near_face]] -= pONE;
	break;

      default:
	break;
      }

      veer_left(&ptet);
    } while ( ! same_positioned_tet(&ptet, &ptet0) );

    col++;
  }

  return mx;
}

pari cusp_homology_relations(Triangulation *manifold)
{
  pari mx = matrix(2 * manifold->num_cusps, manifold->num_generators); 

  /* get meridian and longitude of each cusp */ 
  Tetrahedron *tet; 
  VertexIndex vertex;
  FaceIndex   side;
  Orientation orientation;

  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {
    for (vertex = 0; vertex < 4; vertex++) {
      for (side = 0; side < 4; side++) {
	if (side == vertex ||
	    tet->generator_status[side] != inbound_generator) continue;

	for (orientation = 0; orientation < 2; orientation++) { 
	  /* orientation = right_handed, left_handed */

	  /* meridian and longitude curves */
	  mx[tet->cusp[vertex]->index][tet->generator_index[side]]
	      += integer(tet->curve[M][orientation][vertex][side]);
	  mx[tet->cusp[vertex]->index + manifold->num_cusps][tet->generator_index[side]]
	      += integer(tet->curve[L][orientation][vertex][side]);

	}
      }
    }
  }
  return mx; 
}

struct dual_path_link {
  Tetrahedron* tet;
  FaceIndex prev, next;
};

static void get_path(vector<dual_path_link>& path, Tetrahedron* tet, FaceIndex prev)
{
  path.resize(0);

  FaceIndex next; 
  dual_path_link dpl; 
  int limit = 5000;
  while (--limit > 0) {
    dpl.tet = tet; 
    dpl.prev = prev; 
    next = dpl.next = tet->generator_path; 
    // cout << dpl.tet << ' ' << int(dpl.next) << endl; 
    path.push_back(dpl); 
    if (next == -1) break; 
    prev = EVALUATE(tet->gluing[next],next); // get the glued face on next tet.
    tet = tet->neighbor[next]; 
  }
}

//  Must call choose_generators(manifold, ...) before this function

pari get_parity_matrix(Triangulation* manifold)
{
  int N_gens = get_num_generators(manifold);
  int N = get_num_tetrahedra(manifold); 

  pari parity = matrix(2*N+1, N_gens); // All zeros initially.

  int i;
  Tetrahedron* tet;
  FaceIndex f;
  EdgeIndex e; 
  int j1, j2, n, ti; 
  vector<dual_path_link> P1, P2; 
  const bool print_path=false; 

  for (i=0; i<N_gens; i++) {

    // Find a tet,face pair for generator i.
    for (tet = manifold->tet_list_begin.next;
	 tet!= &manifold->tet_list_end;
	 tet = tet->next) {
      for (f=0; f<4; f++) {
	if (tet->generator_status[f] != not_a_generator &&
	    tet->generator_index[f]  == i) break; // Found it!
      }
      if (f < 4) break; // Found one. 
    }
    if (f==4) { printf("No face for generator %d in get_parity_matrix\n", i);
      return pZERO; } // This should never happen. 

    // Find paths to either side of this face.
    get_path(P1, tet, f);
    get_path(P2, tet->neighbor[f], EVALUATE(tet->gluing[f],f)); 

    // Find the point where the two paths join. 
    j1 = P1.size()-1; j2 = P2.size()-1;
    while (j1 >= 0 && j2 >= 0 && P1[j1].tet==P2[j2].tet) { j1--; j2--; }

    if (j1==P1.size()-1) { printf("Problem in get_parity_matrix\n"); 
      return pZERO; } // This should never happen. 
    
    // Throw away the common end parts, make P1 turn to join P2 reversed. 
    P1.resize(j1+2); 
    P1[j1+1].next = P2[j2+1].prev; 
    for (; j2 >= 0; j2--) P1.push_back(P2[j2]); // Make a loop. 
    // Strictly speaking, should reverse prev and next in P2 but 
    // for present purposes this doesn't matter. 

    // Add up the parity contributions. 
    n = P1.size();
    for (j1=0; j1<n; j1++) {
      e = edge3_between_faces[P1[j1].prev][P1[j1].next];
      ti = P1[j1].tet->index; 

      // print the dual path
      if (print_path) 
	printf("(%d;%d,%d,e%d) ", ti, P1[j1].prev, P1[j1].next, e);

      if (e==9) { printf("Path in get_parity_matrix has backtracking!\n"); 
	return pZERO; } // This should never happen. 

      // This is the parity contribution from Walter's paper.
      if (e==0) { parity[ti][i] += pONE; }
      else if (e==1) { parity[ti+N][i] -= pONE; }
      else {
	parity[ti][i] -= pONE; 
	parity[ti+N][i] += pONE; 
	parity[2*N][i] += pONE; 
      }
    }
    if (print_path) printf("%c\n", 'a'+i); 
  }

  return parity; 
}

// These are log_shapes without any special 
// choice of branch. We put a Pi*I at the end so that
// we can do equations * log_shapes to check them. 

pari p_log_shapes(Triangulation* manifold) 
{
  Tetrahedron* tet; 
  int i, n = get_num_tetrahedra(manifold); 
  pari l_shapes = cvector(2*n + 1); 

  static pari pPII = pPI*pI;
  
  for (tet = manifold->tet_list_begin.next;
       tet!= &manifold->tet_list_end;
       tet = tet->next) {
    i = tet->index;
    l_shapes[i] = log(tet->acc_shape); 
    l_shapes[i+n] = log(pONE - tet->acc_shape); 
  }

  l_shapes[2*n] = pPII; 

  return l_shapes; 
}

void check_equations(Triangulation* manifold, bool with_signs)
{
  pari res = gtrans(get_filled_equations(manifold, with_signs)) * 
    p_log_shapes(manifold); 
  printf("Filled equation results/Pi: "); 
  (res/pPI).print(); 
}

pari cross_ratio(pari const& z);

void check_shapes(Triangulation* manifold)
{
  Tetrahedron* tet; 
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {
    tet->acc_shape.print();
    cross_ratio(tet->acc_corners).print();
  }
}

static bool find_smallest_in_last_row(pari& mat, int& col)
{
  int c, cols = length(mat);
  int row = length(mat[0])-1;
  pari av, min_abs; 

  // find first nonzero. 
  for (c=0; c<cols; ++c)
    if (mat[c][row] != pZERO) break; 
  if (c==cols) return false; // all zero. 
  min_abs = abs(mat[c][row]);
  col = c; 

  // now find the best nonzero one. 
  for (; c<cols; ++c) {
    av = abs(mat[c][row]);
    if (av==pZERO) continue; 
    if (av < min_abs) {
      min_abs = av; 
      col = c; 
    }
  }

  return true; 
}

static bool reduce_other_cols(pari& mat, int col)
{
  int c, cols = length(mat); 
  int row = length(mat[0])-1;

  bool changed = false; 

  pari quot;
  pari mcn, smallest = mat[col][row];

  for (c=0; c<cols; ++c) {
    if (c==col) continue; 
    mcn = mat[c][row]; 
    if (mcn==pZERO) continue; 

    changed = true; 
    quot = gdivround(mat[c][row],smallest);
    mat[c] -= quot * mat[col]; 
  }

  return changed; 
}

static void swap_cols(pari& mat, int c1, int c2)
{
  if (c1==c2) return; 
  pari col1 = mat[c1]; 
  mat[c1] = mat[c2]; 
  mat[c2] = col1; 
}

// Do integer column operations to get
// 0,...,0,1 in the bottom row of the matrix 
// (if possible, return false if not).

static bool make_last_row_en(pari& mat)
{
  int c, cols = length(mat);
  int row = length(mat[0])-1;
  int guard = 1000; 
  
  while (--guard > 0) {
    if (!find_smallest_in_last_row(mat, c)) return false; 

    if (!reduce_other_cols(mat,c)) 
      break; // rest of last row is 0

    if (abs(mat[c][row])==pONE) 
      break; // what we wanted. 
  }

  if (guard <= 0) {
    warn("make_last_row_en in get_equations.cc failed!\n");
    return false; 
  }

  if (mat[c][row] < pZERO) 
    mat[c] = -mat[c]; // make it 1 (not -1).
  swap_cols(mat, c, cols-1); // put it in last column. 

  return mat[cols-1][row] == pONE;
}


// Must call choose_generators before this if with_parity==true. 

pari compute_flattening(Triangulation* manifold, bool with_parity, bool filled, bool with_signs)
{
  pari equations = gtrans(filled ? 
			  get_filled_equations(manifold, with_signs) :
			  get_complete_equations(manifold, with_signs)); 

  int n_shapes = get_num_tetrahedra(manifold); 
  pari flattening = cvector(2*n_shapes); 

  pari small = 1e-6; 

  int i, j, n; 
  if (with_signs) {
    // pari pls = concat(log_shapes(manifold),pI*pPI);
    pari pls = p_log_shapes(manifold);
    pari angle_sums = (equations * pls)/pPI; 

    // angle_sums.print(); 
    pari corr; 
    n = length(angle_sums);
    for (i=0; i<n; i++) {
      corr = ground(gimag(angle_sums[i])); 
      // if (!numerical_zero(angle_sums[i] - pI * corr)) {
      if (gabs(angle_sums[i] - pI * corr) > small) {
	printf("compute_flattening: equations not satisfied!\n");
	return flattening; 
      }
      equations[2*n_shapes][i] -= corr;
    }
  }

  pari ker = kerint(equations); 
  int kd = length(ker); 

  if (kd==0 || !make_last_row_en(ker)) {
    printf("Problem finding combinatorial flattening.\n");
    return flattening;
  }

  // No parity condition. 
  if (!with_parity) {
    for (i=0; i<2*n_shapes; i++) {
      flattening[i] = ker[kd-1][i];
    }
    return flattening; 
  }

  // cout << "Flattening matrix\n";
  // ker.print(); 

  pari parity = get_parity_matrix(manifold); 

  // cout << "Parity matrix\n";
  // parity.print(); 

  // We want flattening to be a linear combination of the 
  // columns of ker, including 1 x the last column, such that 
  // parity * flattening = 0 mod 2. 

  pari parity_ker((parity*ker)%pTWO);

  pari pker = ker_mod_p(parity_ker, pTWO);
  int pkd = length(pker);
  
  // Check for a solution
  if (pkd==0 || !make_last_row_en(pker)) {
    printf("Unable to satisfy flattening parity condition!\n");
    return flattening;
  }

  // Extract column of pker with 1 at the bottom
  pari psol = cvector(kd);
  for (i=0; i<kd; i++) psol[i] = pker[pkd-1][i];
  
  pari f = ker * psol; 
    
  flattening = cvector(2*n_shapes); 
  for (i=0; i<2*n_shapes; i++) flattening[i] = f[i]; 
  return flattening; 
}

/* $Id: get_equations.cc,v 1.2 2009/12/03 22:49:36 matthiasgoerner Exp $ */ 
