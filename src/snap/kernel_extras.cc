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
#include "kernel_extras.hh"
#include "snappea/kernel.h"
#include "snappea/unix_io.h"
#include "snappea/canonize.h"
#include "complex_lll.hh"
#include "tile.hh"
#include "printable.hh"
#include "warn.hh"
#include "snap_io.hh"
#include <cstdio>

using std::swap;
using std::min;
using std::cout;
using std::cerr;
using std::endl;

/*
 *  Here's a description of the standard set of manifolds,
 */

Triangulation *standard_set[MAX_CUSPS + 1][MAX_IN_LEVEL];

char *standard_name[MAX_CUSPS + 1][MAX_IN_LEVEL] = {
	    {"",      "",      "",      "",      "",      ""     },
	    {"",      "",      "",      "",      "",      ""     },
	    {"tr2.0", "tr2.1", "tr2.2", "tr2.3", "tr2.4", "tr2.5"},
	    {"tr3.0", "tr3.2", "tr3.3", "tr3.4", "tr3.5", "tr3.6"},
	    {"tr4.0", "tr4.1", "tr4.2", "tr4.4", "tr4.5", "tr4.6"},
	    {"tr5.0", "tr5.1", "tr5.2", "tr5.3", "tr5.4", ""     },
	    {"tr6.0", "",      "",      "",      "",      ""     }
	};

Boolean standard_set_present = FALSE;
 
int num_standards[MAX_CUSPS + 1] = {0, 0, 6, 5, 6, 5, 1};

#define NUM_A_LENGTHS 18

static double arith_lengths[] = {
  0.862554627662061032185692794,
  0.962423650119206894995517827,
  1.061275061905035652033018916,
  1.087070144995739099785272801,
  1.265948638401894900296366492,
  1.316957896924816708625046347,
  1.418316134968973230200327889,
  1.465715351947290521777344887,
  1.486022124876927092479168720,
  1.534394436502638886658057537,
  1.566799236972411078664056862,
  1.662885891058621075652485039,
  1.725109255324122064371385588,
  1.736596079922649323434519884,
  1.762747174039086050465218649,
  1.852266062700364849100248890,
  1.740021645304850859297651287,
  1.907925539233777287698103877
};



int num_filled_cusps(const Triangulation *manifold)
{
  int i, n, n_cusps;
  Cusp *cusp; 
  n_cusps = get_num_cusps((Triangulation *)manifold); 
  n = 0; 

  for(i=0; i<n_cusps; i++) {
    cusp = find_cusp((Triangulation *)manifold, i);
    if (!cusp->is_complete) n++; 
  }
  return n;
}

string get_filling_name(Triangulation *manifold, int style)
{
  int i, n_cusps;
  string name;
  Boolean complete;
  double m, l; 
  char scratch[40];

  name = get_triangulation_name(manifold); 
  n_cusps = get_num_cusps(manifold);

  for(i=0; i<n_cusps; i++) {

    get_cusp_info(manifold, i, 0, &complete, &m, &l, 0, 0, 0, 0, 0, 0);

    if (!complete) {
      if (fabs(m - floor(m + .5)) < 0.000001 && fabs(l - floor(l + .5)) < 0.000001) {
	if (style==0)
	  sprintf(scratch, "(%d,%d)", (int)floor(m+.5), (int)floor(l+.5));
	else 
	  sprintf(scratch, "_%d_%d", (int)floor(m+.5), (int)floor(l+.5));
      } else {
	if (style==0)
	  sprintf(scratch, "(%lf,%lf)", m, l);
	else 
	  sprintf(scratch, "_%lf_%lf", m, l);
      }
      name += scratch;
    }
    else {
      if (style==0) 
	name += "(,)";
      else if (style==1)
	name += "_0_0";
      else
	name += "_.";
    }
    if (style == 0 && i < n_cusps - 1) name += ' ';
  }

  return name;
}

list<Complex> all_vertices(Triangulation* manifold)
{
  list<Complex> vertex_list;
  Complex this_vertex;
  Tetrahedron *tet;
  int i;

  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next)
    {
      for (i = 0; i < 4; i++) {
	this_vertex = tet->corner[i];

	if(find(vertex_list.begin(), vertex_list.end(), this_vertex) == vertex_list.end())
	  vertex_list.push_back(this_vertex);
      }
    }
  return vertex_list;
}

static MoebiusTransformation word_to_Moebius(GroupPresentation *group, int *word)
{
  O31_matrix result_O31; // not used
  MoebiusTransformation result_Moebius;

  fg_word_to_matrix(group, word, result_O31, &result_Moebius); 
  return result_Moebius;
}

static void get_cusp_matrices(GroupPresentation *group, int index, MoebiusTransformation& meridian, MoebiusTransformation& longitude)
{
  int *word;

  word = fg_get_meridian(group,index);
  meridian = word_to_Moebius(group,word);
  fg_free_relation(word); /* my_free(word); */

  word = fg_get_longitude(group,index);
  longitude = word_to_Moebius(group,word);
  fg_free_relation(word); /* my_free(word);*/
}


line core_geodesic(Triangulation* manifold, GroupPresentation *group, int index)
{
  list<Complex> vertices = all_vertices(manifold); 
  bool found;
  MoebiusTransformation meridian, longitude;
  line mcore, lcore, core; 
  int nm, nl;

  get_cusp_matrices(group, index, meridian, longitude); 

  nm = fixed_points(meridian.matrix, mcore);
  nl = fixed_points(longitude.matrix, lcore);

  core = (nm > nl) ? mcore : lcore;

  if (nm < 1 && nl < 1) {
    warn("Unable to find any fixed points\n");
    return line(Zero, Zero); 
  }

  found = (find(vertices.begin(), vertices.end(), core.end[1]) != vertices.end()); 
  if (!found) {
    found = (find(vertices.begin(), vertices.end(), core.end[0]) != vertices.end()); 
    if (!found) {
      warn("Core geodesic found with neither end a vertex of the triangulation\n"); 
    }
    swap(core.end[0], core.end[1]);
  }
  return core;
}

void print_tet_corners(const Tetrahedron* tet)
{
  int i; 
  printf("["); 
  for (i=0; i<4; i++) {
    print(tet->corner[i]); 
    if (i<3) printf(","); 
  }
  printf("]\n"); 
}

void print_tet_corners(const MoebiusTransformation& mt, const Tetrahedron* tet)
{
  int i; 
  printf("["); 
  for (i=0; i<4; i++) {
    print(mt * tet->corner[i]); 
    if (i<3) printf(","); 
  }
  printf("]\n"); 
}

int generator_num(Tetrahedron* tet, int f)
{
  int gs = tet->generator_status[f]; 
  if (gs==inbound_generator || gs==outbound_generator) 
    return ((gs==outbound_generator) ? 
	    tet->generator_index[f]+1 : -tet->generator_index[f]-1); 
  return 0; 
}

/* Prints the tetrahedra, numbered from 0. For each
   tetrahedron it prints the fundamental group generators associated with 
   the faces of the tetrahedron, followed by the vertices. 
   A generator of 0 means that the face is not one of the boundary
   faces of the fundamental domain given by the union of the 
   tetrahedra. So for each face with 0 generator the triple of 
   vertices should appear again as a 0 generator face of some other
   tetrahedron. 
*/

void print_tetrahedra(Triangulation *manifold)
{
  int i, gn;
  Tetrahedron *tet;
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next)
    {
      printf("Tetrahedron %d: [", tet->index);
      for (i=0; i<4; i++) {
	gn = generator_num(tet,i); 
	printf("%c", gn ? tochar(gn) : '0');
	if (i<3) printf(",");
      }
      printf("]\n");

      print_tet_corners(tet);
    }
}

void print_tetrahedra(const MoebiusTransformation& mt, Triangulation *manifold)
{
  print(mt); 
  printf("\nImage tetrahedra:\n"); 
  Tetrahedron *tet;
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next)
    {
      print_tet_corners(mt,tet);
    }
}

void print_full_shapes(Triangulation *manifold)
{
  int i;
  Tetrahedron *tet;
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next)
    {
      for (i=0; i<3; i++) {
	fwprint(cout, tet->shape[filled]->cwl[ultimate][i].rect, 8); 
	if (i<2) printf(" "); 
      }
      printf("\n"); 
    }
}


void update_conjugacy(WEPolyhedron* p, O31_matrix const& c)
{
  p->conjugacy = inverse(c) * p->conjugacy; 
}

// For each face of the polyhedron prints the following, according to
// which bits of how are set:
// 0x1  Face group element (as an O31 matrix). 
// Ox2  Trace of face pairing element. 
// 0x4  Trace computed via word in the group.
// Ox8  Face group element from word (and conjugacy). 
// 0x10 Prefix, each entry with the face number (0 - n-1). 
// 0x20 Face group element as a word. 

void print_polyhedron(const WEPolyhedron* polyhedron, const GroupPresentation* g, unsigned int how, FILE* fp)
{
  WEFace    *face;
  int       i;
  MoebiusTransformation mt; 
  Complex tr; 
  O31_matrix mx;
  string wd; 

  for (face = polyhedron->face_list_begin.next,
       i = 0;
       face != &polyhedron->face_list_end;
       face = face->next,
       i++)
    {
      if (how & 0x80) 
	if (i%2) continue; // even faces only
      if (how & 0x10) 
	fprintf(fp, "Face %2d: ", i); 
      if (how & 0x1)
	face->group_element->print(fp); 
      if (how & 0x2) { // trace from face element 
	mt = MoebiusTransformation(*face->group_element);
	tr = mt.matrix[0][0] + mt.matrix[1][1];
	fprintf(fp, "%f+%f*i ", tr.real, tr.imag);
      }
      if (how & 0x4) { // trace from face word
	mt = word_to_Moebius(g, face->word);
	tr = mt.matrix[0][0] + mt.matrix[1][1];
	fprintf(fp, "%f+%f*i ", tr.real, tr.imag);
      }
      if (how & 0x8) { // face element from word
	word_to_matrix(g, face->word, mx, polyhedron->conjugacy, 1);
	mx.print(fp); 
      }
      if (how & 0x40) { // generating vertex for the face.
	O31_matrix const& ge(*face->group_element);
	fprintf(fp,"[% 10f,% 10f,% 10f,% 10f] ",
		ge(0,0),ge(1,0),ge(2,0),ge(3,0));
      }
      if (how & 0x20) {
	wd = string(face->word); 
	fprintf(fp, "%s", wd.c_str()); 
      }
      fprintf(fp,"\n"); 
    }
}

Boolean orientation_preserving_isometry(IsometryList *isometry_list, int& isometry_num)
{
  /* check orientations */ 
  Boolean or_pres, or_rev; 
  isometry_list_orientations(isometry_list,&or_pres, &or_rev);
  isometry_num = 0; 
  if (!or_pres) return FALSE; 

  /* find an isometry which is orientation preserving if possible */ 
  int j,img_cusp;
  int cusp_map[2][2];
  for (j=0; j<isometry_list->num_isometries; j++) {
    if (parity[isometry_list->isometry[j]->tet_map[0]] == 0) {
      isometry_num = j; 

      isometry_list_cusp_action(isometry_list,j,0,&img_cusp,cusp_map);
      if (cusp_map[0][0]>0) break; // choose [1,0;0,1] if it's available
    }
  }
  return TRUE; 
}


// perm will be one of 012, 021, 102, 120, 201, 210.
// pnum will then be   0,   1,   2,   3,   4,   5. 

int cs_contrib[6] = {  0,   -2,  0,   -1,  1,   -1 };

int face_mismatch_cs_term(Triangulation* manifold, int report)
{
  int f, f_opp, i, j, v_opp, sign, pnum;
  int perm[3]; 
  int sum = 0; 
  Tetrahedron *tet;
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {

    for (f=0; f<4; f++) {

      f_opp = EVALUATE(tet->gluing[f],f); 

      j = 0; 
      for (i=0; i<4; i++) {
	if (i == f) continue;
	v_opp = EVALUATE(tet->gluing[f],i);
	perm[j++] = (v_opp < f_opp) ? v_opp : v_opp-1;
      }

      pnum = 2 * perm[0] + int(perm[1] > perm[2]); 

      sign = (f & 0x1) ? -1 : 1; 

      if (report) {
	printf("Tet %d, Face %d -> Tet %d, Face %d, Perm [", tet->index, f, tet->neighbor[f]->index, f_opp);
	for (i=0; i<3; i++) {
	  printf("%d", perm[i]);
	  if (i<2) printf(",");
	}
	printf("], Contrib %d\n", sign * cs_contrib[pnum]); 
      }

      sum += sign * cs_contrib[pnum];
    }
  }
  return sum; 
}

void load_standard_set(const string& path, int report)
{
  if (standard_set_present) return; 
  
  int i, j;

  Boolean is_known, requires_init;
  int precision; 
  double chern_simons; 

  if (report) printf("Loading standard set...\n"); 

  for (i = 2; i <= MAX_CUSPS; i++) {
    for (j = 0; j < num_standards[i]; j++) {

      standard_set[i][j] = read_manifold_file(path, standard_name[i][j]);
      if (!standard_set[i][j]) continue;

      get_CS_value(standard_set[i][j], &is_known, &chern_simons, &precision, &requires_init);

      if (report > 1) 
	printf("%s %9.6lf %9.6lf\n", standard_name[i][j],
	       volume(standard_set[i][j], NULL), chern_simons);
    }
  }

  standard_set_present = 1;
}

Boolean belongs_to_standard_set(Triangulation *manifold, 
				IsometryList **isometry_list, 
				int *which)
{
    int             j;
    Boolean         isometric;
    int n_cusps = get_num_cusps(manifold); 

    if (n_cusps < 2 || n_cusps > MAX_CUSPS)
	return(FALSE);

    for (j = 0; j < num_standards[n_cusps]; j++) {

      if (fabs(volume(manifold, NULL) - 
	       volume(standard_set[n_cusps][j], NULL)) > 0.01)
	continue;
      
      compute_isometries(manifold, 
			 standard_set[n_cusps][j], 
			 &isometric, 
			 isometry_list, NULL);

      if (!isometric) {
	if (*isometry_list) free_isometry_list(*isometry_list); 
	*isometry_list = 0; 
	continue; 
      }
      *which = j; 
      return TRUE;
    }
    return FALSE; 
}      

Boolean drill_to_standard(Triangulation *manifold, 
			  Triangulation **drilled, 
			  IsometryList **isometry_list, 
			  Triangulation **standard,
			  int& which_standard,
			  int report)
{
    int                     num_curves,
			    curves_to_try,
			    i;
    DualOneSkeletonCurve    **the_curves;
    Triangulation           *new_triangulation, *drilled_manifold; 

    /* can't handle non-orientable manifolds */ 
    if (get_orientability(manifold) != oriented_manifold) {
      warn("manifold is not orientable\n"); 
      return FALSE; 
    }

    /* make a copy, save the original */ 
    copy_triangulation(manifold, &drilled_manifold);
    
    /* already in standard set ? */
    if (belongs_to_standard_set(drilled_manifold, isometry_list, &which_standard)) {
      *drilled = drilled_manifold; 
      *standard = standard_set[get_num_cusps(*drilled)][which_standard];
      set_triangulation_name(*standard,
			     standard_name[get_num_cusps(*drilled)][which_standard]);
      return TRUE; 
    }


    int cusp_count;
    int dn_cusps = get_num_cusps(drilled_manifold); 
    for (cusp_count = 2; cusp_count <= MAX_CUSPS; cusp_count++) {

      if (report) {
	// print_manifold_info(drilled_manifold);
	printf("not in standard set\n"); 
	printf("finding dual curves\n"); 
      }

      if (dn_cusps < cusp_count) {
	dual_curves(drilled_manifold, 12, &num_curves, &the_curves);

	if (num_curves == 0) {
	  warn("got stuck trying to find dual curves\n");
	  return FALSE; 
	}

	curves_to_try = min(num_curves, 3);

	Complex complete_length, filled_length; 
	MatrixParity parity; 
	    
	for (i = 0; i < curves_to_try; i++) {
    
	  if (report) {
	    get_dual_curve_info(the_curves[i], &complete_length, &filled_length, &parity);
	    printf("drilling curve: "); 
	    print(filled_length); 
	    printf("\n"); 
	  }

	  new_triangulation = drill_cusp(drilled_manifold, the_curves[i], 
					 get_triangulation_name(drilled_manifold));
	  if (new_triangulation == NULL) {
	    warn("failed to drill manifold\n");
	    free_dual_curves(num_curves, the_curves);
	    free_triangulation(drilled_manifold);
	    return FALSE; 
	  }

	  find_complete_hyperbolic_structure(new_triangulation);
	  if (get_complete_solution_type(new_triangulation) != geometric_solution) {
	    warn("bad solution\n");
	    free_triangulation(new_triangulation);
	    continue;
	  }

	  if (belongs_to_standard_set(new_triangulation, isometry_list, &which_standard)) {
	    free_dual_curves(num_curves, the_curves);
	    free_triangulation(drilled_manifold);
	    *drilled = new_triangulation; 
	    *standard = standard_set[get_num_cusps(*drilled)][which_standard];
	    char* name = standard_name[get_num_cusps(*drilled)][which_standard];
	    set_triangulation_name(*standard, name);
	    if (report) printf("drilled = %s\n", name); 
	    return TRUE;
	  }

	  /* That curve didn't get us to a standard one so go 
	     back and try the next one in the list. */ 
	  free_triangulation(new_triangulation);
	}

	/* None of the above curves got us to a standard manifold
	   so we'll go back to drilling out the shortest and see if 
	   the next curve we drill out works. */ 
	new_triangulation = drill_cusp(drilled_manifold, the_curves[0], 
				       get_triangulation_name(drilled_manifold));
	find_complete_hyperbolic_structure(new_triangulation);


	if (report) {
	  get_dual_curve_info(the_curves[0], &complete_length, &filled_length, &parity);
	  printf("returning to drilling shortest curve: "); 
	  print(filled_length); 
	  printf("\n"); 
	}


	free_dual_curves(num_curves, the_curves);
	free_triangulation(drilled_manifold);
	drilled_manifold = new_triangulation;
      }
    }

    free_triangulation(drilled_manifold);
    return FALSE;
}

/* Returns the column position of the standard manifold found. */

int get_filled_standard(Triangulation* manifold, Triangulation** result, int report)
{
  if (!standard_set_present) {
    warn("reference manifolds have not been loaded\n"); 
    return -1; 
  }

  // Drill curves until we get to a manifold in the standard set
  int which_standard; 
  Triangulation *drilled, *standard_tri;
  IsometryList *isom_list = 0;
  if (!drill_to_standard(manifold, &drilled, &isom_list, &standard_tri, which_standard, report)) return -1;

  /* Copy the standard manifold */ 
  copy_triangulation(standard_tri, result); 

  /* Look for an orientation preserving isometry */ 
  int isometry_num; 
  orientation_preserving_isometry(isom_list, isometry_num); 

  /* Refill the drilled cusps via this isometry */ 
  int j, image_cusp; 
  int n_cusps = get_num_cusps(manifold); 

  int dn_cusps = get_num_cusps(drilled); 
  MatrixInt22 cusp_map; 

  for (j=0; j<dn_cusps; j++) {
    
    /* Find the index of the corresponding standard manifold cusp. */ 
    isometry_list_cusp_action(isom_list, isometry_num, j, &image_cusp, cusp_map); 

    /* Set the filling coeffs of the result manifold: the first 
       n_cusps of the drilled manifold correspond to 
       cusps of the original manifold. Therefore they remain complete.
       The others we want to fill by the equivalent of (1,0) dehn 
       filling on the drilled manifold. */ 

    if (j < n_cusps) {
      set_cusp_info(*result, image_cusp, TRUE, 0, 0); 
    } else {
      set_cusp_info(*result, image_cusp, FALSE, 
		  (double)cusp_map[0][0], (double)cusp_map[1][0]); 
    }
  }

  do_Dehn_filling(*result); 

  /* a quick check that everything is OK */ 
  if (1e-10 < fabs(volume(manifold, NULL) - volume(*result, NULL))) {
    warn("filled drilled volume differs from original volume\n"); 
    return -1; 
  }

  return which_standard; 
}

Complex cusp_vertex(Triangulation* manifold, int cusp_num)
{
  // This works for 1-cusped manifolds only. 
  Cusp* cusp = manifold->cusp_list_begin.next;
  while (cusp->next != &manifold->cusp_list_end && cusp_num > 0) {
    cusp = cusp->next; 
    cusp_num--;
  }
  if (cusp_num != 0) {
    cout << "Invalid cusp number in cusp_vertex function.\n";
    return Zero; 
  }
  return cusp->basepoint_tet->corner[cusp->basepoint_vertex];
}

bool cusp_is_complete(Triangulation* manifold, int i)
{ 
  return find_cusp(manifold,i)->is_complete; 
}

// call with length of shortest geodesic
// true means manifold is non-arithmetic, 
// false means we don't know whether the manifold is arith or not. 

#define L_EPS 1e-7

bool length_is_non_arithmetic(double l)
{
  if (l > arith_lengths[1] - L_EPS) return false; 
  int i=0; 
  l *= 2; 
  while (i<NUM_A_LENGTHS) {
    if (l < arith_lengths[i] + L_EPS) break; 
    i++;
  }
  if (i==NUM_A_LENGTHS) return true; // l > biggest in table. 
  if (l > arith_lengths[i] - L_EPS) return false; // l occurs in table. 

  return true; // l not in table. 
}

ShapeInversion **shape_histories(Triangulation *manifold)
{
  int n = get_num_tetrahedra(manifold);
  ShapeInversion **sh = new ShapeInversion* [n];

  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    sh[t->index] = t->shape_history[filled]; 

  return sh;
}

vector<Complex> get_log_shapes(Triangulation* manifold)
{
  int i = 0, n = get_num_tetrahedra(manifold); 
  vector<Complex> shapes(2*n);

  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    {
      shapes[i]   =  t->shape[filled]->cwl[ultimate][0].log;
      shapes[i+n] = -t->shape[filled]->cwl[ultimate][1].log;
      i++; 
    }
  return shapes;
}


void get_edge_equations(Triangulation *manifold, int_matrix& equations)
{
  Tetrahedron *tet;
  EdgeIndex   e;
  int i, r; 
  
  int n_tet = get_num_tetrahedra(manifold); 
  int n_edges = n_tet; 
  int last_col = 2*n_tet;

  equations.set_dims(n_edges, 2*n_tet + 1); 
  equations = 0; 

  /*
    Equations expresses the gluing equations in the following form:
    Each row contains 2*n+1 entries, [a1,...,an,b1,...,bn,c] 
    (where n = n_tet) representing the equation 
    z1^a1...zn^an.(1-z1)^b1.(1-zn)^bn = e^(pi i c).
  */

  /* can't do anything if it's not an oriented manifold */ 
  if (manifold->orientability != oriented_manifold) return; 
  
  /*  number tetrahedra and edge classes */
  number_the_tetrahedra(manifold);
  number_the_edge_classes(manifold);
  
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {
    
    i = tet->index; 
    
    /*
     *  Edge equations.
     */
    for (e = 0; e < 6; e++) { /* Look at each of the six edges. */
      
      r = tet->edge_class[e]->index;
      
      /*
       * edge3[e] = 0, 1, 2 according to whether the contribution
       * to the equation for this edge is z, 1/(1-z), (z-1)/z. 
       * (where z is the parameter for this tetrahedron). 
       */
      switch (edge3[e]) {
      case 0:
	equations[r][i] += 1;
	break;
      case 1:
	equations[r][i+n_tet] -= 1;
	break;
      case 2:
	equations[r][i] -= 1;
	equations[r][i+n_tet] += 1;
	equations[r][last_col] += 1;
	break;
      default:
	break;
      }
    }
  } /* for (tet...) */ 
  
    /* angle sum around edges is 2 PI */ 
  for (r = 0; r < n_edges; r++)
    equations[r][last_col] -= 2; 
}

void get_cusp_equations(Triangulation *manifold, int_matrix& equations)
{
  Tetrahedron *tet;
  VertexIndex v;
  FaceIndex initial_side, terminal_side;
  int init[2], term[2];
  int i, r, contrib, ml; 
  
  int n_tet = get_num_tetrahedra(manifold); 
  int n_cusps = manifold->num_cusps; 
  int last_col = 2*n_tet;
  
  equations.set_dims(2*n_cusps, 2*n_tet + 1); 
  equations = 0;
  
  /* 
     Equations expresses the gluing equations in the following form:
     First we have rows for meridians, then longitudes.
     Each row contains 2*n+1 entries, [a1,...,an,b1,...,bn,c] 
     (where n = n_tet) representing the equation 
     z1^a1...zn^an.(1-z1)^b1.(1-zn)^bn = e^(pi i c).
  */
  
  /* can't do anything if it's not an oriented manifold */ 
  if (manifold->orientability != oriented_manifold) return; 
  
  /*  number tetrahedra and edge classes */
  number_the_tetrahedra(manifold);
  
  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {
    i = tet->index; 

    /*
     *    Cusp equations.
     */
    for (v = 0; v < 4; v++) { /* Look at each ideal vertex. */
      
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
	
	for (ml = 0; ml < 2; ml++)    {       /* meridian,longitude */
	  /*
	   *  Note the intersection numbers of the meridian and
	   *  longitude with the initial and terminal sides (on 
	   *  the right-handed triangular cross section).
	   */
	  init[ml] = tet->curve[ml][0][v][initial_side];
	  term[ml] = tet->curve[ml][0][v][terminal_side];
	  
	  contrib = (int)FLOW(init[ml],term[ml]);
	  
	  r = tet->cusp[v]->index + ml*n_cusps;       
	  
	  switch (edge3_between_faces[initial_side][terminal_side]) {
	  case 0:
	    equations[r][i] += contrib; 
	    break;
	  case 1:
	    equations[r][i+n_tet] -= contrib;
	    break;
	  case 2:
	    equations[r][i] -= contrib;
	    equations[r][i+n_tet] += contrib;
	    equations[r][last_col] += contrib;
	    break;
	  default:
	    break;
	  }
	}
	
      } /* for (initial_side...) */ 
    } /* for (v...) */ 
  }
}

static const double r_eps=1e-9; 

static void ratapp(double x, int& n, int& d)
{
  double der=1.;
  int xint, t, i, c=0; 
  vector<int> cf;

  while (c < 30) {
    xint = int(floor(x+r_eps)); 
    cf.push_back(xint); 
    x -= xint; 
    x = 1.0/x;
    der /= (x*x); 
    if (der < r_eps) break; 
    c++;
  } 

  n=0;
  d=1;
  for (i=cf.size()-1; i>=0; i--) {
    n += cf[i] * d; 
    if (i>0) {
      t = n; n = d; d = t; // swap n & d. 
    }
  }
}

void ratsc(double m, double l, int& mi, int& li, int& d)
{
  int d1, d2, g;

  ratapp(m, mi, d1);
  ratapp(l, li, d2);

  g = gcd(d1,d2); 

  mi *= (d2/g);
  li *= (d1/g); 
  d = d1*(d2/g);
}

void get_surgery_coeffs(Triangulation* manifold, vector<int>& coeffs)
{
  Cusp *cusp; 

  int n_cusps = manifold->num_cusps; 
  coeffs.resize(3*n_cusps); 
  int m, l, d; 

  for (cusp = manifold->cusp_list_begin.next;
       cusp != &manifold->cusp_list_end;
       cusp = cusp->next) {
    if (cusp->is_complete) {
      // For complete cusp, (m, l, denominator) = (1,0,0).
      coeffs[cusp->index] = 1; 
      coeffs[cusp->index + n_cusps] = 0; 
      coeffs[cusp->index + 2*n_cusps] = 0; 
    } else {
      ratsc(cusp->m, cusp->l, m, l, d); 
      coeffs[cusp->index] = m; 
      coeffs[cusp->index + n_cusps] = l;
      coeffs[cusp->index + 2*n_cusps] = d; 
    }
  }
}

void get_filling_equations(Triangulation *manifold, int_matrix& eqns)
{
  int i, j, n = manifold->num_cusps, n_tet = get_num_tetrahedra(manifold);

  eqns.set_dims(n, 2*n_tet + 1);
  eqns = 0; 

  vector<int> coeffs(3*n); 
  get_surgery_coeffs(manifold, coeffs); 
  int_matrix c_eqs; 
  get_cusp_equations(manifold, c_eqs); 

  // meridians then longitudes in complete_eqs. 

  int m, l, d; 
  for (i=0; i<n; i++) {
    m = coeffs[i]; l = coeffs[i+n]; d = coeffs[i+2*n]; 

    for (j=0; j <= 2*n_tet; j++) {
      eqns[i][j] = m * c_eqs[i][j] + l * c_eqs[i+n][j];
    }
    eqns[i][2*n_tet] -= 2 * d; 
  }
}

void get_full_edge_equations(Triangulation *manifold, int_matrix& eqns)
{
  Tetrahedron   *tet;
  EdgeIndex e;
  int i, r; 

  int n_tet = get_num_tetrahedra(manifold); 
  int n_edges = n_tet; 

  eqns.set_dims(n_edges, 3*n_tet); 
  eqns = 0; 

  /*
    Equations expresses the gluing equations in the following form:
    Each row contains 3*n entries, [a1,b1,c1,...,an,bn,cn] (where n = n_tet)
    representing the edge product z1^a1.(1/(1-z1))^b1.(1-1/z1)^c1....
  */

  /* Can't do anything if it's not an oriented manifold */ 
  if (manifold->orientability != oriented_manifold) return; 
    
  /* Number tetrahedra and edge classes. */
  number_the_tetrahedra(manifold);
  number_the_edge_classes(manifold);

  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {

    i = tet->index; 

    /*
     *  Edge equations.
     */
    for (e = 0; e < 6; e++) {   /* Look at each of the six edges. */

      r = tet->edge_class[e]->index;

      /*
       * edge3[e] = 0, 1, 2 according to whether the contribution
       * to the equation for this edge is z, 1/(1-z), (z-1)/z. 
       * (where z is the parameter for this tetrahedron). 
       */
      eqns[r][3*i + edge3[e]] += 1; 
    }
  } /* for (tet...) */ 
}

void get_full_cusp_equations(Triangulation *manifold, int_matrix& eqns)
{
  Tetrahedron *tet;
  VertexIndex v;
  FaceIndex initial_side, terminal_side;
  int init[2], term[2];
  int i, r, contrib, ml; 
  int etype; 

  int n_tet = get_num_tetrahedra(manifold); 
  int n_cusps = manifold->num_cusps; 

  eqns.set_dims(2*n_cusps, 3*n_tet); 
  eqns = 0; 

  /* 
     Equations expresses the gluing equations in the following form:
     Each row contains 3*n entries, [a1,b1,c1,...,an,bn,cn] (where n = n_tet)
     representing the meridian or longitude product 
     z1^a1.(1/(1-z1))^b1.(1-1/z1)^c1....
  */

  /* Can't do anything if it's not an oriented manifold */ 
  if (manifold->orientability != oriented_manifold) return; 
    
  /* Number tetrahedra and edge classes */
  number_the_tetrahedra(manifold);

  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {

    i = tet->index; 

    /*
     *  Cusp equations.
     */
    for (v = 0; v < 4; v++) {   /* Look at each ideal vertex. */
	
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

	  r = 2 * tet->cusp[v]->index + ml;     

	  etype = edge3_between_faces[initial_side][terminal_side];

	  eqns[r][3*i + etype] += contrib; 
	}

      } /* for (initial_side...) */ 
    } /* for (v...) */ 
  }
}

void check_no_generator_transparent(Triangulation* manifold)
{
  Tetrahedron *tet, *nbr_tet;
  FaceIndex   face, nbr_face;
  Permutation gluing; 
  double sum_of_tilts;

  if (!manifold->proto_canonized) {
    printf("Proto_canonized flag not set for this manifold.\n"); 
    return;
  }

  for (tet = manifold->tet_list_begin.next;
       tet != &manifold->tet_list_end;
       tet = tet->next) {
    
    for (face = 0; face < 4; face++) {
      if (tet->generator_status[face] == not_a_generator) continue; 
      
      nbr_tet     = tet->neighbor[face];
      gluing      = tet->gluing[face];
      nbr_face    = EVALUATE(gluing, face);
      
      sum_of_tilts = tet->tilt[face] + nbr_tet->tilt[nbr_face];
      if (sum_of_tilts >= -CONCAVITY_EPSILON) {
	printf("Transparent face has generator!\n");
	return;
      }
    }
  }
}

static void add_vertex(vector<Complex>& vc, Complex const& z)
{
  int i, n=vc.size();
  for (i=0; i<n; i++) {
    if (same_point(vc[i],z)) return;
  }
  vc.push_back(z);
}

void print_ideal_cells(Triangulation* manifold)
{
  Cusp *cusp; 
  Tetrahedron* tet; 
  VertexIndex v; 
  vector<Complex> vl; 

  for (cusp = manifold->cusp_list_begin.next;
       cusp != &manifold->cusp_list_end;
       cusp = cusp->next) {
    if (!cusp->is_finite) continue; 

    vl.resize(0); 

    /* Find all vertices for this "cusp" */ 
    for (tet = manifold->tet_list_begin.next;
	 tet != &manifold->tet_list_end;
	 tet = tet->next) {
      for (v = 0; v < 2; v++) {
	if (tet->cusp[v] != cusp) continue; 

	add_vertex(vl, tet->corner[v]);
	add_vertex(vl, tet->corner[v+2]);
      }
    }

    cout << PSeq(vl) << endl; 
  }
}


// Here follows the code required for 
// seek_acyclic_edge_orientations. 


// this is the orientation of a given edge of this
// tetrahedron, relative to a fixed orientation of the
// edges (EdgeClasses) given by set_edge_orientation_flags. 

static int base_orientation(Tetrahedron* tet, int v1, int v2)
{
  EdgeIndex e = edge_between_vertices[v1][v2];
  return (!edge_direction(tet,e) == (v1 < v2)) ? 1 : -1; 
}

Permutation tet_perm(Triangulation* m, int index, vector<int> const& e_or)
{
  Tetrahedron* tet = find_tetrahedron(m, index); 
  if (!tet) return 0; 
  return tet_perm(tet, e_or); 
}

Permutation tet_perm(Tetrahedron* tet, vector<int> const& e_or)
{
  int i, j, c, e, p[4]; 

  for (i=0; i<4; i++) {
    c = 0; 
    for (j=0; j<4; j++) {
      if (j==i) continue; 
      e = tet->edge_class[edge_between_vertices[i][j]]->index;
      if (base_orientation(tet, i, j)!=e_or[e]) ++c;
    }
    p[i] = c; 
  }

  return CREATE_PERMUTATION(0,p[0],1,p[1],2,p[2],3,p[3]);
}


// returns -1 if edges are cyclic around the face
//          0 if not all edges are oriented and none are forced
//          1 if all edges are oriented and acyclic
//          2 if an edge was forced. 

static int check_face(Tetrahedron* tet, FaceIndex f, vector<int>& e_or)
{

  // get the three vertices of the face
  int v[3], v1, v2; 
  int i, j=0; 
  for (i=0; i<4; i++) {
    if (i==f) continue; 
    v[j++] = i; 
  }

  // add up the edge orientations as we go cyclically around face. 
  int orsum = 0;
  int e, b_or, unor=-1, unor_b;
  v1 = v[2]; 
  for (i=0; i<3; i++) {
    v2 = v[i]; 
    e = tet->edge_class[edge_between_vertices[v1][v2]]->index; 
    b_or = base_orientation(tet, v1, v2);
    orsum += b_or * e_or[e];
    if (!e_or[e]) { unor = e; unor_b = b_or; } // save base orientation of an unoriented edge
    v1 = v2; 
  }

  // the results are.. 
  if (orsum==3 || orsum==-3) return -1; // cyclic
  if (orsum==2 || orsum==-2) { // force the unoriented one. 
    if (unor < 0)
      uFatalError("check_face", "kernel_extras"); // should never happen
    e_or[unor] = ((orsum==2) ? -1 : 1) * unor_b; // reduce sum to 1 or -1; 
    return 2; // an edge orientation was forced.
  } 
  if (unor==-1) return 1; // all oriented and acyclic (orsum should be +/- 1)
  return 0; // some unoriented edge, not forced. 
}


static bool orient_edge(Triangulation* m, vector<int>& e_or, vector<int>& tet_done, int e, int ori)
{

  // set the orientation of the given edge. 
  e_or[e] = ori; 

  // check consistency on the tetrahedra

  bool recheck = true, done; 

  Tetrahedron* tet; 
  FaceIndex f; 
  while (recheck) {

    recheck = false; 

    // check consistency on each tetrahedron

    for (tet = m->tet_list_begin.next;
	 tet != &m->tet_list_end;
	 tet = tet->next) {

      if (tet_done[tet->index]) continue; 

      done=true; // assume oriented and acyclic until proven otherwise
      for (f=0; f<4; f++) {

	switch (check_face(tet, f, e_or)) {

	case -1: // cyclic! 
	  return false; 
	  break;
	case 0: // ok but incomplete
	  done=false; 
	  break; 
	case 1: // ok and complete
	  break; 
	case 2: // an edge was forced
	  recheck=true; 
	  done=false; // tet could now have a cycle
	  break; 
	default:
	  break;
	}
      }

      // if we get here things are still consistent

      if (done) // every face check returned ok and complete. 
	tet_done[tet->index]=1; 
    }
  }

  return true; 
}

// this is called with a partial consitent orientation of the edges.

static bool orient_the_edges(Triangulation* m, vector<int>& e_or, vector<int>& tet_done, bool print_all)
{
  int e, n=get_num_tetrahedra(m); 

  for (e=0; e<n; e++) // find an unoriented edge. 
    if (e_or[e]==0) break; 

  if (e==n) {
    if (!print_all) return true;
    cout << PSeq(e_or) << endl; 
    return false; 
  }

  vector<int> current_e_or = e_or; 
  vector<int> current_tet_done = tet_done; 

  // try orienting edge e +1 (followed by recursively doing the rest of them). 

  if (orient_edge(m, e_or, tet_done, e,  1) &&
      orient_the_edges(m, e_or, tet_done, print_all)) return true; 

  // go back to where we were. 

  e_or = current_e_or; 
  tet_done = current_tet_done; 

  // try orienting edge e -1 (and then the rest of them). 

  if (orient_edge(m, e_or, tet_done, e, -1) &&
      orient_the_edges(m, e_or, tet_done, print_all)) return true; 

  return false; // no consistent orientation starting with this partial. 
}

#if 0
static void clear_flags(Triangulation* m)
{
  Tetrahedron* tet; 
  for (tet = m->tet_list_begin.next;
       tet != &m->tet_list_end;
       tet = tet->next) {

    tet->flag = 0; 

  }
}
#endif

bool seek_acyclic_edge_orientations(Triangulation* m, vector<int>& e_or, bool print_all)
{
  int i, n = get_num_tetrahedra(m); 
  if (e_or.size()!=n) e_or.resize(n); 
  for (i=0; i<n; i++) e_or[i] = 0; // unoriented
  vector<int> tet_done(n,0);

  number_the_edge_classes(m); 

  return orient_the_edges(m, e_or, tet_done, print_all);
}

bool my_canonize(Triangulation** m, n_vector const& cs, bool& copy)
{
  Triangulation* old = *m; 
  copy_triangulation(*m, m);
  copy = true; 

  if (proto_canonize(*m, cs.V)== func_failed) return false; 
  choose_generators(*m, TRUE, FALSE); // Compute the corners. 
  
  int j; 
  bool notcan = !is_canonical_triangulation(*m);
  if (notcan)
    canonical_retriangulation(*m); // This will have corners as well. 

  for (j=0; j<2; j++) 
    (*m)->solution_type[j] = old->solution_type[j]; 
  return true; 
}

bool check_if_isometric(Triangulation* a, Triangulation* b, int report)
{
  /* look for isometries between the two */ 
  IsometryList *isometry_list = 0;
  Boolean isometric; 
  compute_isometries(a,b,&isometric,&isometry_list,NULL); 
  if (!report) {
    if (isometry_list) free_isometry_list(isometry_list); 
    return isometric; 
  }

  if (isometric) cout << "Manifolds are isometric.\n";
  else cout << "Manifolds are not isometric.\n";

  if (!isometric || !isometry_list) return false; 

  // PRINT THE CUSP ACTION

  // In an isometry list the complete cusps are numbered from zero, 
  // while filled cusps are discounted entirely. Therefore we make a 
  // (couple of) tables so we can go from an IsometryList cusp index
  // to a Triangulation cusp index. 

  int i, j = 0, c = get_num_cusps(a); 
  vector<int> a_cusp_num(c); 
  for (i=0; i<c; i++) {
    if (cusp_is_complete(a,i)) a_cusp_num[j++] = i; 
  }

  j = 0; c = get_num_cusps(b); 
  vector<int> b_cusp_num(c);
  for (i=0; i<c; i++) {
    if (cusp_is_complete(b,i)) b_cusp_num[j++] = i; 
  }
  int n_cusps = j; 

  // The rest is straightforward. 

  // Choose an isometry which is orientation preserving, if possible. 
  int isom_num;
  Boolean or_pres = orientation_preserving_isometry(isometry_list, isom_num);

  if (!or_pres) cout << "Orientation reversing isometry\n";

  int mx[2][2]; 
  int image_cusp; 
  for (j=0; j < n_cusps; j++) {
    isometry_list_cusp_action(isometry_list, isom_num, j, &image_cusp, mx); 

    cout << a_cusp_num[j] << " -> " << b_cusp_num[j] << ' '; 
    cout << '[' << mx[0][0] << ',' << mx[0][1] << ';' << 
      mx[1][0] << ',' << mx[1][1] << ']' << endl;
  }

  if (isometry_list) free_isometry_list(isometry_list); 
  return true; 
}

