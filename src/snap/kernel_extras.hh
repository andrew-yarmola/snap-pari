#ifndef _kernel_extras_
#define _kernel_extras_
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
#include "lines.hh"
#include "int_matrix.hh"
#include "n_vector.hh"
#include <list>
#include <vector>
#include <algorithm>

using std::list;
using std::vector;
using std::string;

// Declarations for the standard set of manifolds used by drill_to_standard. 
#define MAX_CUSPS 6
#define MAX_IN_LEVEL 6
extern Triangulation *standard_set[MAX_CUSPS + 1][MAX_IN_LEVEL];
extern char *standard_name[MAX_CUSPS + 1][MAX_IN_LEVEL];
extern int num_standards[MAX_CUSPS + 1];
extern Boolean standard_set_present;
void load_standard_set(const string& search_path, int report = 0);
Boolean drill_to_standard(Triangulation *manifold, Triangulation **drilled, IsometryList **isometry_list, Triangulation **standard, int& which_standard, int report);
int get_filled_standard(Triangulation* manifold, Triangulation** result, int report);


int num_filled_cusps(const Triangulation *manifold);
bool cusp_is_complete(Triangulation* manifold, int i);
string get_filling_name(Triangulation *manifold, int style = 0);

list<Complex> all_vertices(Triangulation* manifold);

line core_geodesic(Triangulation* manifold, GroupPresentation *group, int index);

void print_tetrahedra(Triangulation *manifold);
void print_tetrahedra(const MoebiusTransformation& mt, Triangulation *manifold);
void print_full_shapes(Triangulation *manifold);

void print_polyhedron(const WEPolyhedron* polyhedron, const GroupPresentation* g, unsigned int how, FILE* fp = stdout);
void update_conjugacy(WEPolyhedron* p, O31_matrix const& c);

bool check_if_isometric(Triangulation* a, Triangulation* b, int report);
Boolean orientation_preserving_isometry(IsometryList *isometry_list, int& isometry_num);


int face_mismatch_cs_term(Triangulation* manifold, int report);

bool my_canonize(Triangulation** m, n_vector const& cs, bool& copy);


// Defined in peripheral_curves.c. Not declared in SnapPea because that would
// involve including vector.h in SnapPea.h. 
Boolean get_end_pairing_words(Triangulation* manifold, vector<FGWord>& words);

Complex cusp_vertex(Triangulation* manifold, int cusp=0);

bool length_is_non_arithmetic(double l); 

struct ShapeInversion; 
ShapeInversion** shape_histories(Triangulation *manifold);
vector<Complex> get_log_shapes(Triangulation* manifold);

void ratsc(double m, double l, int& mi, int& li, int& d);

void get_surgery_coeffs(Triangulation* manifold, vector<int>& coeffs);
void get_edge_equations(Triangulation *manifold, int_matrix& equations);
void get_cusp_equations(Triangulation *manifold, int_matrix& equations);
void get_filling_equations(Triangulation *manifold, int_matrix& eqns);
void get_full_cusp_equations(Triangulation *manifold, int_matrix& eqns);
void get_full_edge_equations(Triangulation *manifold, int_matrix& eqns);

void check_no_generator_transparent(Triangulation* manifold);
void print_ideal_cells(Triangulation* manifold);

bool seek_acyclic_edge_orientations(Triangulation* m, vector<int>& e_or, bool print_all=false);

struct Tetrahedron; 
unsigned char tet_perm(Triangulation* m, int index, vector<int> const& e_or);
unsigned char tet_perm(Tetrahedron* tet, vector<int> const& e_or); 

#endif
