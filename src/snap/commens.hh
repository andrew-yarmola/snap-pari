#ifndef COMMENS_H
#define COMMENS_H
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

#include <vector>
#include "n_vector.hh"
#include "snappea/SnapPea.h"

using std::vector;

typedef struct graphnode GraphNode;

int commensurabilities(Triangulation* manifold0, Triangulation* manifold1, 
		       vector<MoebiusTransformation>& TR, bool& or_pres,
		       GraphNode** mapping_graphs = NULL);

vector<int> cusp_hidden_symmetry_classes(int num_cusps, GraphNode* head_graph);

void free_mapping_graphs(GraphNode *);

vector<int> cusp_covering_degrees(Triangulation *m0, Triangulation *m1, GraphNode* head_graph);

vector<int> get_cusp_size_vector(Triangulation* manifold, GraphNode* head_graph, n_vector& sizes);

bool has_finite_vertices(Triangulation *manifold);

void set_same_class(int a, int b, vector<int>& classes);

int normalize_independent_cusps(n_vector& areas, vector<int> const& classes);

int scan_ratio_vector(const char* str, n_vector& v);
void sprint_vector(char* buf, vector<int> const& v, int w=0);
void sprint_ratio_vector(char* buf, n_vector const& v);

void set_to_integer_multiple(n_vector& v);
void remove_common_factor(vector<int>& v);

#endif
