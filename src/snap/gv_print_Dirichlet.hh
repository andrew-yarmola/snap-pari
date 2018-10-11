#ifndef _gv_print_Dirichlet_
#define _gv_print_Dirichlet_
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
#include <vector>

using std::vector;

void gv_print(FILE* fp, WEPolyhedron* polyhedron);
void gv_print_w_axes(FILE* fp, WEPolyhedron* poly);

void gv_print_w_geodesics(FILE* fp, WEPolyhedron* polyhedron, vector<GeodesicWord> const& geodesics, vector<int> const& show);
void gv_print_geodesic_w_ortholines(FILE* fp, WEPolyhedron* poly, GeodesicWord const& gw, vector<Ortholine> const& ort);

#endif
