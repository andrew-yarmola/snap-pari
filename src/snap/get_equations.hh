#ifndef _get_equations_
#define _get_equations_
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
#include "to_pari.hh"

pari cusp_shapes(const Triangulation* manifold);
pari tet_shapes(Triangulation *manifold); 
pari log_shapes(Triangulation *manifold); 
pari surgery_coeffs(Triangulation* manifold);
pari surgery_coeffs2(Triangulation* manifold);
pari edge_equations(Triangulation *manifold, bool with_signs=false, bool full=false);
pari cusp_equations(Triangulation *manifold, bool full=false);
pari get_complete_equations(Triangulation *manifold, bool with_signs=false, bool full=false);
pari get_filled_equations(Triangulation *manifold, bool with_signs=false, bool full=false);
pari edge_homology_relations(Triangulation *manifold);
pari cusp_homology_relations(Triangulation *manifold);
pari holonomies(Triangulation *manifold); 
pari get_parity_matrix(Triangulation* manifold);
pari compute_flattening(Triangulation* manifold, bool with_parity, bool filled, bool with_signs=false);

void check_equations(Triangulation* manifold, bool with_signs=false);

//struct ShapeInversion; 
//pari log_with_history(pari z, ShapeInversion *history);
//pari log_one_minus_z_with_history(pari z, ShapeInversion *history);

void set_accurate_shapes(Triangulation* manifold, const pari& shapes);
void clear_accurate_shapes(Triangulation* manifold);
pari get_accurate_cusp_shapes(const Triangulation* manifold); 

#endif
