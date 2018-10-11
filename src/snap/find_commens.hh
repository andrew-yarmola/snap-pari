#ifndef _find_commens_
#define _find_commens_
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


#include "i_triangulation.hh"

void find_commensurator(i_triangulation& m); 
void print_hidden_symmetries(i_triangulation& m, n_vector const& csm);
void check_commensurability(i_triangulation& m, n_vector const& csm, 
			    i_triangulation& n, n_vector const& csn);

void brute_force_hsymms(i_triangulation& m, int n, int N, 
			vector<int> const& ties, vector<int> const& areas);
bool validate_class_spec(vector<int> const& cs);

#endif
