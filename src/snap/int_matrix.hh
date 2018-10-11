#ifndef _int_matrix_
#define _int_matrix_
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


class int_matrix {
public:
  int *m;
  int rows, cols; 

  int_matrix() : m(0), rows(0), cols(0) {}
  int_matrix(int r, int c) : rows(r), cols(c) 
    { m = new int[r*c]; }
  int_matrix(int_matrix const& M) { copy(M); }
  ~int_matrix() { delete[] m; }

  void operator = (int_matrix const& M);

  void set_dims(int r, int c);
  void copy(int_matrix const& M); 

  int* operator [] (int r) { return &m[r*cols]; }
  const int* operator [] (int r) const { return &m[r*cols]; }

  void operator = (int e);
};

void pretty_print(int_matrix const& M, int w=0, int row=-1);
int field_width(int_matrix const& M);

#endif
