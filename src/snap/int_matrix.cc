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
#include <cstdio>
#include "int_matrix.hh"

void int_matrix::set_dims(int r, int c)
{
  if (rows*cols != r*c) {
    delete[] m; 
    m = new int[r * c]; 
  }
  rows = r; cols = c; 
}

void int_matrix::copy(int_matrix const& M)
{
  set_dims(M.rows, M.cols); 

  int i, n=rows*cols; 
  for (i=0; i<n; i++) {
    m[i] = M.m[i]; 
  }
}

void int_matrix::operator = (int_matrix const& M)
{
  if (&M==this) return; 
  copy(M); 
}
  
void int_matrix::operator = (int e)
{
  int i, n = rows*cols; 
  for (i=0; i<n; i++) m[i] = e; 
}
  
static int print_width(int n)
{
  if (n<0) n = -n; 
  int w = 2; 
  while (n>9) {
    w++; 
    n/= 10; 
  }
  return w; 
}

int field_width(int_matrix const& M)
{
  int i, n = M.rows*M.cols; 
  int w=2, wi; 
  for (i=0; i<n; i++) {
    wi = print_width(M.m[i]); 
    if (wi > w) w = wi; 
  }
  return w; 
}

void pretty_print(int_matrix const& M, int w, int row)
{
  if (w < 2) w = field_width(M); 
  char format[10];
  // printf("field width: %d\n", w); 
  sprintf(format, "%%%dd ", w); 
  int i, j; 
  if (row < 0) {
    for (i=0; i<M.rows; i++) {
      for (j=0; j<M.cols; j++) {
	printf(format, M[i][j]);
      }
      printf("\n"); 
    }
  } else {
    for (j=0; j<M.cols; j++) {
      printf(format, M[row][j]);
    }
  }
}

#if 0
int main()
{
  int_matrix M(2,3); 

  M[0][0] = 12; M[0][1] =  2; M[0][2] =  1; 
  M[1][0] =  4; M[1][1] = -3; M[1][2] =  0; 

  pretty_print(M);

  int_matrix L(2,2);
  L = 0; 
  pretty_print(L); 

  return 0; 
}
#endif


