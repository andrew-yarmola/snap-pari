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
#include "pariwrap.hh"

#ifdef PARI_2_2_OR_LATER
#include "pari_oldnames.hh"
#endif

static initpari ip; 

void show_gc_off()
{ 
  printf("gc_off = %ld\n", gc_off); 
}

pari get_element(pari const& p, int n)
{
  return p[n];
}

bool find_smallest_in_last_row(pari& mat, int& col)
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

bool reduce_other_cols(pari& mat, int col)
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

void swap_cols(pari& mat, int c1, int c2)
{
  if (c1==c2) return; 
  pari col1 = mat[c1]; 
  mat[c1] = mat[c2]; 
  mat[c2] = col1; 
}

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

  return mat[c][row] == pONE;
}

int main()
{
  printf("Number of links at start of main()\n"); 
  pari::count_links(); 
  pari::print_all(); 
  show_gc_off(); 

  // Doing high-precision arithmetic is easy. 
  
  pari a = real(2); // Pari real (floating point), default precision 
  pari b = sqrt(a); // Floating point result, default precision
  
  // Print the results. We should be using C++ output streams 
  // and overloading << but, at the time of writing, Pari didn't 
  // provide us with a function for converting GENs to char*s. 

  printf("\nThe square root of\n");
  a.print();
  printf("is\n");
  b.print();

  // Make, modify, and print a 3-element row vector. This 
  // uses the overloaded [] operator. 

  printf("\n"); 
  pari vec = rvector(3); 
  vec[0] = 4; 
  vec[1] += vec[0]; 
  vec[2] = get_element(vec,0); 
  printf("A vector: "); 
  vec.print();

  pari x = real(100);
  x = real(200);

  pari mat = matrix(2,2); 
  mat[0][0] = vec[0] + pONE; 
  mat[0][1] = vec[0];
  mat[1][0] = vec[1] - pONE; 
  mat[1][1] += mat[0][0]; 
  printf("A matrix: "); 
  mat.print(); 
  mat[0][0] += pTWO; 
  printf("Modified: ");
  mat.print(); 

  printf("element 1,1 = ");
  pari z(get_element(mat[1],1));
  z.print();
  pari e01(mat[0][1]); // use pari copy on pari_slot
  printf("element 0,1 = "); 
  e01.print(); 

  // Creating a polynomial involves an extra step. 

  printf("\n"); 
  printf("Creating a polynomial...\n"); 
  pari p = polynomial(2); // Zero polynomial with space for a quadratic. 
  printf("length %ld, effective length %ld\n", length(p), lgeff(p)); 
  p[0] = integer(2);
  p[1] = integer(-3);
  p[2] = integer(0); 
  printf("length %ld, effective length %ld\n", length(p), lgeff(p)); 

  // Make sure the polynomial knows how long it is (and whether nonzero). 
  set_poly(p); 
  printf("length %ld, effective length %ld\n", length(p), lgeff(p)); 
  printf("The polynomial: "); 
  p.print();

  // Let's take some rubbish off the stack. This is normally 
  // done automatically when garbage_count > 1000 and when garbage
  // collection is next enabled. 

  printf("\n"); 
  show_gc_off(); 
  printf("top %lX avma %lX bot %lX\n", top, avma, bot);
  printf("garbage_count=%ld\n", garbage_count); 
  printf("calling collect_garbage()\n"); 
  pari::collect_garbage(); 
  printf("top %lX avma %lX bot %lX\n", top, avma, bot);
  printf("garbage_count=%ld\n", garbage_count); 
  show_gc_off(); 

  // Let us also resize the stack. This relocates the whole thing. 

  printf("resizing the stack to 6Mb\n"); 
  pari::resize_stack(6000000);
  printf("top %lX avma %lX bot %lX\n", top, avma, bot);

  printf("\n"); 
  pari::count_links(); 

  // Our polynomial should still be there! 

  printf("\nThe polynomial again: "); p.print(); 

  printf("now calling clone_all()\n"); 
  pari::clone_all(); 

  printf("\n"); 
  pari::count_links(); 
  printf("\nThe polynomial again: "); p.print(); 
  printf("The vector again: "); vec.print();
  vec[1] = pONE; 
  printf("The vector after modification: "); vec.print();
  pari::count_links(); 

  printf("assigning to a clone\n");
  x = integer(365); // assign to a clone. 
  pari::count_links(); 

  printf("\n");
  printf("x = "); x.print();
  printf("attempt to get x[0]\n");
  x[0].print(); 
  printf("attempt to assign x[0] = 1\n");
  x[0] = pONE; 
  printf("x = "); x.print();

  printf("\n");
  printf("vec = "); vec.print();
  printf("attempt to get vec[5]\n");
  vec[5].print(); 
  printf("attempt to assign vec[5] = 1\n");
  vec[5] = pONE; 
  printf("vec = "); vec.print();
  
  printf("\n"); 
  pari mx = p_flisexpr("[1,1;-1,1]");
  printf("Matrix: "); mx.print();
  pari res = ker_mod_p(mx, pTWO);
  printf("Mod 2 kernel: "); res.print();

  mx[0][1] = 2; 
  mx[1][1] = 6; 
  printf("Matrix: "); mx.print();
  printf("\ndoing make_last_row_en\n"); 
  if (!make_last_row_en(mx)) {
    printf("make_last_row_en failed!\n");
  } 
  printf("Matrix: "); mx.print();

#if 0
  printf("type of mx %d\n", mx.type()); 
  pari cv = rvector(2); 
  printf("type of row vector %d\n", cv.type());

  printf("3 \\/ -2\n");
  pari r = gdivround(integer(3),integer(-2)); 
  r.print(); 
#endif 

  printf("copy of element 0,1 = ");
  pari y(e01);
  y.print();

  pari bignum = p_flisexpr("817482163498729328749");
  pari rem = bignum % pTWO; 
  printf("remainder = "); rem.print(); 
}


/* $Id: example.cc,v 1.2 2009/12/03 21:52:43 matthiasgoerner Exp $ */ 
