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
#include "tile.hh"

/* The following constructor makes a tiler for a parallelogram with
   sides a and b centered at the origin. */ 

tiler::tiler(const Complex& a, const Complex& b)
: region(4), maps(4)
{
  Complex rot = ((b/a).imag > 0.0) ? I : -I; 

  region[0] = half_plane(-rot * a, -0.5 * b, false);
  region[1] = half_plane(rot * a, 0.5 * b, true); 
  region[2] = half_plane(rot * b, -0.5 * a, false);
  region[3] = half_plane(-rot * b, 0.5 * a, true); 
  
  maps[0] = translation(-b); 
  maps[1] = translation(b); 
  maps[2] = translation(-a); 
  maps[3] = translation(a);
}

list<Complex> filter(const list<Complex>& pts, vector<half_plane>& region)
{
  list<Complex> result; 
  int i,n = region.size();
  list<Complex>::const_iterator pt;
  for (pt = pts.begin(); pt != pts.end(); pt++) {
    for (i=0; i<n; i++)
      if (!region[i](*pt)) break;
    if (i==n) // did not break, so point is in all regions
      result.push_back(*pt);
  }
  return result; 
}

int tiler::which_region(const Complex& z) const
{
  int i, n = region.size(); 
  for (i=0; i<n; i++)
    if (region[i](z)) break;
  if (i==n) return -1; 
  return i;
}

bool tiler::reduce_once(Complex& z) const
{
  int r = which_region(z);
  if (r==-1) return false; 
  double v = region[r].evaluate(z);
  z = (maps[r].inverse())(z); 
  if (v > 1.) return true; 

  // maybe we didn't get any closer? 
  r = which_region(z);
  if (region[r].evaluate(z) + 1e-10 > v) return false;

  return true;
}

void tiler::reduce(Complex& z) const
{ 
  int limit = 100; 
  while (reduce_once(z) && --limit > 0); // empty loop body
  if (limit==0) fprintf(stderr, "limit reached in tiler::reduce\n"); 
} 

void tiler::grow(const Complex& z, list<Complex>& leaves)
{
  int i; 
  Complex w;
  for (i=0; i<maps.size(); i++) {
    w = maps[i](z);
    if (which_region(w)==i) leaves.push_back(w);
  }
}

list<Complex> tiler::tile(const Complex& bp, int layers)
{
  list<Complex> tiles, leaves, new_leaves; 
  leaves.push_back(bp);

  list<Complex>::const_iterator it;
  while (layers-- > 0) {
    for (it = leaves.begin(); it != leaves.end(); it++)
      grow(*it, new_leaves);
    tiles.splice(tiles.end(), leaves);
    leaves.splice(leaves.end(), new_leaves);
  }
  tiles.splice(tiles.end(), leaves);
  return tiles; 
}

list<Complex> tiler::tile(const Complex& bp, double radius)
{
  list<Complex> tiles, leaves, new_leaves; 
  new_leaves.push_back(bp);

  list<Complex>::const_iterator it;
  while (true) {
    for (it = new_leaves.begin(); it != new_leaves.end(); it++)
      if (complex_modulus(*it) <= radius) 
	leaves.push_back(*it);
    new_leaves.erase(new_leaves.begin(), new_leaves.end());

    if (leaves.size()==0) break;

    for (it = leaves.begin(); it != leaves.end(); it++)
      grow(*it, new_leaves);
    tiles.splice(tiles.end(), leaves);
  }
  return tiles;
}

#if 0
int main()
{
  Complex a(1.7, 0.3), b(0.1, 0.9);

  vector<half_plane> pgram(4);

  Complex rot = ((b/a).imag > 0.0) ? I : -I; 

  pgram[0] = half_plane(-rot * a, Zero);
  pgram[1] = half_plane(rot * a, b); 
  pgram[2] = half_plane(rot * b, Zero);
  pgram[3] = half_plane(-rot * b, a); 
  
  vector<translation> mapss(4);

  mapss[0] = translation(-b); 
  mapss[1] = translation(b); 
  mapss[2] = translation(-a); 
  mapss[3] = translation(a);

  tiler a_tiler = tiler(pgram, mapss);

  Complex z(-1.5, 2.5);

  printf("starting with complex number: "); print(z); printf("\n"); 

  a_tiler.reduce_once(z); 
  printf("reduced once: "); print(z); printf("\n"); 

  a_tiler.reduce(z);
  printf("reduced: "); print(z); printf("\n\n"); 


  list<Complex> lattice1 = a_tiler.tile(z, 2);
  list<Complex>::const_iterator i;
  printf("tiling two layers around this\n"); 
  for (i=lattice1.begin(); i!=lattice1.end(); i++) {
    print(*i); printf("\n"); 
  }

  return 0; 
}
#endif
