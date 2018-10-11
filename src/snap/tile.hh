#ifndef _tile_
#define _tile_
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


#include "snappea/complex.h"
#include <list>
#include <vector>

#include "vfunction.hh"

using std::vector;
using std::list;

// half_plane(Complex z) gives the closed half plane through origin with inward normal z
// half_plane(Complex n, Complex pt) gives the closed half-plane through p with normal n

class half_plane : public unary_vfunction<Complex,bool> {
public:
  half_plane(const Complex& n = I) // default is upper half plane
    : normal(n), disp(0.0), open(false) {}
  half_plane(const Complex& n, double d, bool o = false)
    : normal(n), disp(d), open(o) {}
  half_plane(const Complex& n, const Complex& pt, bool o = false)
    : normal(n), disp(n.real * pt.real + n.imag * pt.imag), open(o) {}

  double evaluate(const Complex& t) const
    { return normal.real * t.real + normal.imag * t.imag - disp; }
  bool operator () (const Complex& t) const
    { return open ? evaluate(t) > 1e-10 : evaluate(t) >= 0.; }
private:
  Complex normal;
  double disp;
  bool open; 
};

class translation : public unary_vfunction<Complex,Complex> {
public:
  translation() : t(Zero) {}
  translation(const Complex& tr) : t(tr) {}

  Complex operator () (const Complex& p) const
    { return p + t; }
  translation inverse() const
    { return translation(-t); }
private:
  Complex t;
};

class tiler {
public:
  tiler() {}
  tiler(const Complex& a, const Complex& b); 
  tiler(const vector<half_plane>& r, const vector<translation>& m) 
    : region(r), maps(m) {}

  int which_region(const Complex& z) const;
  bool reduced(const Complex& z) const
    { return which_region(z) == -1; }
  bool reduce_once(Complex& z) const;
  void reduce(Complex& z) const;

  void grow(const Complex& z, list<Complex>& leaves); 
  list<Complex> tile(const Complex& bp, double radius); 
  list<Complex> tile(const Complex& bp, int layers); 

private:
  vector<half_plane> region;
  vector<translation> maps;
};

list<Complex> filter(const list<Complex>& pts, vector<half_plane>& region);

#endif


