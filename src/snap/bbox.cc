/*
** Copyright (C) 2004 Oliver A. Goodman <oag@ms.unimelb.edu.au>
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
#include "bbox.hh"

using std::list;
using std::ostream;

bbox::bbox(list<Complex> const& pl)
{
  if (!pl.size()) {
    ll = Zero; 
    ur = Zero; 
    return;
  }
  ll = pl.front();
  ur = ll; 
  include(pl); 
}

void bbox::include(const Complex& pt)
{
  if (pt.imag > top()) ur.imag = pt.imag;
  if (pt.imag < bottom()) ll.imag = pt.imag;
  if (pt.real > right()) ur.real = pt.real;
  if (pt.real < left()) ll.real = pt.real;
}

void bbox::include(const list<Complex>& pl)
{
  list<Complex>::const_iterator it;
  for (it = pl.begin(); it!=pl.end(); it++) include(*it); 
}

void bbox::expand(double x)
{
  ll -= Complex(x,x); 
  ur += Complex(x,x); 
}

// rows, columns. 
bbox bbox::subbox(int nr, int nc, int r, int c, double aspect, double margin) const
{
  double ws = (right() - left())/nc; 
  double hs = (top() - bottom())/nr;

  double w = ws;
  double h = hs; 

  // reduce h or w to get aspect ratio correct. 
  if (h/w > aspect) h = aspect*w;
  else w = h/aspect; 

  // create margins. 
  h *= (1. - margin); 
  w *= (1. - margin);  

  Complex ll00(left()+0.5*(ws-w), top()-0.5*(hs+h));
  Complex ll = ll00 + Complex(c * ws, -r * hs);
  return bbox(ll, ll + Complex(w,h)); 
}

// Get the transformation that makes inside fit into this bbox. 
void bbox::get_fit(bbox const& inside, double& scale, Complex& translation) const
{
  double vscale = (top()-bottom())/(inside.top()-inside.bottom()); 
  double hscale = (left()-right())/(inside.left()-inside.right()); 

  scale = (vscale < hscale) ? vscale : hscale; // min(vscale, hscale); 

  // the transformation we want is:
  // inside.mid() -> origin, scale, then origin -> mid().
  // this is x-> scale*(x-inside.mid()) + mid();
  // we will actually do scale then translate: x -> scale*x + translation
  // so:
  translation = mid() - scale*inside.mid(); 
}

#if 0 // unfinished
void bbox::get_best_grid(bbox const& inside, int n, int& nr, int& nc) const
{
  int r, c; 
  double sc, best_sc = 0.;
  double asp = inside.aspect(); 
  Complex tr; 
  bbox sub; 
  for (r=1; r<=n; r++) {
    sub = subbox(r, c, 0, 0, asp, 0.);
    get_fit(inside, sc, tr);
    // Unfinished
  }
}
#endif

void bbox::rescale(double scale, Complex const& translation)
{
  ll = scale * ll + translation; 
  ur = scale * ur + translation; 
}

// Checks whether two intervals overlap when projected onto a circle. 

static bool circle_overlap(double l1, double h1, double l2, double h2)
{
  // First compute the difference in angle between the midpoints of
  // the two projected intervals. 

  double adiff = ((l1+h1) - (l2+h2))/2.0;
  while (adiff > PI) adiff -= TWO_PI; 
  while (adiff <-PI) adiff += TWO_PI;
  adiff = fabs(adiff); 

  double radsum = ((h1-l1) + (h2-l2))/2.0;
  return radsum > adiff;
}

bool overlap(bbox const& a, bbox const& b, bool wrap_imag)
{
  if (a.left() > b.right() || a.right() < b.left()) return false; 

  if (!wrap_imag) {
    return (a.bottom() <= b.top() && a.top() >= b.bottom()); 
  }

  return circle_overlap(a.bottom(), a.top(), b.bottom(), b.top()); 
}

ostream& operator << (ostream& out, bbox const& a)
{
  return out << "[" << a.left() << ',' << a.right() << "]x[" << 
    a.bottom() << ',' << a.top() << "]";
}




