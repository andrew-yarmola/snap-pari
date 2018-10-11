#ifndef _bbox_
#define _bbox_
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

#include "snappea/complex.h"
#include <list>

class bbox {
  Complex ll,ur;
public:
  bbox() : ll(Zero), ur(Zero) {}
  bbox(Complex const& LL, Complex const& UR) : ll(LL), ur(UR) {}
  bbox(double l, double r, double b, double t) : ll(l,b), ur(r,t) {}
  bbox(std::list<Complex> const& pl); 

  double top() const { return ur.imag; }
  double bottom() const { return ll.imag; }
  double left() const { return ll.real; }
  double right() const { return ur.real; }
  Complex mid() const { return 0.5*(ll+ur); }
  double aspect() const { return (top()-bottom())/(right()-left()); }

  void include(const Complex& pt);
  void include(const std::list<Complex>& pl);
  void expand(double x); 

  bbox subbox(int nr, int nc, int r, int c, double aspect=1., double margin=.1) const;
  void get_fit(bbox const& inside, double& scale, Complex& translation) const;
  void rescale(double scale, Complex const& translation); 

  friend bool overlap(bbox const& a, bbox const& b, bool wrap_imag); 
  friend std::ostream& operator << (std::ostream& out, bbox const& a); 
  friend bbox operator + (bbox const& b, Complex const& z)
  { return bbox(b.ll+z, b.ur+z); }
};

#endif
