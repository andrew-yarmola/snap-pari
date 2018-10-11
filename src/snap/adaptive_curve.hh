#ifndef _adaptive_curve_
#define _adaptive_curve_
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


#include <list>
#include <vector>
#include "snappea/complex.h"

/* Takes a function object F defining a map from double to Complex. 
   Returns a list of Complexes obtained by evaluating F(x) at a monotone 
   sequence of x values starting with xlo, and ending with xhi. Values 
   are chosen such that successive x values differ by at most d_step 
   while successive values of F(x) are at most distance r_step apart in 
   the complex plane. F must be continuous or this function will fail. 
*/ 

template <class func>
list<Complex> adaptive_curve(func F, double xlo, double xhi, double d_step, double r_step)
{
  list<Complex> out; 
  bool reverse = xlo > xhi; 
  if (reverse) {
    xlo = -xlo; 
    xhi = -xhi; 
  }
  double x = xlo, last_x; 
  Complex z = F(reverse ? -x : x), last_z;
  double end_computed = false; 
  vector<double> waiting_x;
  vector<Complex> waiting_z; 

  r_step *= r_step; // Square it, for comparison with complex_modulus_squared. 

  while (true) {
    out.push_back(z);
    if (end_computed && !waiting_z.size()) break; 

    last_x = x;
    last_z = z; 
    if (waiting_z.size()) {
      x = waiting_x.back();
      waiting_x.pop_back();
      z = waiting_z.back();
      waiting_z.pop_back();
    } else {
      x += d_step; 
      if (x >= xhi) {
	x = xhi; 
	end_computed = true; 
      }
      z = F(reverse ? -x : x);
    }

    while (complex_modulus_squared(z - last_z) > r_step) {
      waiting_x.push_back(x);
      waiting_z.push_back(z);
      x = (x + last_x)/2.0;
      z = F(reverse ? -x : x);

      if (waiting_z.size() > 20) {
	cout << "Function in adaptive_curve appears discontinuous on [" 
	     << xlo << ", " << xhi << "]." << endl; 
	return out; 
      }
    }
  }
  return out; 
}

#endif
