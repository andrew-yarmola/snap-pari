#ifndef _lines_
#define _lines_
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


#include "snappea/Moebius_transformations.h"

using std::ostream;

class line {
public:
  Complex end[2];

  line() {} 
  line(const Complex& a, const Complex& b) { end[0] = a; end[1] = b; }
  line(const MoebiusTransformation& mt); // mt * [0,infinity]
  void print() const;

  friend bool operator == (line const&, line const&); 
  friend bool operator < (line const&, line const&); 
};

void reflection_in_line(const line& l, SL2CMatrix refl);

line orthogonal_through_origin(const line& l);

line ortholine(const line& axis1, const line& axis2);

void find_ortholine(const SL2CMatrix line1, const SL2CMatrix line2, SL2CMatrix ortholine);
void find_ortholine(const line& axis1, const line& axis2, SL2CMatrix ortholine_matrix);

Complex parabolic_fix(const SL2CMatrix m);

void find_line_matrix(const line& l, SL2CMatrix line_matrix);

Complex orthodist(const SL2CMatrix line1, const SL2CMatrix line2);
Complex orthodist(const line& line1, const line& line2);

Complex cosh_orthodist(const line& axis1, const line& axis2);
Complex cosh_orthodist(const SL2CMatrix line1, const SL2CMatrix line2);

int fixed_points(const SL2CMatrix m, line& fix, FILE *fp = NULL);
int real_fixed_points(const SL2CMatrix m, line& fix, FILE *fp = NULL);

line operator * (const MoebiusTransformation& m, const line& l);
bool operator == (const line& a, const line& b); 
int compare(const line& a, const line& b, double eps=1e-10);

inline ostream& operator << (ostream& out, const line& l)
{
  return out << '[' << l.end[0] << ',' << l.end[1] << ']'; 
}

#endif
