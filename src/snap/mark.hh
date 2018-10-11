#ifndef _MARK_H_
#define _MARK_H_
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

#include <iostream>
#include "snappea/complex.h"

class Mark {
  int num;

  void get_ml(int H[]) const;
public:
  Mark(int val=0) : num(val) {}
  Mark(int face, int me, int lo)
  : num(face*0x100 + (me+8)*0x10 + lo+8) {}
  Mark(int face, Mark const& mark)
  : num(face*0x100 + (mark.num & 0xff)) {}

  int value() const { return num; }

  Mark partner() const { return Mark(num + 0x100); } // Fix this! 
  bool set_to(const char* s); 
  int end_index() const { return num/0x100; }
  int lensq() const;
  Complex holo(const Complex H[2]) const; 
  Mark operator - () const { return Mark(num + 0x110 - 2*(num & 0xff)); }

  friend bool operator == (Mark const& a, Mark const& b)
  { return a.num == b.num; }
  friend bool operator != (Mark const& a, Mark const& b)
  { return a.num != b.num; }
  friend bool operator < (Mark const& a, Mark const& b)
  { return a.lensq() < b.lensq(); }
  friend Mark operator + (Mark const& a, Mark const& b)
  { return Mark(a.num+b.num-0x88); }

  friend std::ostream& operator << (std::ostream& out, Mark m);
};

class da_spec {
  Complex d, a; 
  Mark m;
public:
  da_spec() {}
  da_spec(Complex const& _d, Complex const& _a, Mark _m)
    : d(_d), a(_a), m(_m) {}

  Complex const& distance() const { return d; }
  Complex const& angle() const { return a; }
  Mark const& mark() const { return m; }

  bool read(FILE* fp=stdin); 

  friend std::ostream& operator << (std::ostream& out, da_spec const& s);
  friend bool operator == (da_spec const& a, da_spec const& b)
  { return a.a == b.a && a.d == b.d; } // ignores marking
  friend bool operator != (da_spec const& a, da_spec const& b)
  { return !(a==b); }
};

#endif

