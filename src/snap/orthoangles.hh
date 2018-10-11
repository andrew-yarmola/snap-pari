#ifndef _orthoangles_
#define _orthoangles_
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

#include "snappea/SnapPea.h"
#include "mark.hh"

class OrthoAngleSetup; 
class EdgeClass; 
class EdgeIterator;
class TriangleIterator;

class EdgeEnd {
  EdgeClass *edge;
  bool base_end;
public:
  EdgeEnd() : edge(0), base_end(false) {}
  EdgeEnd(EdgeClass* e, bool b_end) : edge(e), base_end(b_end) {}

  int order() const;
  int index() const;

  void first(const Triangulation* T);
  EdgeEnd& operator ++(); 
  EdgeEnd  operator ++(int)
  { EdgeEnd ee(*this); ++(*this); return ee; }
  bool done(const Triangulation* T) const;

  EdgeIterator edge_iterator() const; 
  friend ostream& operator << (ostream& out, EdgeEnd const& e);
};

bool get_edge_orthodistances(Triangulation* manifold, GroupPresentation* G, 
			     vector<Complex>& orth);
void print_oa_cycles(Triangulation* manifold, GroupPresentation* G);
void print_orthoangles(Triangulation* t, GroupPresentation* G);
bool torus_reduce(vector<Complex> const& b, Complex& z);
bool torus_reduce(const Complex b[], Complex& z);

class OrthoAngles {
  OrthoAngleSetup* OAS; 
public:
  OrthoAngles(Triangulation* t, GroupPresentation* G); 
  ~OrthoAngles();

  void get_holonomies(vector<Complex>& H) const; 
  void get_da_spec(EdgeClass *e, vector<Complex>& spec) const; 
  void get_da_spec(EdgeIterator const& it, da_spec& da, bool report=false) const;
  void get_da_spec(EdgeEnd const& ee, da_spec& da) const;
  void get_da_spec(EdgeEnd const& ee, vector<da_spec>& spec) const; 

  void print() const; 
  bool valid() const;
};

#endif

