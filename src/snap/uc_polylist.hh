#ifndef _uc_polylist_
#define _uc_polylist_
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
#include "vfunction.hh"
#include "eq_surface.hh"

using std::list; 

enum printstyle {
  surface, outline, wireframe 
};

enum marking {
  open_marking, closed_marking, unmarked 
};

class uc_vertex : public uc_point {
public:
  int mark_p, mark_n; 

  uc_vertex() {}

  uc_vertex(const uc_point& p, int mark = 0) : uc_point(p), mark_p(mark), mark_n(0) {}
  uc_vertex(const uc_point& p, int mp, int mn) : uc_point(p), mark_p(mp), mark_n(mn) {}

  void print() const; 
};

ostream& operator << (ostream& out, const uc_vertex& c);

/* 
   A polygon<uc_point> represents a polygon whose vertices are uc_point`s. 
   uc_point is required to have:
     a copy constructor uc_point(const uc_point&),
     a mma_print(uc_point, ostream&) fn.,
     and a global function crossing_point(uc_point a, double at, uc_point b, double bt).
   Some subset of the edges of uc_point are marked, which is to say, each edge
   has a flag associated with it. A polygon<uc_point> can be constructed by 
   pushing uc_point`s onto it using push_back(uc_point). By default this marks all but
   the implicit closing edge. The operation can thus be considered to be
   building up a polygonal curve. The close() method marks the closing edge. 
   The mma_print() function, by default, prints only the marked edges. 
   With an argument of "wireframe" it prints all the edges regardless of
   whether they are marked or not. With an argument of "surface" it prints 
   the polygon as a filled 2D polygon (which might be either in the plane 
   or in 3-space).

   Transform applies a function-object uc_point operator()(uc_point) to the vertices
   of the polygon, replacing them by their image vertices. 

   Clip evaluates a function-object double operator()(uc_point) on the vertices
   of the polygon. Vertices which have a positive value are retained while
   those with a negative value are discarded. A new vertex is created on
   each edge which goes between positive and negative values, and a new
   edge is added between consecutive new vertices. 
 */ 

class uc_polygon : public list<uc_vertex> {
public:
  uc_polygon() {}
  uc_polygon(const unary_vfunction<double, uc_point>& F, double lo, double hi, double step);

  uc_polygon& operator = (const uc_polygon& p)
    { list<uc_vertex>::operator = (p); return *this; }

  void push_back(const uc_point& p, marking m = open_marking);
  void push_back(const uc_vertex& p);

  void close();
  void open();

  void mma_print(printstyle s = outline, ostream& out = cout) const;
  void print() const; 

  void transform(const unary_vfunction<uc_point,uc_point>& F); 
  int clip(const unary_vfunction<uc_point,double>& F, int mark_new = 0,
	   const ternary_vfunction<uc_point,uc_point,int,uc_point>* G = 0, double eps = 1e-5);

  uc_polygon crossings(unary_vfunction<uc_point,double> const& F) const;
  uc_polygon corners() const; 
};

bool operator == (uc_polygon const&, uc_polygon const&); 

/* 
   A uc_polylist represents a collections of polygons with vertex type <uc_point>. 
   (No sharing of vertices is currently implemented.)
   uc_point`s are required to have:
     uc_point(const uc_point&) copy constructor,
     mma_print(uc_point, ostream&),
     gv_print(uc_point, ostream&),
     double sort_value(uc_point),
     bool operator == (uc_point a, uc_point b), 
     uc_point crossing_point(uc_point a, double at, uc_point b, double bt).
 */


class uc_polylist : public list<uc_polygon> {
public:
  uc_polylist() {}
  uc_polylist(const binary_vfunction<double, double, uc_point>& F, 
	   double xlo, double xhi, double xstep, 
	   double ylo, double yhi, double ystep);

  // Try not to use this, it will be very slow. 
  uc_polylist& operator = (const uc_polylist& p)
    { list< uc_polygon >::operator = (p); return *this; }

  void gv_print(printstyle s = surface, ostream& out = cout) const;
  void mma_print(printstyle s = surface, ostream& out = cout) const;
  void print() const;

  void add_boundary(const uc_polylist& p);
  void add_mesh(const binary_vfunction<double, double, uc_point>& F, 
	   double xlo, double xhi, double xstep, 
	   double ylo, double yhi, double ystep);

  void boundary()
    { iterator e = --end(); add_boundary(*this); erase(begin(), ++e); }
  void transform(const unary_vfunction<uc_point,uc_point>& F); 
  int clip(const unary_vfunction<uc_point,double>& F, int mark_new = 1,
	   const ternary_vfunction<uc_point,uc_point,int,uc_point>* G = 0, double eps = 1e-5);
  
  void accumulate_edges(const uc_polygon& p);
private:
  void gv_vect_print(bool closed, ostream& out = cout) const;
  void gv_off_print(ostream& out = cout) const;
};


#if 0
// polylist conversion class
// usage: 
//
//  polylist<F> from; 
//  polylist<T> to; 
//
//  to = convert_polylist<F,T>(from);
//
//  unary_vfunction<F,T> fn;
//
//  to = convert_polylist<F,T>(from, fn);
//

template <class F, class T>
class convert_polylist : public polylist<T> {
public:
  convert_polylist(const polylist<F>& m); 
  convert_polylist(const polylist<F>& m, const unary_vfunction<F,T>& H);

  static polygon<T> convert_polygon(const polygon<F>& p);
  static polygon<T> convert_polygon(const polygon<F>& p, const unary_vfunction<F,T>& H);
};

template <class F, class T>
polygon<T> convert_polylist<F,T>::convert_polygon(const polygon<F>& p)
{
  polygon<T> out;
  polygon<F>::const_iterator it; 
  for (it = p.begin(); it != p.end(); ++it)
    out.push_back(pvertex<T>(T(*it), (*it).mark_p)); 

  return out;
}

template <class F, class T>
convert_polylist<F,T>::convert_polylist(const polylist<F>& m)
{
  list< polygon<F> >::const_iterator it; 

  for (it = m.begin(); it != m.end(); it++)
    push_back(convert_polygon(*it)); 
}

template <class F, class T>
polygon<T> convert_polylist<F,T>::convert_polygon(const polygon<F>& p, const unary_vfunction<F,T>& H)
{
  polygon<T> out;
  polygon<F>::const_iterator it; 
  for (it = p.begin(); it != p.end(); ++it)
    out.push_back(pvertex<T>(H(*it), (*it).mark_p)); 

  return out;
}

template <class F, class T>
convert_polylist<F,T>::convert_polylist(const polylist<F>& m, const unary_vfunction<F,T>& H)
{
  list< polygon<F> >::const_iterator it; 

  for (it = m.begin(); it != m.end(); it++)
    push_back(convert_polygon(*it, H)); 
}
#endif

template <class Bbox>
void uc_find_box(const uc_polylist& plist, Bbox& box)
{
  list<uc_polygon>::const_iterator lpci;
  list<uc_vertex>::const_iterator lvci;
  // Iterate over polygons. 
  for (lpci = plist.begin(); lpci != plist.end(); lpci++) {
    // Iterate over vertices.
    for (lvci = (*lpci).begin(); lvci != (*lpci).end(); lvci++) {
      box.include(*lvci);
    }
  }
}


class uc_edge {
public:
  uc_point a, b; 
  int mark_p; 
  int evpi_a, evpi_b; // edge_vertex_ptr index of a and b. 

  uc_edge() {}
  uc_edge(const uc_point& aa, const uc_point& bb, int mn)
    : a(aa), b(bb), mark_p(mn), evpi_a(0), evpi_b(0) {}
};

class uc_end_ptr {
public:
  uc_edge *ptr; 
  int end; // one = a, zero = b

  const uc_point& the_point() const { return (end) ? ptr->a : ptr->b; }
  const uc_point& other_point() const { return (end) ? ptr->b : ptr->a; }
  double sort_value() const { return ::sort_value(the_point()); }
  int the_end() const { return (end) ? ptr->evpi_a : ptr->evpi_b; }
  int other_end() const { return (end) ? ptr->evpi_b : ptr->evpi_a; }

  uc_vertex end_vertex() { return uc_vertex(other_point(),ptr->mark_p); }
};

inline int operator < (const uc_end_ptr& a, const uc_end_ptr& b)
{ return a.sort_value() < b.sort_value(); }

#endif
