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
#include "uc_polylist.hh"

#include <algorithm>
#include <assert.h>
#include <vector>

#define trace false
#define SORTVAL_EPS 1e-6

ostream& operator << (ostream& out, const uc_vertex& c)
{ 
  return out << (uc_point&)c << ' ' << c.mark_p << ' ' << c.mark_n; 
}

void uc_vertex::print() const
{ 
  cout << *this; 
}


uc_polygon::uc_polygon(const unary_vfunction<double,uc_point>& F, 
		       double lo, double hi, double step)
{
  double t;
  for (t=lo; t<=hi; t+=step)
    push_back(uc_vertex(F(t),1));
  open();
}

void uc_polygon::push_back(const uc_point& p, marking m)
{
  int mark = (m==closed_marking || (m==open_marking && size()!=0));
  list<uc_vertex>::push_back(uc_vertex(p,mark)); 
}

void uc_polygon::push_back(const uc_vertex& p)
{
  list<uc_vertex>::push_back(p); 
}

void uc_polygon::open()
{
  front().mark_p = 0;
}

void uc_polygon::close()
{
  front().mark_p = 1;
}

void uc_polygon::transform(const unary_vfunction<uc_point,uc_point>& F)
{
  iterator i;
  for (i=begin(); i!=end(); i++)
    *i = uc_vertex(F(*i),(*i).mark_p);
}

uc_polygon uc_polygon::corners() const
{
  uc_polygon C;
  const_iterator i,j;

  // Loop once around the polygon. j stays one step ahead of 
  // i and returns to the start of the polygon when i is
  // at the end. 
  for (i = begin(), j = i, j++; i!=end(); i++, j++) {
    if (j==end()) j = begin(); 
    if ((*i).mark_p != (*j).mark_p)
      C.push_back(uc_vertex(*i, (*i).mark_p, (*j).mark_p));
  }
  return C; 
}

void uc_polygon::print() const
{
  const_iterator b = begin(); 
  cout << "[";
  while (b != end()) {
    (*b++).print();
    if (b != end()) cout << ", ";
  }
  cout << "]\n";
}

void uc_polylist::print() const
{
  const_iterator i;
  for (i=begin(); i!=end(); i++)
    (*i).print();
}

// The result returned by clip is as follows: 
// Bit 0 is set if there was anything inside F.  
// Bit 1 is set if there was anything outside F.

int uc_polygon::clip(const unary_vfunction<uc_point,double>& F, int mark_new, 
		     const ternary_vfunction<uc_point,uc_point,int,uc_point>* G, double eps)
{
  if (size()==0) return 0; // Can we get out really quick? 

  // Can we get out pretty quick? 
  const_iterator it = begin();
  if (F(*it) > 0.0) { 
    while (it != end() && F(*it) > 0.0) ++it; 
    if (it==end()) return 1; // Everything is inside. 
  } else if (F(*it) < 0.0) {
    while (it != end() && F(*it) < 0.0) ++it; 
    if (it==end()) {
      erase(begin(), end()); 
      return 2; // Everything is outside. 
    }
  }
    
  uc_point cp; 
  uc_vertex prev_pt, pt, crossing_pt; 
  double prev_value, value; 

  uc_polygon n_poly; 

  prev_pt = back(); 
  prev_value = F(prev_pt); 

  // G = 0; // For debugging only. 

  for (it = begin(); it != end(); it++) {
    pt = *it; 
    value = F(pt); 

    if (value > 0.0) {
      if (prev_value < 0.0) { // entering
	cp = G ? (*G)(prev_pt,pt,pt.mark_p) : crossing_point(prev_pt,prev_value,pt,value);
	crossing_pt = uc_vertex(cp, mark_new);
	n_poly.push_back(crossing_pt);
      }
      n_poly.push_back(pt); 

    } else if (value == 0.0) {
      n_poly.push_back(pt); 
      if (prev_value < 0.0) // entering
	n_poly.back().mark_p = mark_new;

    } else /* value < 0.0 */ {
      if (prev_value > 0.0) { // leaving
	cp = G ? (*G)(prev_pt,pt,pt.mark_p) : crossing_point(prev_pt,prev_value,pt,value);
	crossing_pt = uc_vertex(cp, pt.mark_p);
	n_poly.push_back(crossing_pt); 
      }
    }
    prev_pt = pt; 
    prev_value = value; 
  }

  swap(n_poly); 
  if (size()==1) erase(begin());

  return 3;
}

uc_polygon uc_polygon::crossings(const unary_vfunction<uc_point,double>& F) const
{
  uc_polygon result; 

  if (size()==0) return result; // Can we get out really quick? 

  // Can we get out pretty quick? 
  const_iterator it = begin();
  if (F(*it++) >=0.) { 
    while (it != end() && F(*it) >=0.) ++it; 
    if (it==end()) return result; // Everything is inside. 
  } else {
    while (it != end() && F(*it) < 0.) ++it; 
    if (it==end()) {
      return result; // Everything is outside. 
    }
  }

  uc_vertex prev_pt, pt, crossing_pt; 
  double prev_value, value; 

  prev_pt = back(); 
  prev_value = F(prev_pt); 

  for (it = begin(); it != end(); it++) {
    pt = *it; 
    value = F(pt); 

    if (value >=0.) {
      if (prev_value < 0.) { // entering
	crossing_pt = uc_vertex(crossing_point(prev_pt,prev_value,pt,value),pt.mark_p);
	result.push_back(crossing_pt);
      }

    } else if (value < 0.) {
      if (prev_value >=0.) { // leaving
	crossing_pt = uc_vertex(crossing_point(prev_pt,prev_value,pt,value),pt.mark_p);
	result.push_back(crossing_pt); 
      }
    }
    prev_pt = pt; 
    prev_value = value; 
  }

  return result;
}

void uc_polygon::mma_print(printstyle s, ostream& out) const
{
  if (size() < 2) {
    out << "{}";
    return; 
  }
  const_iterator pi; 

  switch (s) {
  case surface:

    out << "Polygon[{"; 
    for (pi = begin(); pi != end(); pi++) {
      if (pi != begin()) out << ", "; 
      ::mma_print(*pi, out); 
    }
    out << "}]"; 
    break;

  case wireframe:

    out << "Line[{"; 
    for (pi = begin(); pi != end(); pi++) {
      if (pi != begin()) out << ", "; 
      ::mma_print(*pi,out); 
    }
    out << ", "; 
    ::mma_print(front(),out);
    out << "}]"; 
    break;

  case outline:

    {
      out << "{";

      bool doing_line = false; 
      bool done_any = false;
      uc_point prev = back();

      for (pi = begin(); pi != end(); pi++) {
	if ((*pi).mark_p) {
	  if (!doing_line) {
	    if (done_any) out << ",";
	    out << "Line[{"; 
	    ::mma_print(prev,out);
	    doing_line = true;
	    done_any = true; 
	  }
	  out << ","; 
	  ::mma_print(*pi,out); 
	} else {
	  prev = *pi; 
	  if (doing_line) {
	    out << "}]"; 
	    doing_line = false;
	  }
	}
      }
      if (doing_line) out << "}]";

      out << "}";
    }
    break;

  default:
    break;
  }
}

uc_polylist::uc_polylist(const binary_vfunction<double, double, uc_point>& F, 
  double x_lo, double x_hi, double x_step, double y_lo, double y_hi, double y_step)
{
  add_mesh(F, x_lo, x_hi, x_step, y_lo, y_hi, y_step); 
}

void uc_polylist::add_mesh(const binary_vfunction<double, double, uc_point>& F, 
  double x_lo, double x_hi, double x_step, double y_lo, double y_hi, double y_step)
{
  double x, y;

  uc_polygon p, empty;

  x_hi -= x_step/10.0; // don`t want x==x_hi accidentally due to roundoff error. 
  y_hi -= y_step/10.0; 

  for (x=x_lo; x<x_hi; x+=x_step) {
    for (y=y_lo; y<y_hi; y+=y_step) { 

      p.push_back(uc_vertex(F(x,y),false));
      p.push_back(uc_vertex(F(x+x_step,y),false));
      p.push_back(uc_vertex(F(x+x_step,y+y_step),false));
      p.push_back(uc_vertex(F(x,y+y_step),false)); 

      push_back(empty);
      back().swap(p);
    }
  }
}

void uc_polylist::transform(const unary_vfunction<uc_point,uc_point>& F)
{
  iterator i;
  for (i=begin(); i!=end(); i++)
    (*i).transform(F);
}

// The result returned by clip is as follows: 
// Bit 0 is set if there was anything inside F.  
// Bit 1 is set if there was anything outside F.
// Ie.
// 0 = nothing either side, 1 = everything inside
// 2 = everything outside,  3 = stuff on both sides. 

int uc_polylist::clip(const unary_vfunction<uc_point,double>& F, int mark_new,
		      const ternary_vfunction<uc_point,uc_point,int,uc_point>* G, double eps)
{
  int res = 0;
  iterator i=begin();
  while (i!=end()) {
    res |= (*i).clip(F,mark_new,G,eps);
    if ((*i).size()==0) 
      erase(i++);
    else ++i;
  }
  return res;
}

void uc_polylist::accumulate_edges(const uc_polygon& p)
{
  bool doing_line = false; 
  uc_point prev = p.back();
  uc_polygon::const_iterator pi; 

  for (pi = p.begin(); pi != p.end(); pi++) {
    if ((*pi).mark_p) {
      if (!doing_line) {
	push_back(uc_polygon());
	back().push_back(prev);
	doing_line = true;
      }
      back().push_back(*pi);
    } else {
      prev = *pi; 
      doing_line = false;
    }
  }
}

void uc_polylist::gv_vect_print(bool closed, ostream& out) const
{
  int v_count = 0;
  list< uc_polygon >::const_iterator lit;
  uc_polygon::const_iterator pli; 

  for (lit = begin(); lit != end(); lit++)
    v_count += (*lit).size();

  // VECT \n NPOLYLINES NVERTICES NCOLORS
  out << "{ VECT\n" << size() << ' ' << v_count << " 0\n"; 
    
  // NV[0] ... NV[NPOLYLINES-1] # number of vertices in each polygon
  for (lit = begin(); lit != end(); lit++)
    out << (closed ? -(*lit).size() : (*lit).size()) << ' ';
  out << '\n'; 

  // NC[0] ... NC[NPOLYLINES-1] # number of colors in each polygon
  for (lit = begin(); lit != end(); lit++)
    out << "0 ";
  out << '\n'; 
  
  // VERT[0] ... VERT[NVERTICES-1]
  for (lit = begin(); lit != end(); lit++)
    for (pli = (*lit).begin(); pli != (*lit).end(); pli++) {
      ::gv_print(*pli,out);
      out << '\n';
    }
  
  out << "\n}\n"; 
}

void uc_polylist::gv_off_print(ostream& out) const
{
  int v_count = 0;
  const_iterator it; 
  for (it = begin(); it != end(); it++) 
    v_count += (*it).size(); 
  
  uc_polygon::const_iterator pi;
  
  // { OFF \n NVERTICES  NFACES  NEDGES 
  out << "{ OFF\n" << v_count << ' ' << size() << " 0\n"; 
  
  // X[0]  Y[0]  Z[0] \n ... X[NVERTICES-1]  Y[NVERTICES-1]  Z[NVERTICES-1]
  for (it = begin(); it != end(); it++)
    for (pi = (*it).begin(); pi != (*it).end(); pi++) {
      ::gv_print(*pi,out); 
      out << '\n';
    }  

  // for each face: NV V[0] V[1] ... V[NV-1]
  v_count = 0; 
  for (it = begin(); it != end(); it++) {
    out << '\n' << (*it).size() << ' '; 
    for (pi = (*it).begin(); pi != (*it).end(); pi++)
      out << v_count++ << ' '; 
  }

  out << "\n}\n"; 
}

void uc_polylist::gv_print(printstyle s, ostream& out) const
{
  const_iterator it;

  switch (s) {
  case wireframe:

    gv_vect_print(true, out); 
    break;

  case outline:
    {
      uc_polylist outl;

      for (it = begin(); it != end(); ++it)
	outl.accumulate_edges(*it);

      outl.gv_vect_print(false, out); 
    }
    break;

  case surface:

    gv_off_print(out);
    break;

  default:
    break;
  }
}

void uc_polylist::mma_print(printstyle s, ostream& out) const
{
  const_iterator it;
  out << "{";

  for (it = begin(); it != end(); it++) {
    if (it != begin()) out << ",\n";
    (*it).mma_print(s, out);
  }
  out << "}\n";
}

// add_boundary makes a list of all the marked edges of p,
// then tries to assemble them into a set of polygonal loops and 
// segments. To create each loop or segment, we begin with an
// arbitrary unused edge (marking is set to zero on used edges). 
// In fact we start at one end of the edge, which gives it an
// orientation. The next edge is the one sharing a vertex with the
// other end of this edge. The previous edge is the one sharing
// a vertex with this end of the edge. We append the initial edge
// to the loop or segment we are creating. Then, while there is a
// next edge, and it is unused, we advance to it and append it to 
// the loop/segment. If the next edge is the initial edge we have
// a loop and are done. If there is a next edge, which is unmarked
// but which is not the initial edge, clearly the edges branch, 
// they do not form a collection only of loops and segments. 
// If there is no next edge we have a segment:
// we then return to the initial edge and prepend previous edges
// until there is no previous edge. We then prepend the start of the
// first (final previous) edge, unmarked so that the result is a
// polygonal segment rather than a loop. 

inline ostream& operator << (ostream& out, const uc_end_ptr& p)
{ return out << p.the_end() << ' ' << p.the_point() << ' ' << p.other_point(); }

static void get_edges(uc_polylist const& P, list< uc_edge >& el)
{
  uc_point prev;

  uc_polylist::const_iterator it; 
  uc_polygon::const_iterator pi; 
  for (it = P.begin(); it != P.end(); it++) {
    
    prev = (*it).back(); 
    for (pi = (*it).begin(); pi != (*it).end(); pi++) {
      if ((*pi).mark_p && prev != (uc_point&)*pi) // Avoid zero length edges.
	el.push_back(uc_edge(prev, *pi, (*pi).mark_p));
      prev = *pi;
    }
  }
}

// next_edge, goes to the other end of this edge, then
// finds the matching vertex. 

static bool next_edge(vector< uc_end_ptr >& epv, int& i)
{
  int j = epv[i].other_end();
  int k = i; 

  // find matching vertex for j
  double s = epv[j].sort_value(); 
  uc_point pt = epv[j].the_point();

  for (i = j+1; i < epv.size() && epv[i].sort_value() <= s + SORTVAL_EPS; i++) 
    {
      if (epv[i].the_point() == pt) return true; 
    }

  for (i = j-1; i >= 0 && epv[i].sort_value() >= s - SORTVAL_EPS; i--)
    {
      if (epv[i].the_point() == pt) return true;
    }     

  i = k; 
  return false;
}

// prev_edge, finds a vertex matching this one then goes to 
// the other end of its edge. 

static bool prev_edge(vector< uc_end_ptr >& epv, int& j)
{
  // find matching vertex for j
  double s = epv[j].sort_value(); 
  uc_point pt = epv[j].the_point();
  int k = j; 

  int i;
  for (i = j+1; i < epv.size() && epv[i].sort_value() <= s + SORTVAL_EPS; i++)
    if (epv[i].the_point() == pt) {
      j = epv[i].other_end();
      return true; 
    }

  for (i = j-1; i >= 0 && epv[i].sort_value() >= s - SORTVAL_EPS; i--)
    if (epv[i].the_point() == pt) {
      j = epv[i].other_end();
      return true;
    }

  j = k; 
  return false;
}

void uc_polylist::add_boundary(const uc_polylist& p)
{
  list< uc_edge > el;
  get_edges(p, el);

  int i, j; 

  // make vector of uc_end_ptrs
  vector< uc_end_ptr > evps(2 * el.size()); 
  list< uc_edge >::iterator it; 
  i = 0; 
  for (it = el.begin(); it != el.end(); it++) {
    evps[i].ptr = evps[i+1].ptr = &(*it); 
    evps[i++].end = 1;
    evps[i++].end = 0; 
  }

  // sort the uc_end_ptrs
  std::sort(evps.begin(), evps.end()); 

  // put in pointers from edges to the uc_end_ptrs
  for (i=0; i<evps.size(); i++) {
    if (evps[i].end) 
      evps[i].ptr->evpi_a = i; 
    else 
      evps[i].ptr->evpi_b = i; 
  }

  // keep count of the number of edges to go
  int edge_count = el.size(); 
  uc_polygon current_poly, empty_poly; 
  bool found; 

  if (trace) {
    cout << edge_count << " marked edges found\n";
  }

  while (edge_count > 0) {

    // skip to the first unused edge
    for (i=0; i < evps.size() && evps[i].ptr->mark_p==0; ++i);
    assert(i != evps.size());

    if (trace) 
      cout << "First edge:" << evps[i] << "\n"; 

    current_poly.push_back(evps[i].end_vertex());
    evps[i].ptr->mark_p = 0; // mark edge as used
    edge_count--; 

    j = i; 
    // seek forward adding vertices to current_poly
    while ((found = next_edge(evps, j)) && evps[j].ptr->mark_p!=0) {

      if (trace) cout << evps[j] << "\n"; 

      current_poly.push_back(evps[j].end_vertex());
      evps[j].ptr->mark_p = 0; // mark edge as used
      edge_count--; 
    }

    if (trace) {
      if (!found) cout << "There was no next edge\n"; 
      else cout << "Next edge (unmarked) was: " << evps[j] << "\n"; 
    }

    if (j!=i) { // it was a line segment so now seek back

      j = i; 

      if (trace) cout << "Was segment: working back to start\n";

      while ((found = prev_edge(evps, j)) && evps[j].ptr->mark_p!=0) {

	if (trace) cout << evps[j] << "\n"; 

	current_poly.push_front(evps[j].end_vertex());
	evps[j].ptr->mark_p = 0; // mark edge as used
	edge_count--; 
      }

      if (trace) {
	if (!found) cout << "There was no previous edge\n"; 
	else cout << "Previous edge (unmarked) was: " << evps[j] << "\n"; 

	cout << "Adding initial vertex: " << evps[j].the_point() << "\n"; 
      }
      
      // add the initial vertex of the segment
      current_poly.push_front(uc_vertex(evps[j].the_point(),0));
    }

    push_back(empty_poly);
    back().swap(current_poly); 
  }
}

