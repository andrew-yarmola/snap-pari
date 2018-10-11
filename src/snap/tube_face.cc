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
#include "tube_face.hh"
#include <algorithm>
#include <fstream>
#include "adaptive_curve.hh"
#include "picfile.hh"
#include "printable.hh"

using std::max;
using std::ofstream;

static bool show_errors=false;

double tube_face::epsilon = 5e-5;
double tube_face::boundary_epsilon = 5e-4; 

#define EPS tube_face::epsilon

tube_face::tube_face(eqsfunc const& G, Mark m, bool transpose)
  : surface(G), mark_num(m)
{
  R_matrix<3> Sph(0.0);
  Sph[0][0] = 1.0; 
  Sph[2][0] =-1.0; 
  Sph[0][2] =-1.0; 
  Sph[2][2] =-G.K_value()*G.K_value(); 
  bq_curve S(Sph, Mark(), G.position().imag); 
  
  edges.make_intervals(S);
  da_spec da(G.geodesic_distance(), G.position(), m);
  clipping_surfaces.push_back(da);
  edges.compute_bounds(box,G,rad);
}

void tube_face::set_to_partner(tube_face const& f, Complex const& oa)
{
  surface = eqsfunc(f.surface.geodesic_distance(), oa); 
  edges = f.edges; 
  edges.x_reflect(); 
  mark_num = f.mark_num.partner();
  edges.compute_bounds(box,surface,rad); 
}

int tube_face::clip(eqsfunc const& G, Mark m, bool ucover)
{
  // Check bounds to see if surfaces overlap at all. 
  if (G.geodesic_distance().real/2.0 > rad) return 0;
  if (!overlap(box, G.bounds(), !ucover)) return 0; 

  // bool ovl = overlap(box, G.bounds(), !ucover);

  bq_curve R(surface.restriction_matrix(G), m, G.position().imag); 

  int res = edges.clip(R);

  if (res) {
    da_spec da(G.geodesic_distance(), G.position(), m);
    clipping_surfaces.push_back(da);
    edges.compute_bounds(box, surface,rad);
  }

//   if (res && !ovl) {
//     cout << "Naughty overlap\n" << box << ' ' << mark_num << endl; 
//     cout << G.bounds() << ' ' << m << endl; 
//   }
  
  if (res == -1) {
    ofstream error_file("tube_errors");
    if (error_file) print_clippings(error_file); 
  }
  return res; 
}

#if 0
int tube_face::clip(tube_face const& f, bool transpose, bool ucover)
{
  return clip(f.surface, f.mark_num, transpose, ucover);
}
#endif

int tube_face::make_face(da_spec const& face, vector<da_spec> const& nbrs, bool report)
{
  // set up the basics
  *this = tube_face(eqsfunc(face.distance(), Zero), face.mark());

  // compute the edges
  da_spec nbr;
  int res=0; 
  int i, n = nbrs.size();
  for (i=0; i<n; i++) {
    nbr = nbrs[i]; 
    res = clip(eqsfunc(nbr.distance(), nbr.angle()), nbr.mark());
    if (res < 0) return res; 
  }

  // position correctly
  translate(face.angle()); 

  if (report) print_edge_markings();

  if (!validate(nbrs, report)) return 1; 

  // remove any components at infinity
  // clean_up(); 

  return 0; 
}

static Mark expected(vector<da_spec> const& spec, int j)
{
  Mark m;
  if (j > spec.size() || j < 0) return m;

  int i;
  for (i=0; i<j; i++) {
    if (spec[i] == spec[j]) break; 
  }
  return spec[i].mark();
}

void tube_face::sort_edges()
{
  edges.sort_intervals(EPS);
}

bool faces_match(tube_face const& a, tube_face const& b, double& urad, double eps, int report)
{
  interval_list a_for = a.edges; 
  interval_list b_rev = b.edges;
  b_rev.x_reflect(); 

  interval_difference(a_for, b_rev, eps, false); 

  if (a_for.L.size()==0 && b_rev.L.size()==0) return true; 

  list<eqs_interval>::iterator it;

  if (report) {
    cout << endl; 
    for (it=a_for.L.begin(); it!=a_for.L.end(); it++)
      cout << '+' << *it << endl; 
    for (it=b_rev.L.begin(); it!=b_rev.L.end(); it++)
      cout << '-' << *it << endl; 
  }

  int j; 
  Complex p; 
  double d; 
  urad = 1000.0;
  for (it=a_for.L.begin(); it!=a_for.L.end(); it++) {
    for (j=0; j<2; j++) {
      p = it->end(j); 
      d = a.surface(p.real, p.imag).distance(); 
      if (d < urad) urad = d; 
    }
  }
  for (it=b_rev.L.begin(); it!=b_rev.L.end(); it++) {
    for (j=0; j<2; j++) {
      p = it->end(j); 
      d = b.surface(p.real, p.imag).distance(); 
      if (d < urad) urad = d; 
    }
  }
  if (report) {
    cout << "Unmatched vertex radius: " << urad << endl; 
  }

  return false; 
}


ostream& operator << (ostream& out, tube_face const& f)
{
  f.edges.print_sorted(out, EPS); 
  return out; 
}

void tube_face::uc_print(ostream& out, vector<Complex> const& H) const
{
  double eps_sq = EPS*EPS; 
  Complex prev, e0, e1; 
  list<eqs_interval>::const_iterator it; 
  for (it = edges.L.begin(); it!=edges.L.end(); it++) {
    if (it!=edges.L.begin() && complex_modulus_squared(it->end(0) - prev) > eps_sq) {
      out << endl; 
    }
    out << "["; 
    e0 = it->end(0);
    surface(e0.real, e0.imag).uc_print(out, H);
    out << ", "; 
    e1 = it->end(1);
    surface(e1.real, e1.imag).uc_print(out, H);
    out << "] " << it->mark() << endl; 
    prev = e1; 
  }
}

void tube_face::print_edge_markings() const
{
  // cout << "Face: " << mark_num.end_index() << " has edges "; 
  edges.print_markings(); 
}

void tube_face::print_clippings(ostream& out) const
{
  int oldprec = out.precision(16);
  out << LPSeq(clipping_surfaces);
  out.precision(oldprec); 

  out << PF(edges); 
}

static bool close_complex(const Complex& a, const Complex& b, double eps)
{
  return complex_modulus_squared(b - a) < eps*eps; 
}

bool tube_face::has_matching(eqs_const_iter const& e, Mark const& mk, vector<Complex> const& H) const
{
  // mark for interval adjacent to e on face with marking mk. 
  Mark adjmk(mk.end_index(), -e->mark()); 

  eqs_const_iter i;
  for (i = edges.L.begin(); i!= edges.L.end(); ++i) {
    if (i->mark() == adjmk) return true; 
  }
  return false; 
}


class EqsIterator {
  list<eqs_interval>& edges;
  list<eqs_interval>::iterator k, path; 
  Complex prv; 
public: 
  EqsIterator(list<eqs_interval>& e) : edges(e) { first(); }
  EqsIterator(list<eqs_interval>& e, eqs_iter const& it) 
    : edges(e), k(it) { path = e.begin(); } 

  void first(); 

  EqsIterator& operator ++();
  EqsIterator  operator ++(int)
  { EqsIterator i(*this); ++(*this); return i; }

  bool path_done() const; 
  bool was_loop() const;

  void start_path() { path = k; }
  bool done() const { return k==edges.end(); }

  void erase_path(); 
  void set_mark(Mark m);

  list<eqs_interval>::iterator operator -> () { return k; }
  list<eqs_interval>::iterator it() { return k; }
};

// usage:
// EqsIterator i(intervals); 
// do something on each path
// while (!i.done()) {
//   for (i.start_path(); !i.path_done(); i++) {
//      do something on this path
//   }
// }
// 
// do something on each interval 
// for (i.first(); !i.done(); ++i) {
// }

void EqsIterator::first()
{
  k = edges.begin();
  path = k; 
  prv = Infinity;
}

EqsIterator& EqsIterator::operator ++()
{
  Mark m = k->mark(); 
  while (!done()) {
    prv = k->end(1); 
    ++k; 
    if (k->mark() != m || path_done()) break; 
  }
  return *this; 
}

bool EqsIterator::path_done() const
{
  if (done()) return true; 
  if (k==path) return false; 
  return !close_complex(prv, k->end(0), EPS);
}

bool EqsIterator::was_loop() const
{
  list<eqs_interval>::iterator l=k; --l; 
  return close_complex(l->end(1), path->end(0), EPS); 
}

void EqsIterator::erase_path()
{
  while (!path_done()) ++(*this); 
  edges.erase(path, k);
  path = k;
}

void EqsIterator::set_mark(Mark m)
{
  list<eqs_interval>::iterator j = k; 
  EqsIterator je = (*this); 
  ++je;
  while (j != je.k) {
    j->set_mark(m);
    ++j;
  }
}




// erase all paths having edges at infinity
void tube_face::clean_up()
{
  EqsIterator i(edges.L); 
  Mark at_infinity; 

  while (!i.done()) {
    if (i->mark() == at_infinity) 
      i.erase_path();
    else ++i;
    if (i.path_done()) i.start_path(); 
  }
}

static inline int pos_mod(int i, int n)
{
  if (i<0) i+=(-i/n+1)*n; 
  return i%n; 
}

typedef vector<eqs_iter> eqs_path;

class EqsPaths {
  list<eqs_interval>& edges;
  vector<eqs_path> paths;
  void compute_paths(); 
public: 

  EqsPaths(list<eqs_interval>& e) : edges(e) 
  { compute_paths(); }

  eqs_path& operator [] (int i) { return paths[i]; }

  eqs_iter it(int l, int p) const
  { int s = paths[l].size(); return paths[l][pos_mod(p,s)]; }
  Mark mark(int l, int p) const 
  { return it(l,p)->mark(); }
  void set_mark(int l, int p, Mark m);
  int path_size(int l) const { return paths[l].size(); }
  int size() const { return paths.size(); }
  void print(vector<da_spec> const& spec) const;
};

void EqsPaths::compute_paths()
{
  eqs_iter k;
  Complex prv = Infinity;
  EqsIterator i(edges);
  while (!i.done()) {
    paths.push_back(eqs_path()); 
    i.start_path();
    while (!i.path_done()) {
      paths.back().push_back(i.it()); 
      ++i;
    }
  }
}

void EqsPaths::set_mark(int l, int p, Mark m)
{
  EqsIterator i(edges, it(l,p));
  i.set_mark(m); 
}

class LoopIterator {
  EqsPaths& L; 
  int l; // loop
  int p; // position
  int d; // direction == +/- 1. 
public:

  LoopIterator(EqsPaths& ivl, int _l, int _p, int _d) 
    : L(ivl), l(_l), p(_p), d(_d) {}

  void operator = (LoopIterator const& i) 
  { l = i.l; p = i.p; d = i.d; }

  // value extraction
  
  Mark mark(int offset) const;
  Mark operator [] (int offset) const { return mark(offset); }
  void set_mark(Mark m, int offset=0); 
  int loop_size() const { return L[l].size(); }

  // iterator movements

  void operator += (int offset);
  bool pinch_off_partner(Mark& ne); 

  LoopIterator& operator ++() { (*this) += 1; return *this; }
  void operator ++ (int) { (*this) += 1; }
  friend bool operator == (LoopIterator const& a, LoopIterator const& b)
  { return a.l==b.l && a.p==b.p && a.d==b.d; }
  friend bool operator != (LoopIterator const& a, LoopIterator const& b)
  { return !(a==b); }
};

Mark LoopIterator::mark(int offset) const
{
  int s = L[l].size();
  return L.mark(l,pos_mod(p+d*offset,s));
}

void LoopIterator::set_mark(Mark m, int offset)
{
  int s = L[l].size();
  return L.set_mark(l,pos_mod(p+d*offset,s), m);
}

void LoopIterator::operator += (int offset)
{
  int s = L[l].size();
  p = pos_mod(p+d*offset,s);
}
  
bool LoopIterator::pinch_off_partner(Mark& ne)
{
  int i,j,s; 

  Mark e  = mark(0); 
  ne = mark(1); 

  // check each loop other than this one. 
  for (i=0; i<L.size(); i++) {
    if (i==l) continue; 

    // look for e. 
    s = L[i].size(); 
    for (j=0; j<s; j++) {
      if (L.mark(i,j)==e) break; 
    }
    if (j==s) continue; // didn't find it. 

    // check for ne. 
    // assuming common orientation of loops this first case 
    // should never occur. 
    if (L.mark(i,j+1) == ne) {
      // found it
      l = i;
      p = j;
      d = -1; 
      return true; 
    } 
    if (L.mark(i,j-1) == ne) {
      l = i; 
      p = j; 
      d = 1; 
      return true; 
    }      
  }
  return false; 
}

class ChangeReporter {
public:
  int face_index(int f, int e) const; 
  int tet_index(int f, int ea, int eb) const; 
};

enum change_type {
  edge_vanish,
  edge_appear,
  bite_vanish,
  bite_appear,
  vertex_open,
  pinch_off,
  bigon_vanish,
  triangle_vanish
};

const int change_num_args[] = {
  1,
  2,
  1,
  1,
  2,
  2,
  0,
  0
};

const char* change_name[] = {
  "edge vanished",
  "edge appeared",
  "bite vanished",
  "bite appeard",
  "vertex opened up",
  "pinch off occurred",
  "bigon vanished",
  "triangle vanished"
};


class AChange {
  change_type t;
  int face; 
  int ea, eb; 
  int num_args() const { return change_num_args[t]; }
public:
  AChange(change_type ct, int f, int a=-1, int b=-1) 
    : t(ct), face(f), ea(a), eb(b) {}
  friend ostream& operator << (ostream& out, AChange C);
};

ostream& operator << (ostream& out, AChange C)
{ 
  out << (C.face/2) << ' ' << change_name[C.t]; 
  if (C.num_args() > 0) out << ' ' << (C.ea+1);
  if (C.num_args() > 1) out << ' ' << (C.eb+1);
  return out; 
}



class intmod {
  int i, n;

  void normalize()
  { if (i<0) i+=(-i/n+1)*n; i = i%n; }

public:
  intmod(int j, int m) : i(j), n(m) { normalize(); }

  intmod& operator ++() { i = (i+1)%n; }

  intmod operator += (int b) { i+=b; normalize(); }

  operator int() { return i; }

  friend intmod operator + (intmod const& a, int b)
  { return intmod(a.i+b,a.n); }

  friend bool operator == (intmod const& a, intmod const& b)
  { return a.n==b.n && a.i==b.i; }
  friend bool operator != (intmod const& a, intmod const& b)
  { return !(a==b); }
  friend ostream& operator << (ostream& out, intmod const& a)
  { return out << "mod(" << a.i << ',' << a.n << ")"; }
};

class half_pinch {
  Mark x, y; 
public:
  int j;

  half_pinch() : j(-1) {}
  half_pinch(Mark _x, Mark _y, int _j) 
    : x(_x), y(_y), j(_j) {}

  friend bool match(half_pinch const& a, half_pinch const& b)
  { return a.y==b.x && a.x==b.y; }
};

static bool have_pinch(vector<half_pinch>& halves, half_pinch const& hp, int& j)
{
  int i, n=halves.size();
  for (i=0; i<n; i++)
    if (match(hp,halves[i])) break; 
  if (i==n) {
    halves.push_back(hp);
    return false; 
  }
  j = halves[i].j;
  for (; i<n-1; i++) 
    halves[i]=halves[i+1];
  halves.resize(n-1);
  return true;
}

// return values: 
// -1 - no match was possible.
// 0  - match without changes.
// nc - match with nc changes. 

static int try_match(LoopIterator const& i0, vector<da_spec> const& spec, intmod j0, int index, bool report)
{
  int n = spec.size();
  int nc = 0; // number of changes. 
  LoopIterator i = i0; 

  // check we've got a match at the starting point. 
  Mark ex = expected(spec,j0);
  if (i[0] != ex) return -1; 

  // we want to work our way around spec and find the
  // corresponding edges of L (as referenced by i). 

  intmod j=j0;
  int lim=0;

  // check for a reversed triangle
  if (i.loop_size()==3 && n==3 && 
      i[1]==expected(spec,2) &&
      i[2]==expected(spec,1)) return -1; 

  bool fix_markings=false; 
  Mark ex1, ex2, at_infinity; 
  int bite_limit = 16; 
  int pinch_limit = 8; 

  int jj;
  Mark imk; 
  vector<half_pinch> halves; 
  vector<AChange> changes; 

  while (lim < n) {

    ex  = expected(spec,j);
    ex1 = expected(spec,j+1);
    ex2 = expected(spec,j+2); 

    // the easy case, next edge matches as expected. 
    if (i[1] == ex1) {
      if (i[1] != spec[j+1].mark()) fix_markings=true; 
      ++i; ++j; ++lim; continue; 
    }

    // pinch off
    if (i.pinch_off_partner(imk) ) {
      if (--pinch_limit<0) {
	if (report) cout << "pinch limit reached\n"; 
	return -1; 
      }
      if (have_pinch(halves,half_pinch(i[0],imk,j),jj))
	changes.push_back(AChange(pinch_off,index,j,jj));
      ++nc; continue;
    }

    // a vertex of the face has opened up
    if (i[1] == at_infinity && i[2] == ex1) {
      changes.push_back(AChange(vertex_open,index,j,j+1));
      i+=2; ++j; ++lim; ++nc; continue; 
    }

    // edge vanish
    if (i[1] == ex2) {
      changes.push_back(AChange(edge_vanish,index,j+1));
      ++i; j+=2; lim+=2; ++nc; continue; 
    }

    // edge has appeared
    if (i[2] == ex1) {
      changes.push_back(AChange(edge_appear,index,j,j+1));
      i+=2; ++j; ++lim; ++nc; continue; 
    }

    // a bite has been put back
    if (ex2 == ex && n > 2) {
      changes.push_back(AChange(bite_vanish,index,j+1));
      j+=2; lim+=2; ++nc; continue; 
    }

    // someone has taken a bite
    if (i[2] == ex) {
      if (--bite_limit<0) {
	if (report) cout << "bite limit reached!\n"; 
	return -1; 
      }
      changes.push_back(AChange(bite_appear,index,j));
      i+=2; ++nc; continue; 
    }

    return -1; 
  }

  if (j!=j0) return -1; 

  if (i!=i0) {
    // check again for bite or pinch-off
    if (i[2] == ex) {
      i+=2; ++nc; 
    } else if (i.pinch_off_partner(imk) ) {
      if (have_pinch(halves,half_pinch(i[0],imk,j),jj))
	changes.push_back(AChange(pinch_off,index,j,jj));
      ++nc;
    }
  }
  if (i!=i0) return -1; 

  int k;
  if (report) 
    for (k=0; k<changes.size(); k++)
      cout << changes[k] << endl; 

  if (nc==0 && fix_markings) {
    lim = 0;
    for (lim=0; lim < n; ++lim) {
      if (i[0] != spec[j].mark())
	i.set_mark(spec[j].mark());
      ++i; ++j;
    }
  }

  return nc; 
}

// return codes as for try_match. 

static int seek_match(EqsPaths& L, vector<da_spec> const& spec, int index, bool report)
{
  int n = spec.size();

  // Single face vanishing?
  if (L.size()==0) {
    if (n>3) {
      if (report) cout << "NON-GENERIC: vanishing of " << n << "-face\n";
      return -1;
    } 
    if (report) cout << "edge " << (index/2) << " of order " << n << " vanishes\n"; 
    return 1; 
  }

  int dir = (index%2) ? 1 : -1; 

  // Now look for a correspondence. 
  intmod j(0,n);
  int k; 
  int l, s, p;
  int res; 
  for (k=0; k<n; ++j, k++) {
    for (l=0; l<L.size(); l++) {
      s = L[l].size(); 
      for (p=0; p<s; p++) {

	if (L.mark(l,p) != expected(spec,j)) continue; 
	if ((res=try_match(LoopIterator(L,l,p,dir),spec,j,index,report)) > -1) 
	  return res; 
      }
    }
  }

  return -1; 
}

bool tube_face::validate(vector<da_spec> const& nbrs, bool report)
{
  EqsPaths L(edges.L); 
  // if (report) L.print(nbrs); 
  return seek_match(L, nbrs, mark_num.end_index(), report)==0; 
}

static string mark_index(vector<da_spec> const& spec, Mark m)
{
  int i, n=spec.size();
  char buf[20];
  string out; 
  bool first_write=true; 
  for (i=0; i<n; i++) {
    if (expected(spec,i) != m) 
      continue; 
    if (!first_write) out += "/";
    sprintf(buf, "%d", i+1);
    out += buf; 
    first_write = false; 
  }
  if (first_write) return "*";
  return out;
}

void EqsPaths::print(vector<da_spec> const& spec) const
{
  int l, p, n=size(), s; 
  
  cout << '(';
  for (l=0; l<n; l++) {
    if (l) cout << ' ';

    cout << '(';
    s = paths[l].size(); 
    for (p=0; p<s; p++) {
      if (p) cout << ',';
      cout << mark_index(spec,mark(l,p));
    }
    cout << ')';
  }
  cout << ')' << endl;
}



#if 0
void tube_face::clean_up() 
{
  double eps_sq = EPS*EPS; 
  Complex prev, e0, e1; 
  Mark at_infinity; 
  bool to_erase = false; 
  list<eqs_interval>::iterator it, i0; 

  for (i0 = it = edges.L.begin(); it!=edges.L.end(); it++) {
    if (it!=edges.L.begin() && complex_modulus_squared(it->end(0) - prev) > eps_sq) {
      if (to_erase) edges.L.erase(i0,it);
      i0 = it; 
      to_erase = false; 
    }
    if (it->mark()==at_infinity) to_erase = true; 
    prev = it->end(1); 

  }
  if (to_erase) edges.L.erase(i0,it);
}
#endif


void tube_face::translate(Complex const& z)
{
  surface = surface + z; 
}

#if 0
/* This is the blueprint for traversing the closed loops in 
   a face whose edges have been sorted. */ 

void tube_face::traverse_loops() const
{
  if (!edges.L.size()) return; 

  list<eqs_interval>::const_iterator k, loop;
  Complex vx; 
  bool was_loop;

  k=edges.L.begin();
  while (k!=edges.L.end()) {
    loop = k; 

    // Traverse one loop or curve. 
    while (true) {

      vx = k->end(1); k++;
      if (k==edges.L.end() || !close_complex(vx, k->end(0), EPS)) break; 

      // This interval and the previous one belong to the same curve. 
    }
    was_loop = close_complex(vx, loop->end(0), EPS);

  }
}
#endif 

static bool first_non_short(eqs_const_iter& k, list<eqs_interval> const& L)
{
  Complex z, z0 = k->end(0);
  
  while (true) {
    if (!close_complex(z0, (z=k->end(1)), tube_face::boundary_epsilon)) 
      break; 
    ++k; 
    if (k == L.end() || !close_complex(z, (z0 = k->end(0)), EPS)) 
      return false; 
  }
  return true; 
}

static bool next_non_short(eqs_const_iter& k, list<eqs_interval> const& L, Complex& z)
{
  Complex z0;
  z = k->end(1);

  while (true) {
    ++k; 
    if (k == L.end() || !close_complex(z, (z0 = k->end(0)), EPS)) 
      return false; 
    if (!close_complex(z0, (z = k->end(1)), tube_face::boundary_epsilon)) 
      break; 
  }
  return true; 
}


uc_polylist tube_face::get_boundary() const
{
  uc_polylist boundary; 
  if (!edges.L.size()) {
    return boundary; 
  }

  list<eqs_interval>::const_iterator k, loop;
  Complex z, prev_z, last_on_cpt; 

  uc_polygon b_cpt, empty; 
  Mark prev_mark_num; 

  k=edges.L.begin();
  while (k!=edges.L.end()) {
    loop = k; 

    if (!first_non_short(k,edges.L)) 
      continue; // all edges short on this component

    prev_mark_num = k->mark(); 
    prev_z = k->end(1); 

    // Get a closed loop if possible. 
    while (true) {

      if (!next_non_short(k,edges.L,last_on_cpt)) break; 

      z = 0.5 * (prev_z + k->end(0)); 
      b_cpt.push_back(uc_vertex(surface(z.real,z.imag),
				prev_mark_num.value(), k->mark().value()));

      prev_mark_num = k->mark(); 
      prev_z = k->end(1); 
    }

    if (close_complex(last_on_cpt, loop->end(0), EPS)) { // Was loop

      first_non_short(loop, edges.L);

      z = 0.5 * (prev_z + loop->end(0)); 
      b_cpt.push_back(uc_vertex(surface(z.real, z.imag), 
				prev_mark_num.value(), loop->mark().value()));
    }

    // Save the loop. 
    boundary.push_back(empty); 
    boundary.back().splice(boundary.back().end(), b_cpt); 

  }
  return boundary;
}


void tube_face::picfile_print(picfile& pic, color const& col) const
{
  list<Complex> zero_trans;
  zero_trans.push_back(Zero); 
  picfile_print(pic, col, zero_trans); 
}

void tube_face::picfile_print(picfile& pic, color const& col, list<Complex> const& trans) const
{
  if (!edges.L.size()) return; 

  list<Complex> segment, curve, this_curve, empty; 
  list<eqs_interval>::const_iterator k;

  list<Complex>::const_iterator it, i;
  Complex loop, vx; 
  bool was_loop; 
  k=edges.L.begin();
  while (k!=edges.L.end()) {
    loop = k->end(0);

    // Get a closed loop if possible. 
    while (true) {
      segment = adaptive_curve(linkpic_func(*k, surface), 
			       k->x(0), k->x(1), .2, 0.05);
      curve.splice(curve.end(), segment); 
      vx = k->end(1); 
      k++;
      if (k==edges.L.end() || !close_complex(vx, k->end(0), EPS)) break; 
    }
    was_loop = close_complex(vx, loop, EPS);

    // Print it out. 
    for (i=trans.begin(); i!=trans.end(); i++) {
      // Translate the curve. 
      for (it = curve.begin(); it!=curve.end(); it++) 
	this_curve.push_back(*it + *i); 
      
      if (was_loop && col != black)
	pic.print_polygon(this_curve, col, true); // true = outline
      else 
	pic.print_line(this_curve); 
      this_curve = empty; 
    }
    curve = empty; 
  }
}

void tube_face::picfile_print(picfile& pic, Complex const& center, double size) const
{
  if (!edges.L.size()) return; 

  // Get the curves of the intervals. 
  list<Complex> all_curves, segment, curve; 
  list<eqs_interval>::const_iterator k;
  for (k=edges.L.begin(); k!=edges.L.end(); k++) {
    segment = k->get_polyline(); 
    all_curves.splice(all_curves.end(), segment); 
  }

  // Find a bounding box for the curves. 
  bbox box; 
  box.include(all_curves); 
  double extent = 0.0001; // Avoid division by zero if somehow box is zero size. 
  extent = max(max(box.top(), -box.bottom()), max(-box.left(), box.right()));
  double scale = 0.45 * size/extent;

  list<Complex>::iterator it;
  Complex loop, vx; 
  k=edges.L.begin();
  while (k!=edges.L.end()) {
    loop = k->end(0);
    while (true) {
      segment = k->get_polyline(); 
      curve.splice(curve.end(), segment); 

      vx = k->end(1); 
      k++;
      if (k==edges.L.end() || !close_complex(vx, k->end(0), EPS)) break; 
    }
    // Rescale and reposition curve. 
    for (it = curve.begin(); it!=curve.end(); it++) 
      *it = scale * *it + center; 

    if (close_complex(vx, loop, EPS))
      pic.print_polygon(curve, color(1.0, 0.85, 0.85), true); // true = outline
    else 
      pic.print_line(curve); 
    curve = segment; // segment is empty. 
  }
}

