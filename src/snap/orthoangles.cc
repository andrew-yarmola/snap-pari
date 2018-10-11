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
#include "kernel_extras.hh"
#include "orthoangles.hh"
#include "snappea/kernel.h"
#include "complex_lll.hh"
#include "printable.hh"

using std::cout;
using std::cerr;
using std::endl;

static inline Complex trace(MoebiusTransformation const& m)
{
  return m.matrix[0][0] + m.matrix[1][1]; 
}

// Requires the unsimplified fundamental group. 
// returns cosh orthodistances. 

bool get_edge_orthodistances(Triangulation* manifold, 
			     GroupPresentation* G, 
			     vector<Complex>& orth)
{
  int i, n = get_num_tetrahedra(manifold); 
  orth.resize(n);

  if (get_num_cusps(manifold)!=1) {
    cerr << "Sorry, get_edge_orthodistances not implemented yet for multi-cusped triangulations\n";
    return false; 
  }

  vector<FGWord> core_ort_words; 
  get_end_pairing_words(manifold, core_ort_words); 
  if (core_ort_words.size()!=n) {
    cerr << "unexpected error in get_edge_orthodistances\n"; 
    return false; 
  }

  FGWord lw = fg_longitude(G, 0);
  MoebiusTransformation l = word_to_Moebius(G, lw); 

  // We use the following formula for cosh orthodistance between 
  // the axis of l and the axis of e.l.e^(-1). 
  // 
  // (cosh)orthodist = 2*(tr(le) tr(l^(-1)*e) - tr^2(e))/(tr^2(l) - 4) - 1. 

  Complex trl = trace(l);
  Complex L = trl*trl - Four; 
  if (L == Zero) {
    for (i=0; i<n; i++) orth[i] = Infinity; 
    return true; // All infinite. 
  }

  Complex tre;
  MoebiusTransformation e, li = inverse(l); 
  for (i=0; i<n; i++) {
    e = word_to_Moebius(G, fg_word_from_original(G, core_ort_words[i]));
    tre = trace(e); 
    orth[i] = Two * (trace(l*e) * trace(li*e) - tre * tre)/L - One; 
  }
  return true; 
}

// The formula sqrt(x^2-1) returning sinh(w) when x=cosh(w), requires
// a choice of sign. We choose the sign such that w has positive real
// part, i.e. such that exp(w) (== sinh(w) + cosh(w)) has absolute value >=1. 

static Complex pd_sinh(Complex const& csh)
{
  Complex snh = complex_sqrt(csh*csh - One);
  return (complex_modulus_squared(csh + snh) < 1.) ? -snh : snh;
}

// This is just the cosine rule in 3-dimensional hyperbolic geometry. 

static Complex orth_angle(Complex const& a, Complex const& b, Complex const& c)
{
  return complex_acosh((b * c - a)/(pd_sinh(b)*pd_sinh(c)));
}

// We supply the orthodistances visible in a certain 
// tetrahedron, seen from the point of view of a particular
// vertex v0 of one of its ideal vertex triangles T. We suppose
// that the orthoangle of v0 is zero, that we already know
// the orthoangle of v1 and wish to determine the orthoangle
// of v2. We have to give all 6 orthodistances of the hextet
// in the order d0 d1 d2 d01 d02 d12. These determine the 
// orthoangles up to sign, so prev_oa is only required in 
// order to choose the right sign. 

static Complex next_orthoangle(Complex d[6], Complex prev_oa = Zero)
{
  Complex oa_01 = orth_angle(d[3],d[0],d[1]);

  // sort out the sign of oa_01
  if (prev_oa != Zero) {
    if (prev_oa != oa_01) {
      oa_01 = -oa_01;
      if (prev_oa != oa_01) {
	cerr << "next_orthoangle can't match previous angle\n";
      }
    }
  }

  Complex oa_02 = orth_angle(d[4],d[0],d[2]);
  Complex oa_12 = orth_angle(d[5],d[1],d[2]); 

  int i, j;
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      if (oa_02 == oa_01 + oa_12) break; // success!
      oa_12 = -oa_12; // try flipping oa_12
    }
    if (j<2) break; // success!
    oa_02 = -oa_02; 
  }
  if (i==2) {
    cerr << "next_orthoangle: orthoangles don't form a triangle!\n";
  }

  // cout << oa_01 << ',' << oa_02 << ',' << oa_12 << endl; 

  return oa_02;
}

static void get_the_orthodistances(vector<Complex> const& orth, EdgeIterator const& i,
				   Complex d[6])
{
  d[0] = orth[i.tet->edge_class[edge_between_vertices[i.V[0]][i.V[1]]]->index];
  d[1] = orth[i.tet->edge_class[edge_between_vertices[i.V[0]][i.V[3]]]->index];
  d[2] = orth[i.tet->edge_class[edge_between_vertices[i.V[0]][i.V[2]]]->index];
  d[3] = orth[i.tet->edge_class[edge_between_vertices[i.V[1]][i.V[3]]]->index];
  d[4] = orth[i.tet->edge_class[edge_between_vertices[i.V[1]][i.V[2]]]->index];
  d[5] = orth[i.tet->edge_class[edge_between_vertices[i.V[3]][i.V[2]]]->index];
}

void get_edge_da_spec(vector<Complex> const& orth, EdgeClass* e, vector<Complex>& da_spec, bool report)
{
  da_spec.resize(0); 
  
  Complex d[6], oa = Zero, first_oa, od; 

  // Now just work around the edge printing them out. 
  EdgeIterator it(e); 
  int ord = e->order; 

  // Work out the first pair. 
  get_the_orthodistances(orth, it, d);
  oa = orth_angle(d[3],d[0],d[1]);
  od = complex_acosh(d[1]);
  first_oa = oa; 

  for (; it.index<ord; ++it) {

    if (od.real < 0.) od = -od; 
    if (report) cout << "[" << oa << "," << od << "]" << endl; 

    da_spec.push_back(od); 
    da_spec.push_back(oa); 

    get_the_orthodistances(orth, it, d);
    oa = next_orthoangle(d, oa); 
    od = complex_acosh(d[2]); 
  }
  
  if (oa != first_oa) {
    cerr << "we didn't get back to our starting ortho-angle\n"; 
  }
}

void print_oa_cycles(Triangulation* manifold, GroupPresentation* G)
{
  // Get all the edge orthodistances
  vector<Complex> orth;
  if (!get_edge_orthodistances(manifold, G, orth)) {
    return; 
  }
  if (orth[0]==Infinity) {
    cerr << "Sorry, can't compute orthoangles for complete structure.\n";
    return;
  }

  OrthoAngles OA(manifold, G);
  if (!OA.valid()) return;

  OA.print(); 

  EdgeClass *e; 
  vector<Complex> da_spec; 

  int j;
  vector<Complex> da(2); 

  // Run through the edges
  for (e = manifold->edge_list_begin.next; 
       e!= &manifold->edge_list_end;
       e = e->next) {

    cout << "Edge: " << e->index << endl; 
    OA.get_da_spec(e, da_spec);
    for (j=0; j<da_spec.size(); j+=2) {
      da[1] = da_spec[j]; 
      da[0] = da_spec[j+1]; 
      cout << PSeq(da) << endl;  
    }
    cout << endl; 
    get_edge_da_spec(orth, e, da_spec, true); // true = print it out
  }
}

static const FaceIndex other_face_at_edge3[4][3] =
  {{1,2,3},
   {0,3,2},
   {3,0,1},
   {2,1,0}};

class TriangleIterator {
  Tetrahedron* tet; 
  VertexIndex v;
public:
  TriangleIterator() : tet(0), v(0) {}
  TriangleIterator(Tetrahedron* t, VertexIndex vx) : tet(t), v(vx) {}

  int f3(int f) const { return edge3_between_faces[v][f]; }
  int f4(int f) const { return other_face_at_edge3[v][f]; }
  int index() const 
  { return tet->index*4 + v; }
  int vertex_edge_index(VertexIndex vx) const
  { return tet->edge_class[edge_between_vertices[v][f4(vx)]]->index; }
  int face_edge_index(FaceIndex f) const
  { return tet->edge_class[edge_between_faces[v][f4(f)]]->index; }

  TriangleIterator neighbor(FaceIndex f, FaceIndex& nf) const; 
  void get_gluing(FaceIndex f, VertexIndex v[3]) const; 

  void first(Triangulation* T)
  { tet = T->tet_list_begin.next; v = 0; }
  bool done(Triangulation* T) const
  { return tet == &T->tet_list_end; }

  TriangleIterator& operator ++(); 
  TriangleIterator  operator ++(int)
  { TriangleIterator ti(*this); ++(*this); return ti; }

  // what |= 0x1 for neighbors, 0x2 for vertices. 
  void print(ostream& out, int what=0) const; 

  EdgeEnd edge_end(VertexIndex v) const; 

  // friend class EdgeEnd;
  friend bool operator ==(TriangleIterator const& a, TriangleIterator const& b)
  { return a.tet==b.tet && a.v == b.v; }
  friend bool operator !=(TriangleIterator const& a, TriangleIterator const& b)
  { return !(a==b); }
};

TriangleIterator TriangleIterator::neighbor(FaceIndex f, FaceIndex& nf) const
{
  FaceIndex F = f4(f); // adjust incoming face index (0-2) -> (0-3)
  Permutation gluing = tet->gluing[F]; 
  TriangleIterator nt(tet->neighbor[F], EVALUATE(gluing, v));
  nf = nt.f3(EVALUATE(gluing, F));
  return nt; 
}

void TriangleIterator::get_gluing(FaceIndex f, VertexIndex V[3]) const
{
  FaceIndex F = f4(f); 
  Permutation gluing = tet->gluing[F]; 
  TriangleIterator nt(tet->neighbor[F], EVALUATE(gluing, v));

  VertexIndex u;
  for (u=0; u<3; u++) {
    V[u] = nt.f3(EVALUATE(gluing, f4(u)));
  }
}

TriangleIterator& TriangleIterator::operator ++()
{
  ++v;
  if (v>3) {
    v = 0; 
    tet = tet->next; 
  }
  return *this; 
}

EdgeEnd TriangleIterator::edge_end(VertexIndex vx) const
{
  EdgeIndex e = edge_between_vertices[f4(vx)][v];
  return EdgeEnd(tet->edge_class[e], edge_base(tet,e)==v);
}

int EdgeEnd::order() const
{
  return edge->order;
}

int EdgeEnd::index() const
{
  return 2*edge->index + (base_end ? 1:0); 
}

EdgeEnd& EdgeEnd::operator ++()
{
  if (base_end) edge = edge->next; 
  base_end = !base_end;
}

void EdgeEnd::first(const Triangulation* T)
{
  edge = T->edge_list_begin.next; base_end = false; 
}

bool EdgeEnd::done(const Triangulation* T) const
{
  return edge == &T->edge_list_end; 
}

EdgeIterator EdgeEnd::edge_iterator() const
{
  return EdgeIterator(edge, base_end ? 1:0);
}

ostream& operator << (ostream& out, EdgeEnd const& e)
{
  return out << 'E' << e.edge->index << ',' << (e.base_end ? 'b':'t');
}



// print as: (tet,vertex) 
// or: (t,v) (t0,v0):v00,v01,v02 (t1,v1):v10,v11,v12 (t2,v2):v20,v21,v22

void TriangleIterator::print(ostream& out, int print_what) const
{
  out << '(' << tet->index << ',' << int(v) << ')';

  const int print_nbrs = 0x1; 
  const int print_vertices = 0x2; 

  if (print_what & print_nbrs) { 
    cout << ' ';
    FaceIndex f, nf;
    int j; 
    VertexIndex perm[3];
    TriangleIterator nbr;
    for (f=0; f<3; f++) {
      nbr = neighbor(f, nf);
      nbr.print(out, false);
      cout << ':';
      get_gluing(f, perm);
      for (j=0; j<3; j++) {
	cout << int(perm[j]);
	if (j<2) cout << ',';
      }
      if (f<2) cout << ' ';
    }
  }

  if (print_what & print_vertices) { 
    cout << ' '; 
    int i; 
    for (i=0; i<3; i++) {
      cout << edge_end(i);
      if (i<2) cout << ' ';
    }
  }
}

TriangleIterator get_triangle(EdgeIterator const& ei, VertexIndex& v)
{
  TriangleIterator ti(ei.tet, ei.V[1]);
  v = ti.f3(ei.V[0]);
  return ti; 
}

class OrthoAngleSetup {
  Triangulation* t;
  vector<Complex> _rel_orthoangle; // per triangle-face
  vector<Complex> _abs_orthoangle; // per edge-end
  vector<Complex> D;   // orthodistances, per edge
  vector<Complex> holonomies; 

  Complex rel_oa(TriangleIterator const& it, FaceIndex f) const
  { return _rel_orthoangle[3*it.index() + f]; }
  void set_rel_oa(TriangleIterator const& it, FaceIndex f, Complex const& z)
  { _rel_orthoangle[3*it.index() + f] = z; }
  bool triangle_done(TriangleIterator const& it) const; 
  bool do_triangle(TriangleIterator const& it, FaceIndex f);
  bool find_undone_triangle(TriangleIterator& it) const; 

  bool unassigned(EdgeEnd const& ee) const
  { return _abs_orthoangle[ee.index()] == Infinity; }
  Complex abs_oa(EdgeEnd const& ee) const 
  { return _abs_orthoangle[ee.index()]; }
  bool set_abs_oa(EdgeEnd const& , Complex const& z);
  Complex orthodistance(EdgeEnd const& ee) const;
  bool add_holonomy(Complex const& z); 
  bool find_unassigned(EdgeEnd& ee) const; 

  bool compute_relative_orthoangles(); 
  bool compute_absolute_orthoangles(); 
  void reduce_orthoangles(); 

  bool check_f3f4() const; 
  void print_triangle(TriangleIterator const& ti) const; 
public:
  OrthoAngleSetup(Triangulation* T, vector<Complex> const& orth);

  void print() const; 
  void get_holonomies(vector<Complex>& H) const; 
  void get_da_spec(EdgeEnd const& ee, vector<Complex>& da_spec) const; 
  void get_da_spec(EdgeIterator const& ei, da_spec& da, bool report=false) const;
  void get_da_spec(EdgeEnd const& ee, da_spec& da) const;
  void get_da_spec(EdgeEnd const& ee, vector<da_spec>& spec) const; 
  bool valid() const { return holonomies.size()==2; } 
};

OrthoAngleSetup::OrthoAngleSetup(Triangulation* T, vector<Complex> const& orth)
  : t(T), D(orth), 
    _rel_orthoangle(get_num_tetrahedra(T)*12, Infinity), 
    _abs_orthoangle(get_num_tetrahedra(T)*2,  Infinity)
{
  if (D.size()!=get_num_tetrahedra(T)) {
    cerr << "Incorrect size vector of orthodistances\n"; 
    return; 
  }
  if (!compute_relative_orthoangles()) {
    cerr << "Trouble finding relative orthoangles\n";
    return; 
  }
  if (!compute_absolute_orthoangles()) {
    cerr << "Trouble finding absolute orthoangles\n";
    return; 
  }
  reduce_orthoangles(); 
}

bool OrthoAngleSetup::check_f3f4() const
{
  TriangleIterator i; 
  int v, f; 
  for (v=0; v<4; v++) {
    for (f=0; f<3; f++) {
      i = TriangleIterator(0,v); 
      if (i.f3(i.f4(f)) != f) return false; 
    }
  }
  return true; 
}

void OrthoAngleSetup::print_triangle(TriangleIterator const& ti) const
{
  const int print_nbrs_and_vertices = 0x3; 
  ti.print(cout, print_nbrs_and_vertices); 
  cout << endl; 
  int t_side = 3*ti.index(), j, i=0;
  vector<Complex> side(3);
  for (j=t_side; j<t_side+3; ++j) 
    side[i++] = _rel_orthoangle[j];
  cout << PSeq(side) << endl; 
}

void OrthoAngleSetup::print() const
{
  if (!check_f3f4()) cout << "f3f4 broken\n";

  cout << "Relative orthoangles:\n";
  TriangleIterator ti; 
  for (ti.first(t); !ti.done(t); ++ti)
    print_triangle(ti); 

  cout << "\nAbsolute orthoangles:\n" << PSeq(_abs_orthoangle); 
  cout << "\nHolonomies:\n" << PSeq(holonomies) << endl;
}

// A triangle is done when all its sides have a relative orthoangle. 

bool OrthoAngleSetup::triangle_done(TriangleIterator const& it) const
{
  FaceIndex f;
  for (f=0; f<3; f++) {
    if (rel_oa(it, f)==Infinity) return false; 
  }
  return true; 
}

// do_triangle computes the relative orthoangles on this triangle, 
// keeping the sign of the orthoangle on side f fixed. 

bool OrthoAngleSetup::do_triangle(TriangleIterator const& it, FaceIndex f)
{
  FaceIndex i,j;
  Complex z[3]; 
  for (i=0; i<3; i++) {
    z[i] = orth_angle(D[it.face_edge_index(i)], 
		      D[it.vertex_edge_index((i+1)%3)], 
		      D[it.vertex_edge_index((i+2)%3)]);
  }

  // if side f is already set it should already equal +/- z[f]. 
  Complex zf = rel_oa(it,f); 
  if (zf == Infinity) {
    zf = z[f];
  } else if (z[f] != zf) {
    if (z[f] != -zf) return false;
    z[f] = zf;
  }

  int f1 = (f+1)%3; 
  int f2 = (f+2)%3; 

  // try to make sum zero by flipping signs of z[f1] and z[f2]. 
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      if (zf + z[f1] + z[f2] == Zero) break; 
      z[f1] = -z[f1]; 
    }
    if (j < 2) break; 
    z[f2] = -z[f2];
  }
  if (i==2) return false; // couldn't make the sum zero

  // set my relative orthoangles
  for (i=0; i<3; i++) 
    set_rel_oa(it, i, z[i]); 

  // set neighbor relative orthoangles
  FaceIndex nf; 
  TriangleIterator nt;
  Complex nz;
  for (i=0; i<3; i++) {
    nt = it.neighbor(i, nf);
    nz = rel_oa(nt, nf);
    if (nz != Infinity) {
      if (nz != -z[i]) return false;
    } else {
      set_rel_oa(nt, nf, -z[i]); 
    }
  }
  
  return true; 
}

bool OrthoAngleSetup::find_undone_triangle(TriangleIterator& it) const
{
  Tetrahedron* tet; 
  VertexIndex v; 
  for (it.first(t); !it.done(t); ++it) 
    if (!triangle_done(it)) return true; 

  return false; 
}

bool OrthoAngleSetup::compute_relative_orthoangles()
{
  TriangleIterator it, nt, nt2; 
  vector<TriangleIterator> todo; // really means nbrs todo.

  int current; 
  FaceIndex f, nf, nf2; 
  while (find_undone_triangle(it)) {

    // visit first triangle and put it on queue
    if (!do_triangle(it,0)) return false; 
    todo.resize(0);
    todo.push_back(it); 
    current=0; 

    while (current < todo.size()) {
      
      // pull triangle off the queue
      it = todo[current++]; 

      // visit its neighbors
      for (f=0; f<3; f++) {
	nt = it.neighbor(f,nf); 

	nt2= nt.neighbor(nf,nf2);
	if (nt2 != it || nf2 != f) {
	  cout << "I am not my neighbors neighbor!\n";
	}

	if (triangle_done(nt)) continue;
	if (!do_triangle(nt,nf)) return false; 
	todo.push_back(nt);
      }
    }
  }

  return true; 
}

inline static int round_to_int(double x)
{
  return (int)floor(x + .5); 
}

inline static double frac_part(double x)
{
  return x - floor(x + .5); 
}

inline static bool small(double x, double eps)
{
  return fabs(x) < eps; 
}

bool torus_reduce(vector<Complex> const& b, Complex& z)
{
  double C[2];
  if (!get_dep(b, z, C)) return false; // b is lin dep
  int i;
  z = Zero; 
  for (i=0; i<2; i++) {
    C[i] -= floor(C[i] + .5 - 1e-8);
    z += C[i] * b[i]; 
  }
  return true; 
}

bool torus_reduce(const Complex b[], Complex& z)
{
  vector<Complex> B(2);
  B[0]=b[0]; B[1]=b[1];
  return torus_reduce(B,z);
}

static bool get_int_dep(vector<Complex> const& b, Complex z, vector<int>& c)
{
  const float eps = 1e-8; 

  if (b.size() > 2) return false; // expect a basis

  c.resize(b.size());

  // 0-dimensional case.
  if (b.size() ==0) {
    if (complex_small(z,eps)) return true; 
    return false; 
  }

  // 1-dimensional case.
  if (b.size() ==1) {
    Complex q = z/b[0]; 
    if (!(small(q.imag, eps) && small(frac_part(q.real), eps))) 
      return false;
    c[0] = round_to_int(q.real); 
    return true; 
  }

  // 2-dimensional case. 
  double C[2];
  if (!get_dep(b, z, C)) return false; 
  int i;
  for (i=0; i<2; i++) {
    if (!small(frac_part(C[i]),eps)) return false; 
    c[i] = round_to_int(C[i]); 
  }
  return true; 
}

bool OrthoAngleSetup::add_holonomy(Complex const& z)
{
  // update lattice generators. 
  holonomies.push_back(z);
  lll_reduce(holonomies);
  remove_zeros(holonomies);

  // check the lattice is still discrete. 
  vector<int> c; 
  get_int_dep(holonomies, z, c); 

  // if it's not discrete c will have a big entry.
  int i, n = c.size(), ac, cmax=0; 
  for (i=0; i<n; i++) {
    ac = (c[i]>0) ? c[i] : -c[i]; 
    if (ac > cmax) cmax = ac; 
  }

  return cmax < 10; 
}

bool OrthoAngleSetup::set_abs_oa(EdgeEnd const& ee, Complex const& z)
{
  int i = ee.index();
  if (_abs_orthoangle[i]==Infinity) {
    _abs_orthoangle[i] = z; 
    return true; 
  } 
  return add_holonomy(z - _abs_orthoangle[i]); 
}

bool OrthoAngleSetup::find_unassigned(EdgeEnd& ee) const
{
  for (ee.first(t); !ee.done(t); ++ee)
    if (unassigned(ee)) return true; 
  return false; 
}

bool OrthoAngleSetup::compute_absolute_orthoangles()
{
  EdgeEnd ee, nbr;
  vector<EdgeEnd> todo; 
  int current; 
  EdgeIterator it; 
  TriangleIterator ti;
  VertexIndex v; 
  Complex nbr_z; 
  bool queue_nbr; 

  while (find_unassigned(ee)) {

    // visit first end and put it on the queue
    set_abs_oa(ee, Zero); 
    todo.resize(0);
    todo.push_back(ee); 
    current = 0; 

    while (current < todo.size()) {

      // pull end off the queue
      ee = todo[current++];  

      // visit its neighbors
      for (it = ee.edge_iterator(); it.index < ee.order(); it++) {
	ti = get_triangle(it,v); 
	nbr = ti.edge_end((v+1)%3); 
	nbr_z = abs_oa(ee) + rel_oa(ti, (v+2)%3); 
	queue_nbr = unassigned(nbr); 
	if (!set_abs_oa(nbr, nbr_z)) return false; 
	if (queue_nbr) todo.push_back(nbr); 
      }
    }
  }
  return true; 
}

void OrthoAngleSetup::get_holonomies(vector<Complex>& H) const
{
  if (holonomies.size()!=2) {
    cerr << "Invalid holonomies\n";
    return; 
  }
  H = holonomies;
}

Complex OrthoAngleSetup::orthodistance(EdgeEnd const& ee) const
{
  Complex d = complex_acosh(D[ee.index()/2]);
  if (d.real < 0) d = -d; 
  return d; 
}

void OrthoAngleSetup::get_da_spec(EdgeIterator const& ei, da_spec& da, bool report) const
{
  VertexIndex v; 
  TriangleIterator ti = get_triangle(ei,v); 
  EdgeEnd ee = ti.edge_end(v);
  EdgeEnd e2 = ti.edge_end((v+1)%3);
  Complex d, a;
  d = orthodistance(e2); 
  a = rel_oa(ti, (v+2)%3);

  vector<int> c(2,0); 
  if (holonomies.size()==2) {
    Complex hol = -abs_oa(e2) + abs_oa(ee) + a;
    get_int_dep(holonomies, hol, c);
  }

  Mark mark(e2.index(), c[0], c[1]);
  da = da_spec(d, a, mark); 

  if (report) {
    cout << "[" << d << ',' << a << "] " << mark << ' ';
    cout << '(' << (ee.index()/2) << ',' << (e2.index()/2) << ',';
    cout << ti.face_edge_index((v+2)%3) << ')' << endl; 
  }
}

void OrthoAngleSetup::get_da_spec(EdgeEnd const& ee, vector<Complex>& spec) const
{
  spec.resize(0); 
  EdgeIterator ei; 
  da_spec da; 

  // walk around this end
  for (ei = ee.edge_iterator(); ei.index < ee.order(); ei++) {
    get_da_spec(ei, da);
    spec.push_back(da.distance());
    spec.push_back(da.angle());
  }
}

void OrthoAngleSetup::get_da_spec(EdgeEnd const& ee, da_spec& da) const
{
  da = da_spec(orthodistance(ee), abs_oa(ee), Mark(ee.index(), 0));
}

void OrthoAngleSetup::get_da_spec(EdgeEnd const& ee, vector<da_spec>& spec) const
{
  spec.resize(ee.order());
  EdgeIterator ei; 
  for (ei = ee.edge_iterator(); ei.index < ee.order(); ei++) {
    get_da_spec(ei, spec[ei.index]);
  }
}

void OrthoAngleSetup::reduce_orthoangles()
{
  if (holonomies.size() != 2) {
    cerr << "invalid holonomies: " << PSeq(holonomies) << endl; 
    return; 
  }

  int i, n = _abs_orthoangle.size(); 
  for (i=0; i<n; i++) {
    torus_reduce(holonomies, _abs_orthoangle[i]); 
  }
}

void print_orthoangles(Triangulation* t, GroupPresentation* G)
{
  direct_edges(t);

  vector<Complex> orth;
  if (!get_edge_orthodistances(t, G, orth)) {
    cout << "Couldn't get the orthodistances!\n";
    return; 
  }

  OrthoAngleSetup OAS(t, orth);

  OAS.print(); 
}

// da_spec

// Wrapper class for OrthoAngleSetup

OrthoAngles::OrthoAngles(Triangulation* t, GroupPresentation* G)
{
  number_the_edge_classes(t); 
  direct_edges(t); 

  vector<Complex> orth;
  if (!get_edge_orthodistances(t, G, orth)) {
    cerr << "Couldn't get the orthodistances!\n";
    return; 
  }
  OAS = new OrthoAngleSetup(t, orth); 
}

void OrthoAngles::get_holonomies(vector<Complex>& H) const
{
  return OAS->get_holonomies(H);
}

void OrthoAngles::get_da_spec(EdgeClass* e, vector<Complex>& da_spec) const
{
  EdgeEnd ee(e, true); 
  OAS->get_da_spec(ee, da_spec);
}

void OrthoAngles::get_da_spec(EdgeIterator const& it, da_spec& da, bool report) const
{
  OAS->get_da_spec(it, da, report);
}

void OrthoAngles::get_da_spec(EdgeEnd const& ee, da_spec& da) const
{
  OAS->get_da_spec(ee, da); 
}

void OrthoAngles::get_da_spec(EdgeEnd const& ee, vector<da_spec>& spec) const
{
  OAS->get_da_spec(ee, spec);
}

void OrthoAngles::print() const
{
  OAS->print();
}

bool OrthoAngles::valid() const
{ 
  return OAS != 0 && OAS->valid(); 
}

OrthoAngles::~OrthoAngles()
{
  delete OAS;
}
