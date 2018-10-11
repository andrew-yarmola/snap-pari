#ifndef _polytope_T_
#define _polytope_T_
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
#include <iostream>

using std::list;
using std::vector;
using std::ostream;

template <class VS>
struct vertex;

template <class VS>
struct face {
  typedef typename list<vertex<VS> >::iterator vertex_p; 
  typedef typename VS::V VEC;

  list<vertex_p> I;
  int index;
  bool keep;
  VEC F;

  face(int f, VEC const& _F) : index(f), F(_F) {}
  void print_incidences(ostream& out) const;
};

template <class VS>
struct vertex {
  typedef typename list<face<VS> >::iterator face_p;
  typedef typename VS::V VEC;

  list<face_p> I; // incident faces, must remain sorted
  int index; 
  int cut_value;
  VEC V; 

  vertex(int i, VEC const& v) : index(i), V(v) {}
};

template <class VS>
struct edge {
  typedef typename list<vertex<VS> >::iterator vertex_p; 
  typedef typename list<face<VS> >::iterator face_p;
  typedef typename list<face_p>::const_iterator IF_p; 
  typedef typename list<vertex_p>::const_iterator IV_p; 
  typedef typename VS::V VEC;
  typedef typename VS::SC SC;

  vertex_p a,b;

  edge(vertex_p ap, vertex_p bp) : a(ap), b(bp) {}
  void cut(vertex_p v);
  VEC crossing(VEC const& F, bool& bad_cut) const; 
};

template <class VS>
struct nface {
  typedef typename list<face<VS> >::const_iterator face_cp;
  typedef typename list<vertex<VS> >::const_iterator vertex_cp; 
  typedef typename list<vertex_cp>::const_iterator CIV_p; 
  typedef typename VS::V VEC;

  list<face_cp> IF;
  list<vertex_cp> I;

  bool intersect_with(face_cp fp);
  bool get_dh_barycenter(VEC& c) const; 
};

template <class VS>
ostream& operator << (ostream& out, nface<VS> const& f);

template <class VS>
struct polytope {
  /* some versions of g++ cannot handle the following typedefs which
     are using the same name twice */
  //   typedef face<VS> face; 
  //   typedef vertex<VS> vertex; 
  //   typedef edge<VS> edge; 


  typedef list<nface<VS> > nface_list; 

  typedef typename list<face<VS> >::iterator face_p;
  typedef typename list<vertex<VS> >::iterator vertex_p; 
  typedef typename list<face<VS> >::const_iterator face_cp;
  typedef typename list<vertex<VS> >::const_iterator vertex_cp; 
  typedef typename list<edge<VS> >::const_iterator edge_p; 

  typedef typename list<face_p>::const_iterator IF_p; 
  typedef typename list<vertex_p>::const_iterator IV_p; 

  typedef typename list<face_cp>::const_iterator CIF_p; 
  typedef typename list<vertex_cp>::const_iterator CIV_p; 

  typedef typename VS::V VEC;
  typedef typename VS::M MAT;
  typedef typename VS::SC SC;
  
  list<face<VS> > FL; 
  list<vertex<VS> > VL; 
  list<edge<VS> > EL; 

  int dim; 
  int num_cuts; 

  bool bad_cut, bad_edge_cut; 

  // LI_set_queue* pending;

  polytope() : dim(0), num_cuts(0), bad_cut(false), bad_edge_cut(false) {} 
  // pending=new LI_set_queue; 
  polytope(int d); 
  // ~polytope() { delete pending; }

  void initialize(int d);
  void initialize(MAT const& FM); 

  face_p find_face(int i); 
  face_p add_face(VEC const& F);

  int cut(VEC F, face_cp *new_face=0, list<face<VS> > *keep_dead=0);

  void get_lattice(vector<nface_list>& L) const;

  void print_cuts() const; 
  void print_dhv() const; 
  void print() const; 
};

#if 0
template <class VS>
struct LI_set_queue {
  typedef face<VS> face; 
  typedef typename VS::V VEC;
  typedef typename VS::M MAT;
  typedef typename list<face<VS> >::const_iterator face_cp;

  vector<VEC> reduced;
  list<face> LI, LD;

  bool add(face const& F, face_cp* fp);
};
#endif

template<class IT>
inline bool face_less(IT a, IT b)
{ return a->index < b->index; }

template<class VS>
ostream& operator << (ostream& out, edge<VS> const& e);
template<class VS>
ostream& operator << (ostream& out, vertex<VS> const& v);
template<class VS>
ostream& operator << (ostream& out, face<VS> const& v);
template<class VS>
ostream& operator << (ostream& out, polytope<VS> const& p);

/* ---------------------------------- */

#include <cstdio>

using std::cout;
using std::endl; 


// check if face f occurs in the incidence list between 
// p and e, moving p up to where f occurs or the first 
// face greater than f. 

template<class FP, class IFP>
bool has_face(FP f, IFP& p, IFP const& e)
{
  while (p!=e && face_less(*p,f)) ++p; 
  return p!=e && *p==f; 
}

// compare a and b relative to r; want to determine relations between
// the set of faces common to r,a and common to r,b. 
// result = 0 if equal (a==b)
//          1 if r,a has more faces ie. (b < a)
//          2 if r,b has more faces ie. (a < b) 
//          3 if each has faces not in the other. (a,b not comparable)

template<class VS>
int compare_edges(vertex<VS> const& a, vertex<VS> const& b, vertex<VS> const& r)
{
  typedef typename polytope<VS>::IF_p IF_p;
  IF_p pr, pa=a.I.begin(), pb=b.I.begin(); 
  int result = 0; 
  for (pr=r.I.begin(); pr!=r.I.end() && result < 3; pr++) {
    if (has_face(*pr,pa,a.I.end())) {
      if (!has_face(*pr,pb,b.I.end())) result |= 1;
    } else {
      if (has_face(*pr,pb,b.I.end())) result |= 2;
    }
  }
  return result; 
}


template<class VS>
bool has_common_faces(edge<VS> const& e)
{
  typedef typename polytope<VS>::face_p FP; 
  typedef typename polytope<VS>::IF_p IF_p; 
  list<FP> const& fla = e.a->I;
  list<FP> const& flb = e.b->I;

  IF_p pa;
  IF_p pb=flb.begin();
  for (pa=fla.begin(); pa!=fla.end(); pa++) {
    if (has_face(*pa,pb,flb.end())) return true;
  }
  return false; 
}

// assume e is lexicographically higher than any 
// element of EL so for x in EL have (xa < xb) < (ea < eb)
// so that xa < eb is guaranteed.

template<class VS>
void insert_edge(list<edge<VS> >& EL, edge<VS> const& e)
{
  // if (!has_common_faces(e)) return; // no common faces so not an edge.

  typename list<edge<VS> >::iterator xp = EL.begin();
  int res; 
  while (xp!=EL.end()) {

    if (e.a == xp->a) {
      res = compare_edges(*e.b, *xp->b, *e.a);
    } else if (e.a == xp->b) {
      res = compare_edges(*e.b, *xp->a, *e.a);
    } else if (e.b == xp->b) {
      res = compare_edges(*e.a, *xp->a, *e.b);
    } else { // no vertices in common so not-comparable
      res = 3; 
    }

    if (res==2 || res==0) return; // faces(e) <= faces(*xp)
    if (res==1) { // faces(*xp) < faces(e)
      EL.erase(xp++); 
    }
    else ++xp; 
  }

  // if we get here e was bigger than or incomparable with anything in EL. 
  EL.push_back(e);
}

template<class VS>
void get_edges(list<vertex<VS> >& VL, list<edge<VS> >& EL0)
{
  typename polytope<VS>::vertex_p v1, v2; 
  list<edge<VS> > EL;

  // v1,v2 iterate over all possible edges in lexicographic order. 
  for (v1=VL.begin(); v1!=VL.end(); ++v1) {
    for (v2=v1, ++v2; v2!=VL.end(); ++v2) {
      insert_edge(EL,edge<VS>(v1,v2)); 
    }
  }
  EL0.splice(EL0.end(), EL); 
}

template<class VS>
typename polytope<VS>::face_p polytope<VS>::find_face(int i) 
{
  face_p fp=FL.begin();
  while (fp!=FL.end() && fp->index < i) ++fp;
  if (fp->index != i) { cout << "face " << i << " does not exist!\n"; }
  return fp; 
}

template<class VS>
typename polytope<VS>::face_p polytope<VS>::add_face(VEC const& F)
{
  int fi = FL.back().index + 1; // last face has biggest index
  FL.push_back(face<VS>(fi, F)); 
  face_p fp=FL.end();
  return --fp; 
}

template<class VS>
polytope<VS>::polytope(int d)
  : num_cuts(0), bad_cut(false), bad_edge_cut(false) //, pending(0)
{
  initialize(d); 
}

template<class VS>
void polytope<VS>::initialize(int d)
{
  initialize(id_mat(d));
}

template<class VS>
void polytope<VS>::initialize(MAT const& FM)
{
  int d = rows(FM); 

  bool inv; 
  MAT VM = inverse(FM, &inv); 
  if (!inv) {
    printf("non-invertible matrix of initial faces supplied to polytope\n");
    return; 
  }

  dim = d; 

  int i, j; 
  for (i=0; i<d; i++)
    FL.push_back(face<VS>(i, row(FM,i)));
  for (i=0; i<d; i++)
    VL.push_back(vertex<VS>(i, col(VM,i)));

  vertex_p vi;
  face_p fi; 
  for (i=0,vi=VL.begin(); vi!=VL.end(); i++,vi++) {
    for (j=0,fi=FL.begin(); fi!=FL.end(); j++,fi++) {
      if (i==j) continue; 
      vi->I.push_back(fi);
      fi->I.push_back(vi);
    }
  }

  // every pair of vertices has an edge in the standard polytope. 
  vertex_p v1, v2; 
  for (v1=VL.begin(); v1!=VL.end(); ++v1) {
    for (v2=v1, ++v2; v2!=VL.end(); ++v2) {
      EL.push_back(edge<VS>(v1,v2)); 
    }
  }
}

template<class VS>
ostream& operator << (ostream& out, edge<VS> const& e)
{
  return cout << '(' << e.a->index << ',' << e.b->index << ')'; 
}

template<class VS>
void edge<VS>::cut(vertex_p v)
{
  list<face_p> const& fla = a->I;
  list<face_p> const& flb = b->I;

  // set up the incidences of the new vertex
  IF_p pa;
  IF_p pb=flb.begin();
  for (pa=fla.begin(); pa!=fla.end(); pa++) {
    if (has_face(*pa,pb,flb.end())) { // face *pa incident with edge, hence new vertex.
      v->I.push_back(*pa);
      (*pa)->I.push_back(v); 
    }
  }
  v->cut_value = 0; 

  // cut this edge. 
  if (a->cut_value == -1) a=v;
  else b=v; 
}

template <class VS>
typename VS::V edge<VS>::crossing(VEC const& F, bool& bad_cut) const
{
  SC Fa = dotprod(F,a->V), Fb = dotprod(F,b->V);
  if (size(Fb-Fa) < VS::ne_EPS) { // Was 2*VS::ne_EPS 
    cout << "undefined edge cut\n"; 
    VEC V(VS::half * a->V + VS::half * b->V); 
    normalize(V);
    return V;
  }
  SC t = Fa/(Fa-Fb);
  VEC V((VS::one-t) * a->V + t * b->V);

  // check for bad edge cuts. 
  if (!bad_cut) {
    double err = (VS::eq_EPS/size(Fa-Fb)) *
      size(norm(b->V - a->V)/norm(V)); 
    if (err > VS::ne_EPS) bad_cut = true; 
  }

  normalize(V);
  return V;
}

template <class VS>
int polytope<VS>::cut(VEC F, face_cp *new_face, list<face<VS> > *keep_dead)
{
  vertex_p vp; 
  face_p fp;
  typename list<vertex_p>::iterator ivp;
  typename list<face_p>::iterator ifp;

  bool tracing = (num_cuts==-1); 
  if (tracing) cout << "Tracing is on!\n";

  normalize(F); 

#if 0
  if (pending) {
    face_cp fp; 
    if (pending->add(F, &fp)) { 
      // Have initial simplex. 
    } else { 
      // No initial simplex yet. 
    }
    if (new_face) *new_face = fp; 
  }
#endif

  // set cut_values
  bool in=false, out=false; 
  SC x; 
  for (vp=VL.begin(); vp!=VL.end(); ++vp) {
    x = dotprod(F,vp->V);
    if (size(x) < VS::eq_EPS) {
      vp->cut_value = 0; 
    } else {
      if (size(x) < VS::ne_EPS && !bad_cut) {
	cout << "Warning: cut number " << (num_cuts+1) << " is bad, ";
	cout << "size = " << size(x) << endl; 
	cout << vp->V << endl; 
	bad_cut=true; 
      }
      if (x > VS::zero) { vp->cut_value = 1; in=true; }
      else { vp->cut_value = -1; out=true; }
    }
  }


  // check for non-trivial cut
  if (!in) {
    cout << "cut would result in empty polyhedron (not done)\n";
    return -1;
  }
  if (!out) return 1; // trivial cut. 

  // REMOVE DEAD FACES 

  // decide which faces to keep 
  // and remove dead vertices from incidence lists
  for (fp=FL.begin(); fp!=FL.end(); fp++) {
    fp->keep = false; 
    for (ivp=fp->I.begin(); ivp!=fp->I.end();) {
      if ((*ivp)->cut_value > 0) { 
	fp->keep = true; 
	ivp++; 
      } else if ((*ivp)->cut_value < 0) {
	fp->I.erase(ivp++);
      } else {
	++ivp; 
      }
    }
  }

  // remove dead faces from incidence lists. (and get max vertex number)
  int vn = 0;
  for (vp=VL.begin(); vp!=VL.end(); ++vp) {
    if (vp->index > vn) vn=vp->index; // get max vertex number 
    if (vp->cut_value==1) continue; 
    for (ifp=vp->I.begin(); ifp!=vp->I.end();) {
      if (!(*ifp)->keep) vp->I.erase(ifp++);
      else ifp++; 
    }
  }

  // erase dead faces themselves. 
  for (fp=FL.begin(); fp!=FL.end();) {
    if (!fp->keep) {
      if (keep_dead)
	keep_dead->splice(keep_dead->end(), FL, fp++);
      else 
	FL.erase(fp++); 
    }
    else fp++;
  }

  // cut edges, creating new vertices & remove dead edges
  list<vertex<VS> > NVL;
  typename list<edge<VS> >::iterator ep = EL.begin(); 
  bool had_bad_edge_cut = bad_edge_cut; 
  while (ep!=EL.end()) {
    if (ep->a->cut_value * ep->b->cut_value == -1) {
      NVL.push_back(vertex<VS>(++vn, ep->crossing(F, bad_edge_cut))); // make new vertex.
      vp = NVL.end(); --vp; // get iterator pointing to it. 
      ep->cut(vp); // set up its face list and make it an end of this edge. 
      ++ep; 
    } else if (ep->a->cut_value + ep->b->cut_value <= 0) {
      EL.erase(ep++);
    } else {
      ++ep; 
    }
  }

  if (bad_edge_cut && !had_bad_edge_cut) {
    cout << "Warning: bad edge cut during cut " << (num_cuts+1) << endl; 
  }

  // remove dead vertices now that edges don't point to them anymore. 
  vp=VL.begin();
  while (vp!=VL.end()) {
    if (vp->cut_value==-1) {
      VL.erase(vp++);
    } else {
      ++vp;
    }
  }

  // find new edges from zero vertex set
  vertex_p vp2; 
  list<edge<VS> > NEL;
  // iterate over possible edges from zero vertices to zero vertices. 
  for (vp=VL.begin(); vp!=VL.end(); ++vp) {
    if (vp->cut_value > 0) continue; 
    for (vp2=vp, ++vp2; vp2!=VL.end(); ++vp2) {
      if (vp2->cut_value > 0) continue; 
      insert_edge(NEL,edge<VS>(vp,vp2)); 
    }
  }
  // iterate over possible edges from zero vertices to new zero vertices. 
  for (vp=VL.begin(); vp!=VL.end(); ++vp) {
    if (vp->cut_value > 0) continue; 
    for (vp2=NVL.begin(); vp2!=NVL.end(); ++vp2) {
      insert_edge(NEL,edge<VS>(vp,vp2)); 
    }
  }
  // iterate over possible edges between new zero vertices. 
  for (vp=NVL.begin(); vp!=NVL.end(); ++vp) {
    for (vp2=vp, ++vp2; vp2!=NVL.end(); ++vp2) {
      insert_edge(NEL,edge<VS>(vp,vp2)); 
    }
  }
  EL.splice(EL.end(), NEL); 

  // add the new face
  fp = add_face(F); 
  if (new_face) *new_face=fp; // tell caller about it if required. 

  // add incidence with new face to each existing and new zero vertex
  for (vp=VL.begin(); vp!=VL.end(); ++vp) {
    if (vp->cut_value > 0) continue; 
    vp->I.push_back(fp); 
    fp->I.push_back(vp); 
  }
  for (vp=NVL.begin(); vp!=NVL.end(); ++vp) {
    vp->I.push_back(fp); 
    fp->I.push_back(vp); 
  }

  // splice in the new zero vertices
  VL.splice(VL.end(), NVL); 

  num_cuts++; 

  if (dim <= 3) return 0; 

  // check if we've created any bigons!
  static bool bigon_warning_given = false; 
  if (bigon_warning_given) return 0; 
  for (fp=FL.begin(); fp!=FL.end(); fp++) {
    if (fp->I.size()==2) {
      cout << "Warning: after " << num_cuts << " cuts polytope contains a bigon\n";
      bigon_warning_given = true; 
      break; 
    }
  }

  return 0; // non-trivial cut
}

template <class VS>
ostream& operator << (ostream& out, vertex<VS> const& v)
{
  out << v.index << '=' << v.V << ' ';
  typename polytope<VS>::IF_p fp;
  for (fp=v.I.begin(); fp!=v.I.end(); ++fp) {
    if (fp!=v.I.begin()) out << ',';
    out << (*fp)->index;
  }
  return out; 
}

template <class VS>
ostream& operator << (ostream& out, face<VS> const& v)
{
  out << v.index << '=' << v.F << ' ';
  typename polytope<VS>::IV_p vp;
  for (vp=v.I.begin(); vp!=v.I.end(); ++vp) {
    if (vp!=v.I.begin()) out << ',';
    out << (*vp)->index;
  }
  return out; 
}

template <class VS>
void face<VS>::print_incidences(ostream& out) const
{
  out << '(';
  typename polytope<VS>::IV_p vp;
  for (vp=I.begin(); vp!=I.end(); ++vp) {
    if (vp!=I.begin()) out << ',';
    out << (*vp)->index;
  }
  out << ')';
}

template <class VS>
ostream& operator << (ostream& out, polytope<VS> const& p)
{
  typename polytope<VS>::face_cp fp;
  // print faces. 
  for (fp=p.FL.begin(); fp!=p.FL.end(); fp++)
    out << '(' << (*fp) << ')' << ' ';
  out << endl; 
  // print vertices
  typename polytope<VS>::vertex_cp vp; 
  for (vp=p.VL.begin(); vp!=p.VL.end(); vp++) 
    out << '(' << (*vp) << ')' << ' ';
  out << endl; 
  // edges
  typename polytope<VS>::edge_p ep; 
  for (ep=p.EL.begin(); ep!=p.EL.end(); ep++)
    out << (*ep) << ' '; 
  out << endl; 
  return out; 
}

template <class VS>
void polytope<VS>::print() const
{ 
  face_cp fp;
  // print faces. 
  cout << "faces\n";
  for (fp=FL.begin(); fp!=FL.end(); fp++)
    cout << '(' << (*fp) << ')' << ' ';
  cout << endl; 
  // print vertices
  cout << "vertices\n";
  vertex_cp vp; 
  for (vp=VL.begin(); vp!=VL.end(); vp++) 
    cout << '(' << (*vp) << ')' << ' ';
  cout << endl; 
  // edges
  cout << "edges\n";
  edge_p ep; 
  for (ep=EL.begin(); ep!=EL.end(); ep++)
    cout << (*ep) << ' '; 
  cout << endl; 
}

template <class VS>
void polytope<VS>::print_cuts() const
{
  vertex_cp vp; 
  for (vp=VL.begin(); vp!=VL.end(); vp++) {
    cout << vp->index << " -> " << vp->cut_value << endl;
  }
}

template <class VS>
void polytope<VS>::print_dhv() const
{
  vertex_cp vp; 
  VEC v(VS::vector(dim-1)); 
  for (vp=VL.begin(); vp!=VL.end(); vp++) {
    dehomogenized_copy(v, vp->V); 
    cout << vp->index << '=' << v << endl;
  }
  cout << endl; 
  face_cp fp;
  for (fp=FL.begin(); fp!=FL.end(); fp++) {
    fp->print_incidences(cout);
    cout << ' ';
  }
  cout << endl; 
}


// compare a and b; want to determine relations between
// the set of vertices of a and b. 
// result = 0 if equal (a==b)
//          1 if a has more vertices ie. (b < a)
//          2 if b has more vertices ie. (a < b) 
//          3 if each has vertices not in the other. (a,b not comparable)

template <class VS>
int compare(nface<VS> const& a, nface<VS> const& b)
{
  typename polytope<VS>::CIV_p pa=a.I.begin(), pb=b.I.begin(); 
  int res=0; 
  while (res != 3 && (pa!=a.I.end() || pb!=b.I.end())) {
    if (pa==a.I.end()) { res |= 2; break; } // a finished first so b bigger
    if (pb==b.I.end()) { res |= 1; break; } // b finished first so a bigger
    if ((*pa)->index < (*pb)->index) { res |= 1; ++pa; continue; } // a bigger
    if ((*pb)->index < (*pa)->index) { res |= 2; ++pb; continue; } // b bigger
    ++pa; ++pb; 
  }
  return res; 
}

template <class VS>
void remove_non_maximal(list<nface<VS> >& L) 
{
  int res; 
  typename list<nface<VS> >::iterator i, j;
  for (i=L.begin(); i!=L.end();) {
    for (j=i, j++; j!=L.end();) {
      res = compare(*i, *j); // 1 if j<i, 2 if i<j, 0 if i==j
      if (res==0||res==1) L.erase(j++); 
      else if (res==2) { L.erase(i++); break; }
      else j++; 
    }
    if (j==L.end()) i++; 
  }
}

template <class VS>
bool nface<VS>::intersect_with(face_cp fp)
{
  typename list<vertex_cp>::iterator pa=I.begin(); 
  typename polytope<VS>::IV_p pb=fp->I.begin(); 
  int isize=I.size(); 
  while (pa!=I.end()) {
    if (pb==fp->I.end()) 
      { I.erase(pa, I.end()); break; } // b finished so remove rest of a
    if ((*pa)->index < (*pb)->index) { I.erase(pa++); continue; }
    if ((*pb)->index < (*pa)->index) { ++pb; continue; }
    ++pa; ++pb; 
  }
  if (I.size() < isize) {
    IF.push_back(fp); 
    return true; 
  }
  return false; 
}

template <class VS>
void polytope<VS>::get_lattice(vector<nface_list>& L) const
{
  // make sure L is the right size. 
  L.resize(dim);

  // create the top dimensional cell. 
  int n = dim-1; 
  vertex_cp vp; 
  L[n].push_back(nface<VS>()); 
  for (vp=VL.begin(); vp!=VL.end(); vp++)
    L[n].back().I.push_back(vp); 

  // add the faces (codim 1 cells)
  --n; 
  face_cp fp;
  IV_p ivp; 
  for (fp=FL.begin(); fp!=FL.end(); fp++) {
    L[n].push_back(nface<VS>());
    L[n].back().IF.push_back(fp); 
    for (ivp=fp->I.begin(); ivp!=fp->I.end(); ++ivp)
      L[n].back().I.push_back(*ivp); 
  }

  // get nfaces of all dimensions. 
  typename nface_list::iterator nfp; 
  for (;n>0; --n) { 

    // get n-1 cells by intersecting n-cells with faces. 
    for (nfp=L[n].begin(); nfp!=L[n].end(); ++nfp) {
      
      fp = nfp->IF.back(); // will do intersection with faces not done yet
      while (1) {
	++fp; 
	L[n-1].push_back(*nfp); // put a copy of an n-face at the end. 
	
	// find a non-trivial intersection
	while (fp != FL.end() && !L[n-1].back().intersect_with(fp)) ++fp; 
	if (fp==FL.end()) { 
	  L[n-1].pop_back(); 
	  break; 
	}
      }
    }
    remove_non_maximal(L[n-1]); 
  }
}

template <class VS>
ostream& operator << (ostream& out, nface<VS> const& f)
{
  typename polytope<VS>::CIV_p i; 
  out << '('; 
  for (i=f.I.begin(); i!=f.I.end(); i++) {
    if (i!=f.I.begin()) out << ',';
    out << ((*i)->index); 
  }
  return out << ')'; 
}

template <class VS>
void print(vector<list<nface<VS> > > const& L)
{
  int i, n=L.size();
  typename list<nface<VS> >::const_iterator fi; 
  for (i=n-1; i>=0; i--) {
    if (i==n-1) cout << i << "-cell:  ";
    else cout << i << "-faces: "; 
    for (fi=L[i].begin(); fi!=L[i].end(); ++fi)
      cout << (*fi) << ' ';
    cout << endl; 
  }
}

template <class VS>
bool nface<VS>::get_dh_barycenter(VEC& c) const
{
  if (!I.size()) return false; 

  // get the (dehomogenized) dimension right. 
  int d = dim(I.back()->V) - 1;
  if (dim(c) != d) c = VEC(d); 

  c *= VS::zero; 

  // get average of vertices. 
  VEC dhv(VS::vector(d)); 
  CIV_p ivp; 
  for (ivp=I.begin(); ivp!=I.end(); ++ivp) {
    if (size((*ivp)->V[d]) < VS::ne_EPS) return false; 
    dehomogenized_copy(dhv, (*ivp)->V); 
    c += dhv; 
  }
  c /= VS::scalar(I.size()); 
  return true; 
}

template <class VS>
void print_barycenters(vector<list<nface<VS> > > const& L)
{
  int i, n=L.size();
  typename VS::V c; 
  typename list<nface<VS> >::const_iterator fi; 
  for (i=n-1; i>=0; i--) {
    if (i==n-1) cout << i << "-cell:  ";
    else cout << i << "-faces: "; 
    for (fi=L[i].begin(); fi!=L[i].end(); ++fi) {
      fi->get_dh_barycenter(c); 
      cout << c << ' ';
    }
    cout << endl; 
  }
}

#if 0
template <class VEC>
int max_entry(VEC const& V, double& max_size)
{
  int i, p=0, d=dim(V); 
  double sz; 

  max_size=0.; 
  for (i=0; i<d; ++i) {
    sz = size(reduced[r][i]);
    if (sz > max_size) 
      { max_size=sz; p = i; }
  }
  return p;
}

template <class VS>
bool LI_set_queue<VS>::add(VEC const& F, face_cp* fp)
{
  int d = dim(F);

  // Check if we're already full
  if (LI.size()==d) {
    LD.push_back(face(LD.size()+d,F));
    fp = LD.end(); --fp; 
    return true; 
  }

  VEC RV(F); 
  int n=reduced.size(), p, r;
  double max_size;

  // Reduce RV
  for (r=0; r<n; ++r) {
    p = max_entry(reduced[r], max_size); 
    RV -= (RV[p]/reduced[r][p])*reduced[r];
  }

  max_entry(RV, max_size);
  if (max_size < VS::eq_EPS) {
    LD.push_back(face(LD.size()+d,F));
    fp = LD.end(); --fp; 
    return false; 
  }

  reduced.push_back(RV);
  LI.push_back(face(LI.size(),F));
  fp = LI.end(); --fp; 
  return LI.size()==d;
}
#endif

#endif
