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
#include <iomanip>
#include "closed_geodesics.hh"
#include <algorithm>
#include "printable.hh"

using std::set;
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::setw;
using std::vector;
using std::min;
using std::max;

#define EPSILON O31_line::epsilon

class DirichletDomain {
  const WEPolyhedron* polyhedron; 
  vector<O31_vector> base_translate; 

  static O31_vector e0; 

  O31_vector normal(int i) const // inward normal
  { return e0 - base_translate[i]; }
public:
  DirichletDomain(const WEPolyhedron* p); 

  bool in_boundary(interval const& i) const;
  void compute_equivalents(list<interval>& equivalents, bool end) const;
  bool next_interval(interval& P) const;
  int normalize(O31_vector& point, O31_matrix& tr, FGWord& nw, int report=0) const;
  bool crossing_interval(interval& I) const; 
  void all_crossings(O31_line const& l, GSpec const& g, list<interval>& ivls) const; 
  void all_crossings(interval const& I, list<interval>& ivls) const;
  int find_line(list<interval> const& li, O31_line l, interval& found) const;
  int get_lifts(list<interval> const& il, double radius, set<line_wd>& S) const;  
};

static bool same(interval const& a, interval const& b, double eps)
{
  return close(a.end(0),b.end(0),eps) &&
    close(a.end(1),b.end(1),eps); 
}

static bool contains(list<interval> const& il, interval const& I)
{
  list<interval>::const_iterator i;
  for (i=il.begin(); i!=il.end(); ++i)
    if (same(*i, I, 10*EPSILON)) return true; 
  return false; 
}

static bool contains(list<interval> const& il, O31_vector const& pt)
{
  list<interval>::const_iterator i;
  for (i=il.begin(); i!=il.end(); ++i)
    if (close(i->end(1), pt, 10*EPSILON)) return true; 
  return false; 
}

inline static int signof(double x)
{
  if (x > EPSILON) return 1; 
  if (x <-EPSILON) return -1;
  return 0; 
}

bool DirichletDomain::in_boundary(interval const& I) const
{
  int i, n = base_translate.size();
  for (i=0; i<n; i++) 
    if (signof(I.midpoint() * normal(i))==0) return true;
  return false; 
}

/*
 * Compute_equivalents(equiv, end). 
 *
 * With end=true, two intervals are considered directly 
 * equivalent if their (final) endpoints are on the boundary of the
 * domain, and a face pairing carries one onto the other. 
 * With end=false, two intervals contained in the boundary 
 * of the domain are directly equivalent if a face pairing carries 
 * one to the other. 
 *
 * This function computes the equivalence class of the supplied 
 * intervals (in whichever sense). 
 *
 * For a single interval whose endpoint does not hit an edge 
 * of the domain, there will be just one (end) equivalent interval, 
 * namely the one it is carried onto by the face pairing of the
 * face continaing the endpoint. 
 */

void DirichletDomain::compute_equivalents(list<interval>& equivalents, bool end) const
{
  int i; 
  WEFace* face; 
  interval I;
  int guard=200; 
  list<interval>::iterator ivl = equivalents.begin(); 

  // We'll compute the immediate equivalents of 
  // each interval in turn, appending them to the list
  // until we catch up with the end (hopefully)!

  while (ivl != equivalents.end() && --guard > 0) {

    // For each face.. 
    for (face = polyhedron->face_list_begin.next, i=1;
	 face != &polyhedron->face_list_end;
	 face = face->next, i++) {

      // If face contains end/interval..
      if (signof((end ? ivl->end(1) : ivl->midpoint()) * normal(i))==0) {

	// Compute its partner under the face pairing.
	I = *ivl; 
	I.transform(*face->mate->group_element, face->mate->word);
	
	// Add it to the list.
	if (!contains(equivalents, I))
	  equivalents.push_back(I); 
      }
    }
    ++ivl; 
  }
  if (guard == 0) {
    cout << "Trouble in compute_equivalents!\n"; 
  }
}

/* 
 * Next_interval(I).
 *
 * I is an interval which is assumed to have endpoints 
 * on the boundary of the domain. To get the next interval
 * we use compute_equivalents to find an interval whose
 * outgoing endpoint points into the domain. 
 */

bool DirichletDomain::next_interval(interval& I) const
{
  list<interval> equivs;
  equivs.push_back(I); 

  compute_equivalents(equivs, true);

  int j, n = base_translate.size();
  list<interval>::iterator i; 

  // Look through the intervals.. 
  for (i=equivs.begin(); i!=equivs.end(); ++i) {

    // If end is on a face and points outward, this isn't the one for us.  
    for (j=1; j<n; ++j) {
      if (signof(i->end(1) * normal(j))) continue;
      if (signof(i->direction(1) * normal(j)) < 0) break; // Outward.
    }
    if (j==n) break; // End was inward pointing. 
  }

  // Did we find an inward pointing end?
  if (i == equivs.end()) {
    cout << "Couldn't find an inward pointing normal\n";
    return false; 
  }

  I = *i; // Yes

  // Next interval is extension of this one, chopped to the domain. 
  return crossing_interval(I);
}


bool DirichletDomain::crossing_interval(interval& I) const
{
  double lo = -1e10;
  double hi =  1e10;

  O31_vector N; 
  double si, ci;
  int i, sgn; 

  for (i = 1; i<base_translate.size(); ++i) {
    N = normal(i); 
    sgn = I.L.crossing_parameter(N, si, ci);
    if (sgn==0) continue; // face and line don't meet.
    
    if (sgn > 0) { // line entering domain
      if (si > lo) { lo = si; } 
    } else { // line leaving domain
      if (si < hi) { hi = si; } 
    }
  }

  I.set_range(asinh(lo), asinh(hi)); 

  return lo < hi;
}

void DirichletDomain::all_crossings(O31_line const& l, GSpec const& g, list<interval>& ivls) const
{
  interval t(l,0.,0.,FGWord(),g); 

  // get point on t closest to the origin
  O31_vector p = l.projection(e0);

  // translate t so it meets the domain
  if (normalize(p, t.L.T(), t.word) > 1) {
    // p was on an edge or vertex of the domain so perturb it
    O31_vector q(1.,.1,.05,.031);
    q.normalize();
    p = t.L.projection(q);
    normalize(p, t.L.T(), t.word); 
  }

  if (!crossing_interval(t)) {
    cout << "trouble in all_crossings(line)\n"; 
    return; 
  }

  list<interval> all_t0;
  all_t0.push_back(t);
  compute_equivalents(all_t0, false);

  double s[2], c[2], x, y; 

  int count=0;
  while (count < 50) {

    ivls.push_back(t);
  
    if (!next_interval(t)) break;

    if (contains(all_t0,t)) break; 
      
    ++count; 
  }

  if (count==50) {
    cout << "geodesic didn't close up after 50 crossings\n";
  } 
}

void DirichletDomain::all_crossings(interval const& I, list<interval>& ivls) const
{
  interval t(I); 

  // translate t so it starts inside the domain
  O31_vector p = t.end(0);
  normalize(p, t.L.T(), t.word); 

  if (!crossing_interval(t)) {
    cout << "trouble in all_crossings(interval)\n"; return; 
  }

  double s[2], c[2]; 

  double lo, hi;
  int count=0;
  while (count < 50) {

    lo = max(I.lo, t.lo);
    hi = min(I.hi, t.hi);

    ivls.push_back(interval(t.L, lo, hi, t.word, t.g));
  
    if (hi==I.hi) break; 

    // t = interval(I.L, hi, hi, I.word);

    if (!next_interval(t)) break; 

    ++count; 
  }

  if (count==50) {
    cout << "geodesic didn't close up after 50 crossings\n";
  } 
}


O31_vector DirichletDomain::e0(1.,0.,0.,0.);

DirichletDomain::DirichletDomain(const WEPolyhedron* p)
  : polyhedron(p)
{
  base_translate.resize(polyhedron->num_faces + 1); 

  base_translate[0] = e0; 

  int i = 1;
  WEFace* face;

  for (face = polyhedron->face_list_begin.next;
       face != &polyhedron->face_list_end;
       face = face->next)
    base_translate[i++] = *face->group_element * e0; 
}

void dirichlet_normalize(const WEPolyhedron* polyhedron, O31_vector& point, O31_matrix& tr, FGWord& nw, int report)
{
  DirichletDomain D(polyhedron);

  D.normalize(point, tr, nw, report); 
}

// dirichlet_normalize(point, tr, nw, report); 
//
// Finds an element N of the fundamental group such that N*point is inside polyhedron. 
// Replaces tr by N*tr, and nw by a word Nw*nw where Nw is a word representing N
// with respect to the (unsimplified) fundamental group generators. 

int DirichletDomain::normalize(O31_vector& point, O31_matrix& tr, FGWord& nw, int report) const
{

  WEFace *face, *nearest_face; 

  int num_base_translates = base_translate.size();

  FGWord prev_nw; 
  double cosh_min_dist, cosh_dist; // , org_dist, prev_org_dist = 1e6;
  int nearest_face_index; 
  O31_matrix ntr; 
  int i, maxtr = 200, n_closest; 
  // vector<int> closest; 

  if (report) cout << "Point: " << point << "\n"; 

  while (--maxtr) {

    cosh_min_dist = point[0]; 
    nearest_face = 0; 
    nearest_face_index = 0; 
    
    n_closest = 1; 
    // closest.resize(1); 
    // closest[0] = 0; 

    if (report) cout << "Cosh distances pt -> basepoint & translates:\n[" << cosh_min_dist; 

    // Find which base_translate is closest
    for (i = 1, face = polyhedron->face_list_begin.next;
	 face != &polyhedron->face_list_end; // i < num_base_translates; 
	 i++, face = face->next) {

      cosh_dist = -(base_translate[i] * point); 

      if (report) cout << ", " << cosh_dist;

      if (cosh_dist < cosh_min_dist - EPSILON) {
	nearest_face = face;
	nearest_face_index = i; 
	cosh_min_dist = cosh_dist; 
	n_closest = 1;
	// closest.resize(1); 
	// closest[0] = i; 
      } else if (cosh_dist < cosh_min_dist + EPSILON) {
	++n_closest;
	// closest.push_back(i);
      }
    }

    if (report) cout << "]\nNearest face: " << nearest_face_index << endl;

    if (nearest_face==0) {
      // if (closest.size() > 1) 
      // cout << "Point on boundary in dirichlet_normalize: " << PSeq(closest) << endl;
      return n_closest; 
    }

    // Apply transformation taking nearest_face to its mate
    ntr = *nearest_face->mate->group_element; 
    point = ntr * point; 
    tr = ntr * tr; 
    nw = nearest_face->mate->word * nw;

    if (report) cout << "Transformation * point: " << point << endl;
  }

  cout << "dirichlet_normalize: point outside domain after 200 steps\n";
  return 0; 
}

ostream& operator << (ostream& out, const interval& ivl)
{
  ios::fmtflags old = out.setf(ios::fixed); 
  out << "[";
  fwprint(out, ivl.lo, 0);
  out << ","; 
  fwprint(out, ivl.hi, 0);
  out.flags(old); 
  return out << "] " << ivl.word; 
}

bool operator == (interval const& a, interval const& b)
{
  return hyp_close(a.end(0), b.end(0), EPSILON) &&
    hyp_close(a.end(1), b.end(1), EPSILON);
}


#if 0
static bool face_contains(WEFace* face, interval const& i)
{
  O31_vector n = face->group_element->col(0);
  n[0] -= 1.0;
  return i.L.is_in_plane(n); 
}

static void compute_interval_image(interval const& i, WEFace* face, interval& img)
{
  // Map interval to its image in the paired face. 
  img.L    = *face->mate->group_element * i.L;
  img.word = (face->mate->word) * i.word; 
  if (&img == &i) return; 
  img.lo = i.lo; 
  img.hi = i.hi; 
  img.g = i.g; 
}

static void get_orbit(const WEPolyhedron *poly, interval i, list<interval>& l)
{
  WEFace *face, *onface[2]; 
  int num_faces = 0; 

  // count how many faces contain this edge; 
  for (face = poly->face_list_begin.next;
       face != &poly->face_list_end;
       face = face->next) {
    if (face_contains(face, i)) {
      onface[num_faces] = face; 
      num_faces++;
      if (num_faces==2) break; 
    }
  }

  if (num_faces==0) {
    l.push_back(i); 
    return; 
  }

  cout << "equiv ivls(2)\n";

  if (num_faces==1) {
    l.push_back(i);
    compute_interval_image(i, onface[0], i); 
    if (!face_contains(onface[0]->mate, i)) {
      cerr << "interval in face, but image not in paired face!\n"; 
      return; 
    }
    l.push_back(i); 
    return; 
  }

  // num_faces == 2. 

  // find the edge that contains both faces.
  WEEdge *edge, *common_edge=NULL; 
  for (edge = poly->edge_list_begin.next;
       edge != &poly->edge_list_end;
       edge = edge->next) {
    if (edge->f[left]==onface[0]) {
      if (edge->f[right]==onface[1]) {
	common_edge = edge; 
	break; 
      }
    }
    if (edge->f[left]==onface[1]) {
      if (edge->f[right]==onface[0]) {
	common_edge = edge; 
	break; 
      }
    }
  }
  if (!common_edge) {
    cerr << "interval appears to be in two faces without a common edge!\n"; 
    return; 
  }

  // get images of the interval in all edges identified with this one. 
  int j, n = common_edge->e_class->num_elements; 

  for (j=0, edge=common_edge; j<n; j++, edge=edge->neighbor[left]) {

    if (!face_contains(edge->f[left], i) ||
	!face_contains(edge->f[right], i)) {

      cerr << "interval in edge but image not in identified edge!\n";
      return; 
    }
    
    l.push_back(i);
    compute_interval_image(i, edge->f[left], i); 
  }
  if (edge!=common_edge) {
    cerr << "e_class->num_elements seems to be incorrect!\n";
  }
}

static void expand_to_full_set(const WEPolyhedron *poly, list<interval>& l)
{
  list<interval> l1;
  l1.splice(l1.end(), l); // transfer l to l1.

  list<interval>::const_iterator it;
  for (it = l1.begin(); it!=l1.end(); ++it)
    get_orbit(poly, *it, l);
}
#endif 

bool all_crossing_lifts(const WEPolyhedron *poly, O31_line const& L, GSpec const& g, list<interval>& l, bool full_set)
{
  DirichletDomain D(poly);
  D.all_crossings(L, g, l);
  int n0 = l.size(); 
  if (full_set) { 
    // expand_to_full_set(poly, l); 
    D.compute_equivalents(l,false);
    // if (l.size() > n0) cout << "added " << (l.size()-n0) << " equivalents\n";
  }
  return true;
}

bool all_crossing_lifts(const WEPolyhedron *poly, interval const& I, list<interval>& l, bool full_set)
{
  DirichletDomain D(poly);
  D.all_crossings(I, l);
  if (full_set) D.compute_equivalents(l,false);
  return true; 
}


// li should be a full set of intervals: see comments for all_crossing_lifts.

int DirichletDomain::find_line(list<interval> const& li, O31_line l, interval& found) const
{
  // translate l so it meets the domain
  O31_vector p = l.projection(e0);
  FGWord conj;
  normalize(p, l.T(), conj); 

  int res; 
  list<interval>::const_iterator it;
  for (it=li.begin(); it!=li.end(); it++) {
    res = compare(it->L, l);
    if (res) {
      found = *it; 
      return res; 
    }
  }
  return 0; 
}

static O31_line to_O31_line(GeodesicWord const& g)
{
  line l;
  if (fixed_points(g.mt.matrix,l) < 2) {
    cout << "big trouble converting geodesic to O31_line\n";
    l = line(Zero, Infinity);
  }
  return O31_line(l); 
}

bool find_geodesic(const WEPolyhedron* domain,
		   vector<GeodesicWord> const& geodesics, 
		   O31_line const& L, 
		   int start_g_num, int& gn, int& ori, FGWord& conj)
{
  // Make sure there is some point to carrying on. 
  int n = geodesics.size(); 
  if (n==0) return false; 

  DirichletDomain D(domain); 

  /* Find all images of L which intersect with the Dirichlet domain. */
  list<interval> il; 
  D.all_crossings(L, 0, il);
  D.compute_equivalents(il,false);

  if (!il.size()) return false; 

  /* At this point we could discard from il all images of L which
     do not minimize distance to the origin (basepoint of domain). */ 

  int i;
  O31_line G;
  list<interval>::const_iterator it; 
  for (i=start_g_num; i<n; i++) {
    G = to_O31_line(geodesics[i]); 
    for (it=il.begin(); it!=il.end(); ++it) {
      if ((ori = compare(it->the_line(), G, 10*EPSILON))) break; 
    }
    if (ori) break; 
  }
  if (i==n) return false; 

  gn = i; 
  conj = it->the_word(); 
  return true; 
}


// The list of intervals here is assumed to be composed of all intervals 
// for g1, followed by those for g2 etc. up to those for gn. 
// li should be a full set of intervals: see all_crossing_lifts. 

bool is_symmetry_of_geodesic_lifts(const WEPolyhedron* poly, list<interval> const& li, O31_matrix T, int report)
{
  DirichletDomain D(poly);

  list<interval>::const_iterator ia; 
  int gn = -1; 
  ia = li.begin();
  O31_line L;
  interval found; 
  int res, overall=true; 

  while (true) {

    // Skip to next unchecked geodesic (by number). 
    while (ia!=li.end() && ia->gnum()==gn) ia++; 
    if (ia==li.end()) break; 
    gn = ia->gnum(); 

    L = T * ia->the_line(); 

    res = D.find_line(li, L, found);
    if (!res) overall = false; 

    if (report) {
      if (res) cout << setw(4) << found.gnum() * res;
      else cout << "   ?"; 
    } 
  }
  return overall; 
}

double outradius(list<interval> const& l)
{
  double cosh_d_max = 1.0, c_d; 

  list<interval>::const_iterator it;
  for (it = l.begin(); it != l.end(); it++) {
    c_d = ((*it).end(0))[0];
    if (c_d > cosh_d_max) cosh_d_max = c_d; 
    c_d = ((*it).end(1))[0];
    if (c_d > cosh_d_max) cosh_d_max = c_d; 
  }
  return acosh(cosh_d_max); 
}



static int lexcmp(O31_vector const& a, O31_vector const& b, double eps)
{
  int i;
  for (i=0; i<4; ++i) {
    if (a[i] > b[i] + eps) return 1;
    if (a[i] < b[i] - eps) return -1; 
  }
  return 0; 
}

bool operator < (line_wd const& a, line_wd const& b)
{
  if (a.d < b.d - EPSILON) return true; 
  if (a.d > b.d + EPSILON) return false; 
  
  int res = lexcmp(a.l.end(0), b.l.end(0), EPSILON);
  if (res) { return res < 0; }
  
  return (lexcmp(a.l.end(1), b.l.end(1), EPSILON) < 0);
}

line_wd operator * (O31_matrix const& M, line_wd const& ln)
{
  return line_wd(M * ln.L(), ln.gnum);
}

int DirichletDomain::get_lifts(list<interval> const& il, double radius, set<line_wd>& S) const
{
  set<line_wd> leaf;

  // get the initial set of geodesics. 
  list<interval>::const_iterator j;
  for (j=il.begin(); j!=il.end(); ++j)
    leaf.insert(line_wd(j->L, j->gnum()));

  line_wd L, M;
  int i; 
  WEFace* face; 
  double d;

  while (leaf.size() > 0) {

    L = *leaf.begin();
    leaf.erase(leaf.begin()); 

    S.insert(L);
    d = L.distance() + EPSILON;

    // check "neighbors" of L
    for (face = polyhedron->face_list_begin.next, i=1;
	 face != &polyhedron->face_list_end;
	 face = face->next, i++) {

      M = *face->group_element * L; 
      if (M.distance() < d || M.distance() > radius) continue; 

      leaf.insert(M); 
    }
  }

  return S.size();
}

int get_lifts(const WEPolyhedron* poly, list<interval> const& il, double radius, set<line_wd>& S)
{
  DirichletDomain D(poly);
  return D.get_lifts(il, radius, S); 
}
