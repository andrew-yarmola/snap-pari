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

#include "tube.hh"
#include "connected_faces.hh"
#include "snappea/kernel.h"

void connected_faces::set_pairings(vector<Ortholine> const& ol)
{
  int i, n_ort = ol.size();
  pairing_matrices.resize(n_ort); 
  for (i=0; i<n_ort; i++) {
    pairing_matrices[i] = O31_matrix(ol[i]); 
  }
}   

static void copy_to_point_vector(uc_polygon const& F, vector<point>& vp)
{
  uc_polygon::const_iterator k; 

  vp = vector<point>(); // Resize to zero. 
  vp.reserve(F.size()); // Ensure there will be no need to realloc. 
  for (k=F.begin(); k!=F.end(); k++)
    vp.push_back(point(*k)); 
}

static void copy_to_point_vector(uc_polygon const& F, vector<point>& vp, O31_matrix const& M)
{
  uc_polygon::const_iterator k; 

  vp = vector<point>(); // Resize to zero. 
  vp.reserve(F.size()); // Ensure there will be no need to realloc. 
  for (k=F.begin(); k!=F.end(); k++)
    vp.push_back(M * point(*k)); 
}

static int compare_point_vectors(vector<point>& a, vector<point>& b, double& max_dis)
{
  bool match;
  double dis; 
  max_dis = 0.0; 

  if (a.size()==0 || b.size()==0) return 0; 

  int i=0, j, n = a.size();

  // Find distance of the closest pair of adjacent points in a. 
  double cutoff = alt_distance(a.front(), a.back());  
  for (j=0; j<n-1; j++) {
    dis = alt_distance(a[j],a[j+1]);
    if (dis < cutoff) cutoff = dis; 
  }
  cutoff /= 5.0;
  // Each point of b must be significantly closer to a point of a than
  // the adjacent point of a is. 

  while (i<n) {
    match = false; 
    for (j=i; j<b.size(); j++) {
      dis = alt_distance(a[i],b[j]);
      if (dis < cutoff) {
	// Move element of b to the corresponding position of the matched element of a. 
	swap(b[i],b[j]); 
	i++;
	match = true;

	// Keep track of how bad the worst match was. 
	if (dis > max_dis) max_dis = dis; 
	break;
      }
    }
    if (!match) {
      // Move unmatched element of a to the end. 
      // Reduce n so we don`t try to match it again. 
      n--; 
      swap(a[i],a[n]);
      max_dis = cutoff; 
    }
  }
  return n; 
}

static bool faces_match(uc_polygon const& f, uc_polygon const& bf, O31_matrix const& M, int report)
{
  int m; 
  vector<point> back_face, face;
  double dist;

  // Skip empty faces. 
  if (f.size()==0 && bf.size()==0) {
    if (report) cout << "are empty.\n"; 
    return true;
  }

  copy_to_point_vector(f, face);
  copy_to_point_vector(bf, back_face, M); 

  // Compare face and back_face. 
  m = compare_point_vectors(face,back_face,dist); 
  if (m == face.size() && m == back_face.size()) {
    if (report) {
      if (dist < tube_face::boundary_epsilon) cout << "match\n";
      else cout << "match, with error " << dist << endl;
    }
    return true; 
  } 

  if (report) {
    cout << back_face.size()-m << " unmatched on first, ";
    cout << face.size()-m << " unmatched on second, ";
    cout << "error exceeds " << dist << endl;
  }
  return false; 
}



bool connected_faces::compute_connected_faces(vector<tube_face>& faces)
{
  int i, j, k, l, n_ort = faces.size()/2;

  vector<int> subface_count(n_ort);
 
  int bcount = 0; // Number of face boundaries so far. 

  uc_polylist b, b0, b1; // Boundaries in b0 pair with those in b1.

  for (i=0; i<n_ort; i++) {
    b = faces[2*i].get_boundary(); 
    b0.splice(b0.end(), b); 

    b = faces[2*i+1].get_boundary();
    b1.splice(b1.end(), b); 

    // Check the count. 
    if (b1.size() != b0.size()) {
      cout << "Face and partner seem to split into different numbers of subfaces.\n";
      cout << "i = " << i << " sizes " << b0.size() << ' ' << b1.size() << endl; 
      return false; 
    }

    // Store the subface count. 
    subface_count[i] = b0.size() - bcount; 
    bcount = b0.size(); 
  }

  // Make FL and face_num vectors of the appropriate sizes. 
  FL = vector<uc_polygon>(2*bcount);
  face_num = vector<int>(bcount); 

  uc_polygon tmp; 
  uc_polylist::const_iterator p0 = b0.begin(), p1 = b1.begin(); 

  for (i = 0, j = 0; i < n_ort; i++) {

    for (k=0; k<subface_count[i]; k++) {
      face_num[j/2+k] = i; 
    }

    if (subface_count[i] == 0) continue; 

    if (subface_count[i] == 1) {
      FL[j] = *p0; 
      p0++; j++; 
      FL[j] = *p1; 
      p1++; j++; 
      continue; 
    } 

    // OK so face splits. 
    // Copy b0 faces. 
    for (k=0; k<subface_count[i]; k++) {
      FL[j+2*k] = *p0; p0++; 
    }

    // Insert b1 faces in the appropriate place.
    for (k=0; k<subface_count[i]; k++) {
      // Extract a b1 face. 
      tmp = *p1; p1++; 

      // Find out where it goes. 
      for (l=0; l<subface_count[i]; l++) {
	// Does it belong with this one? 
	if (faces_match(tmp, FL[j+2*l], pairing_matrices[i], 0)) {
	  FL[j+2*l+1].splice(FL[j+2*l+1].end(), tmp); 
	  break; 
	}
      }
      if (l==subface_count[i]) {
	cout << "Not all subfaces were matched.\n"; 
      }
    }
    j += 2*subface_count[i]; 
  }

  // Get starting indices in FL of the subfaces coming from each face. 
  vector<int> subface_start(2*n_ort); 
  j = 0; 
  for (i=0; i<n_ort; i++) {
    subface_start[2*i] = j; 
    subface_start[2*i+1] = j+1; 
    j += 2*subface_count[i]; 
  }

  // Renumber the edges. 
  uc_polygon::iterator vx; 
  for (j=0; j<2*bcount; j++) {
    for (vx=FL[j].begin(); vx!=FL[j].end(); vx++) {
      (*vx).mark_p = subface_mark_number((*vx).mark_p, *vx,subface_start,subface_count);
      (*vx).mark_n = subface_mark_number((*vx).mark_n, *vx,subface_start,subface_count);
    }
  }

  return true;
}

static void uc_print(uc_vertex const& v, vector<Complex> const& H)
{
  v.uc_print(cout, H);
  cout << ' ' << Mark(v.mark_p) << ' ' << Mark(v.mark_n) << '\n';
}

void connected_faces::print_connected_faces() const
{
  vector<Complex> H(2); 
  T->get_holonomies(H[0],H[1]); 

  int i, n_fac = FL.size();

  uc_polygon::const_iterator j; 

  for (i=0; i<n_fac; i++) {
    
    cout << "Face " << i << '\n'; 

    for (j=FL[i].begin(); j!=FL[i].end(); j++)
      uc_print(*j, H);
  }
}


/* Once a tube has been successfully computed we can extract the combinatorial 
   information and use it to construct an ideal triangulation for the manifold
   with the core geodesic of the tube drilled out. The relevant starting information 
   is all contained in the connected faces. The numbering scheme for connected 
   faces is necessarily different for that used for faces. The relation between
   the two schemes is given by a vector<int>. The following code is all ultimately
   concerned with drilling. In it we assume that the connected faces have been
   stored in a component called FL and the relation between the numbering schemes
   in a vector<int> component called face_num. 
*/ 

static int which_corner(const point& p, uc_polygon const& face, uc_vertex& v)
{
  uc_polygon::const_iterator vx;
  int vn;
  for (vx = face.begin(), vn = 0; vx != face.end(); vx++, vn++) {
    double d = alt_distance(p,*vx); 
    if (d < tube_face::boundary_epsilon) {
      v = *vx; 
      return vn;
    }
  }
  return -1; // p not found in face. 
}

face_corner connected_faces::adjacent_corner(point const& p, Mark const& nbr, uc_vertex& v) const
{
  // Get face and transformation. 
  int fnum = nbr.end_index(); 

  int vnum = which_corner(trans(-nbr) * p, FL[fnum], v);
  return face_corner(fnum, vnum); 
}

O31_matrix connected_faces::pairing_matrix(int face) const
{
  int fn = face/2; 
  if (!(face % 2)) return pairing_matrices[fn]; // Even. 
  else return inverse(pairing_matrices[fn]); // Odd. 
}

O31_matrix connected_faces::cf_pairing_matrix(int face) const
{
  int fn = face_num[face/2]; 
  if (!(face % 2)) return pairing_matrices[fn]; // Even. 
  else return inverse(pairing_matrices[fn]); // Odd. 
}

face_corner connected_faces::paired_corner(const point& p, int face, uc_vertex& v) const
{
  int paired_face = face ^ 1;
  int vnum = which_corner(cf_pairing_matrix(face) * p, FL[paired_face], v);
  return face_corner(paired_face, vnum);
}

static void reverse_face(uc_polygon& F)
{
  reverse(F.begin(), F.end());
  uc_polygon::iterator i; 
  for (i=F.begin(); i!= F.end(); i++)
    swap((*i).mark_p, (*i).mark_n); 
}

void connected_faces::orient_connected_faces()
{
  vector<Complex> H(2); 
  T->get_holonomies(H[0],H[1]); 

  int i, j, n_fac = FL.size();

  if (!n_fac) return; 

  bool same_or, all_done; 
  face_corner nbr; 
  uc_vertex pv; 
  
  vector<int> done(n_fac); 
  for (i=0; i<n_fac; i++) done[i] = 0; 
  uc_polygon::const_iterator v; 

  // 1 says face has been oriented, 2 says all its neighbors have been oriented as well. 
  done[0] = 1; 
  while (1) {
    all_done = true; 

    // Find a face which has been oriented, but whose neighbors haven`t. 
    for (i=0; i<n_fac; i++) {
      if (done[i]==1) break; 
      if (!done[i]) all_done = false; 
    }
    if (i==n_fac) break; 

    // Orient its neighbors. 
    for (v=FL[i].begin(), j=0; 
	 v!=FL[i].end(); 
	 ++v, ++j) {

      nbr = adjacent_corner(*v, v->mark_p, pv);
      if (nbr.vertex==-1) { 
	cout << "Couldn't find neighbor: "; 
	uc_print(*v, H); 
      }
      same_or = (Mark(i,-Mark(v->mark_p)) == Mark(pv.mark_n)); 

      if (done[nbr.face]) {
	if (!same_or) {
	  cerr << "Unable to orient faces of the tube.\n"; 
	  return; 
	}
      } else {
	if (!same_or)
	  reverse_face(FL[nbr.face]);
	done[nbr.face] = 1; 
      }
    }
    done[i] = 2; 
  }

  if (!all_done)
    cerr << "Unable to orient: tube appears to be disconnected.\n"; 

}

static void print_face_corners(face_corner const& f, face_corner const& p, 
			       face_corner const& a, face_corner const& n)
{
  cout << f << ' ' << p << ' ' << a << ' ' << n << endl; 
}

void connected_faces::print_face_corners() const
{
  int n_fac = FL.size(); 
  if (!n_fac) {
    cout << "Problem computing vertices of the tube faces.\n";
    return; 
  }

  int i, j; 
  face_corner fc, fc0;
  uc_polygon::const_iterator v0; 
  uc_vertex v, av; 

  int nv; 
  for (i=0; i<n_fac; i+=2) {

    // For each pair of faces print the pairing. 
    nv = FL[i].size(); 
    cout << i << ":(";
    for (j=0; j<nv; ++j) {
      if (j>0) cout << ','; 
      cout << j; 
    }
    cout << ") -> " << (i+1) << ":(";
    for (j=0, v0=FL[i].begin(); j<nv; ++j, ++v0) {
      if (j>0) cout << ','; 
      fc = paired_corner(*v0, i, v); 
      cout << fc.vertex; 
    }
    cout << ")\n"; 
  }

  fc_map done; 
  int nf; 
  for (i=0; i<n_fac; ++i) {
    nv = FL[i].size(); 
    for (j=0, v0=FL[i].begin(); j<nv; ++j, ++v0) {
      fc = face_corner(i,j); 
      if (done.find(fc)!=done.end()) continue; 

      cout << '('; 

      // Find adjacency cycle starting at v/fc.
      fc0 = fc; 
      v = *v0; 
      nf = 0; 
      do {
	done[fc] = fc; 
	if (fc!=fc0) cout << ','; 
	cout << '(' << fc.face << ',' << fc.vertex << ')'; 
	if (v.distance() > tube::ideal_cutoff) cout << '*'; 
	fc = adjacent_corner(v, v.mark_p, av);
	v = av; 
      } while (fc!=fc0 && ++nf < 10); 

      cout << ')' << endl; 
    }
  }
}

uc_vertex connected_faces::vertex(face_corner const& fc) const
{
  uc_polygon::const_iterator v = FL[fc.face].begin(); 
  for (int i=0; i < fc.vertex; i++) v++;
  return *v; 
}

// Returns the complex length of an edge of the tube as seen from the
// core geodesic. The edge is the incoming edge of fc. 

Complex connected_faces::edge_holonomy(face_corner fc) const
{
  if (fc.type != fc_back) return Zero;
  if (fc.face < 0) return Zero; 

  face_corner prev(fc);
  prev.vertex--;
  if (prev.vertex < 0) prev.vertex = FL[fc.face].size()-1; 
  return Complex(vertex(fc)) - Complex(vertex(prev));
}

void connected_faces::recover_holonomies(Triangulation* m, Complex H[2]) const
{
  if (!TB) {
    cout << "Triangulation builder must be initialized before call to recover_holonomies.\n";
    return; 
  }

  face_vertex tri; 
  int v, f, i;
  Complex z; 
  Tetrahedron* tet; 

  H[0] = Zero; 
  H[1] = Zero; 

  // Iterate over tetrahedra
  for (tet = m->tet_list_begin.next;
       tet != &m->tet_list_end;
       tet = tet->next) {

    // Iterate over faces
    for (f = 0; f < 4; f++) {

      tri.face = 4*tet->index + f; 

      // Iterate over tet face vertices (corners)
      for (v = 0; v < 4; v++) {
	if (v==f) continue; 

	tri.vertex = v; 

	z = 0.5 * edge_holonomy(TB->dual_edge(tri)); 

	// Add up meridian and longitude contributions
	for (i=0; i<2; i++) {
	  H[i] += double(tet->curve[i][right_handed][v][f]) * z; 
	}
      }
    }
  }
}

bool connected_faces::initialize_triangulation_bldr(int report) 
{
  int n_fac = FL.size(); 
  if (!n_fac) {
    cerr << "Connected faces are required for this function.\n";
    return false; 
  }

  if (TB) { delete TB; }
  TB = new triangulation_builder();

  int i, j, nj; 
  face_corner fc, paired, adjacent, next, rev_adj, paired_adj;
  uc_polygon::const_iterator v, nv; 
  uc_vertex d, pc; 

  for (i=0; i<n_fac; i++) {

    // v and nv point at a corner and the next corner around.
    // j and nj give the indices of these corners. 

    for (v=FL[i].begin(), nv = v, ++nv, j=0, nj=1; 
	 v!=FL[i].end(); 
	 ++v, ++nv, ++j, ++nj) {

      if (nv == FL[i].end()) { nv = FL[i].begin(); nj = 0; }
    
      fc = face_corner(i,j); 
      next = face_corner(i,nj); 
      adjacent = adjacent_corner(*v, v->mark_p, d);
      paired = paired_corner(*v, i, pc); 

      if (v->distance() < tube::ideal_cutoff) {
	if (report)
	  ::print_face_corners(fc,paired,adjacent,next); 

	TB->add_face_corner(fc,paired,adjacent,next); 

      } else {

	rev_adj = adjacent_corner(*v, v->mark_n, d);
	paired_adj = adjacent_corner(pc, pc.mark_p, d); 

	if (report) {
	  ::print_face_corners(fc, paired(fc_front), adjacent(fc_front), fc(fc_front));
	  ::print_face_corners(fc(fc_front), paired, fc(fc_top), next);
	  ::print_face_corners(fc(fc_top), fc(fc_ideal), rev_adj, adjacent(fc_top));
	  ::print_face_corners(fc(fc_ideal),fc(fc_top),paired_adj(fc_ideal),rev_adj(fc_ideal));
	}

	TB->add_face_corner(fc, paired(fc_front), adjacent(fc_front), fc(fc_front));
	TB->add_face_corner(fc(fc_front), paired, fc(fc_top), next);
	TB->add_face_corner(fc(fc_top), fc(fc_ideal), rev_adj, adjacent(fc_top));
	TB->add_face_corner(fc(fc_ideal), fc(fc_top), paired_adj(fc_ideal), rev_adj(fc_ideal));
      }

    }
  }
  return true; 
}

bool connected_faces::get_triangulation_data(TriangulationData& TD, int report)
{
  if (!initialize_triangulation_bldr(report > 1)) return false; 
  TB->compute_cell_complex(); 
  return TB->get_triangulation_data(TD, report); 
}

void connected_faces::print_triangulation_data()
{
  if (!initialize_triangulation_bldr()) return; 
  TB->compute_cell_complex(); 
  TB->print_cell_complex(); 
}


static bool point_in_face(point const& p, uc_polygon const& F)
{
  uc_polygon::const_iterator k; 

  for (k = F.begin(); k != F.end(); k++)
    if (alt_distance(p, *k) < tube_face::boundary_epsilon) return true; 

  return false; 
}

int connected_faces::subface_mark_number(Mark const& mark, point const& p, 
					 vector<int> const& subface_start, 
					 vector<int> const& subface_count) const
{
  int i, fnum = mark.end_index();
  int from = subface_start[fnum]; 

  if (subface_count[fnum/2]==1) 
    return Mark(from, mark).value();

  int to = from + 2*subface_count[fnum/2]; 

  point pim = trans(-mark) * p; 

  for (i=from; i<to; i+=2)
    if (point_in_face(pim, FL[i])) 
      return Mark(i, mark).value();

  cout << "Edge renumbering failed.\n"; 
  return Mark(from, mark).value();
}

Triangulation* connected_faces::drill()
{
  TriangulationData TD;
  if (!get_triangulation_data(TD, false)) {
    cout << "drilling failed\n";
    return 0;
  }
  TD.name = "unknown";

  // Store the new manifold. 
  Triangulation *new_m; 
  data_to_triangulation(&TD, &new_m);
  if (!new_m) {
    cout << "drilling returned invalid triangulation data.\n"; 
  }
  
  // Clean up TD. 
  my_free_array(TD.tetrahedron_data); 
  my_free_array(TD.cusp_data); 
  
  // Determine the refilling curve. 
  Complex H[2]; 
  recover_holonomies(new_m, H); 
  double det = (H[0].real * H[1].imag - H[0].imag * H[1].real);
  double me = -H[1].real * TWO_PI / det;
  double lo =  H[0].real * TWO_PI / det;
  
  // Adjust so that refilling curve is (1,0). If there is a 
  // problem changing it just print out the new refilling curve. 
  int ime = int(floor(me + .5)); 
  int ilo = int(floor(lo + .5)); 
  long int a, b; 
  if (fabs(ime-me) > 1e-2 || fabs(ilo - lo) > 1e-2 || 
      euclidean_algorithm(ime, ilo, &a, &b)!=1) {
    cout << '[' << H[0] << ", " << H[1] << ']' << endl; 
    cout << "re-filling " << '(' << me << ", " << lo << ')' << endl; 
    return new_m; 
  }

  int nc = get_num_cusps(new_m); 
  MatrixInt22 *pcb = new MatrixInt22 [nc];
  pcb[0][0][0] = ime; 
  pcb[0][0][1] = ilo; 
  pcb[0][1][0] = -b; 
  pcb[0][1][1] = a; 
  int i; 
  for (i=1; i<nc; i++) {
    pcb[i][0][0] = 1; 
    pcb[i][0][1] = 0; 
    pcb[i][1][0] = 0; 
    pcb[i][1][1] = 1; 
  }
  if (change_peripheral_curves(new_m, pcb)!=func_OK)
    cout << "re-filling " << '(' << me << ", " << lo << ')' << endl; 

  delete [] pcb;
  
  return new_m;
}

