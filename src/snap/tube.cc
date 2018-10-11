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
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "tube.hh"
#include "snappea/SnapPea.h"
#include "snappea/triangulation.h"
#include "tile.hh"
#include "complex_lll.hh"
#include "orthoangles.hh"
#include "connected_faces.hh"
#include "printable.hh"

using std::setw;
using std::max;
using std::sort;

double tube::ideal_cutoff = 5.;
static double match_epsilon = 5e-5; 

class clipping_cache {
  Complex H[2];
  Complex org; 
  vector<Mark> holo; 
  vector<eqsfunc> surfaces;

public:
  vector<int> transpose_faces; // tweaking

  void set_holonomies(Complex const& m, Complex const& l); 
  int add_face(vector<tube_face>& faces, vector<Ortholine>& ortholines, Ortholine ol, bool report);
  void add_unclipped_face(vector<tube_face>& faces, vector<Ortholine>& ortholines, Ortholine const& ol);
  void print_clippings(Ortholine const& o, int i0) const;

private:
  void add_surfaces(Ortholine const& ol);
  bool transpose_on_face(int n) const;
  Mark marking(int face, int trn) const; 
};

tube::tube()
  : H(2)
{ 
  CF = new connected_faces(this); 
  CC = new clipping_cache; 
}

tube::~tube()
{ 
  delete CF; 
  delete CC;
}

void tube::clear() 
{
  faces.resize(0);
  ortholines.resize(0); 

  delete CF;
  CF = new connected_faces(this); 
  delete CC;
  CC = new clipping_cache; 
}

void tube::set_holonomies(Complex const& me, Complex const& lo)
{
  H[0] = me; 
  H[1] = lo;

  lll_reduce(H);

  CC->set_holonomies(H[0], H[1]);
}

void tube::transpose_face(int n) 
{ 
  CC->transpose_faces.push_back(n); 
}
void tube::transpose_no_faces() 
{ 
  CC->transpose_faces.clear(); 
}

int tube::add_face(Ortholine const& ol, bool report)
{
  return CC->add_face(faces, ortholines, ol, report);
}

void tube::print_clip_list(vector<Ortholine> const& ort) const
{
  cout << H[0] << ' ' << H[1] << endl;

  int i, n = ort.size();

  for (i = 0; i < n; i++) {
    CC->print_clippings(ort[i], i); 
  }
}

void tube::add_unclipped_face(Ortholine const& ol)
{
  CC->add_unclipped_face(faces, ortholines, ol);
}

// CLIPPING CACHE STUFF

static bool zero_on_torus(Complex const& a)
{
  const double eps = 1e-6;
  if (fabs(a.real) > eps) return false; 
  double diff = fabs(a.imag);
  if (diff < eps) return true;
  if (fabs(diff - TWO_PI) < eps) return true;
  return false; 
}

void clipping_cache::set_holonomies(Complex const& mer, Complex const& lon)
{
  H[0] = mer; 
  H[1] = lon; 

  int i,j;
  holo.resize(0);
  holo.reserve(25);
  for (i=-2; i<3; ++i)
    for (j=-2; j<3; ++j)
      holo.push_back(Mark(0,i,j));

  sort(holo.begin(), holo.end());

  int count=0;
  for (j=holo.size()-1; j>=0; --j) {
    for (i=0; i<j; i++)
      if (zero_on_torus(holo[i].holo(H)-holo[j].holo(H))) break;
    if (i<j) {
      holo.erase(holo.begin()+j);
      count++;
    }
  }
#if 0
  cout << "Holonomies [" << H[0] << "," << H[1] << "]\n";
  cout << PSeq(holo) << endl;
  cout << "Num erased = " << count << endl;
#endif
}

void clipping_cache::add_surfaces(Ortholine const& ol)
{
  int i,j,n=holo.size();
  int i0=surfaces.size(); 

  // Make sure there is enough space.
  if (surfaces.capacity()-surfaces.size() < n)
    surfaces.reserve(surfaces.size() + 5*n);

  surfaces.resize(i0 + 2*n); 

  // Append equidistant surface objects for n translates
  // of the two surfaces corresponding to this ortholine. 
  for (i=0; i<2; i++) {
    surfaces[i0+n*i] = eqsfunc(ol.distance(), ol.position(i%2));
    for (j=1; j<n; j++) {
      surfaces[i0+n*i+j] = surfaces[i0+n*i] + holo[j].holo(H);
    }
  }
}

Mark clipping_cache::marking(int face, int trn) const
{
  while (trn >= holo.size()) {
    trn -= holo.size(); 
    face += 1; 
  }

  return Mark(face,0,0) + holo[trn]; 
}

bool clipping_cache::transpose_on_face(int n) const
{
  if (!transpose_faces.size()) return false; 
  vector<int>::const_iterator i; 
  for (i=transpose_faces.begin(); i!=transpose_faces.end(); i++) 
    if (n==*i) break;
  return i!=transpose_faces.end(); 
}

int clipping_cache::add_face(vector<tube_face>& faces, vector<Ortholine>& ortholines, Ortholine ol, bool report)
{
  int i, j, n_ort = ortholines.size(), n = holo.size(); 

  int ns0 = surfaces.size(); 

  // Set up normalization such that faces 0,1 are symmetric wrt 0.
  if (!ns0) {
    Complex d = ol.position(1) - ol.position(0), d1;
    torus_reduce(H,d);
    for (i=0; i<n; i++) {
      d1 = d + holo[i].holo(H); 
      if (complex_modulus_squared(d1) < complex_modulus_squared(d)) 
	d = d1; 
    }
    org = ol.position(0) + 0.5 * d;
    torus_reduce(H,org);
  }

  // Translate the ortholine ends into the torus fundamental domain.
  Complex a = ol.position(0) - org; torus_reduce(H,a);
  Complex b = ol.position(1) - org; torus_reduce(H,b); 
  ol.set_to(ol.distance(), a, b); 

  // Add surfaces for the new faces and their translates. 
  add_surfaces(ol); 

  int this_n; 
  tube_face nf;
  int res;
  bool new_stuff = false;
  Mark clip_marking; 

  // Create two new faces. 
  faces.reserve(faces.size()+2); 
  for (i=0; i<2; i++) {

    this_n = n*i + 2*n*n_ort; 

    faces.push_back(tube_face(surfaces[this_n], marking(faces.size(), 0), false));

    // Clip against existing surfaces. 
    for (j=0; j<surfaces.size(); j++) {

      if (j < ns0 && faces[j/n].is_empty()) {
	// No face exists for this surface so don't bother clipping. 
	j+=n-1; 
	continue; 
      }

      if (j==this_n) continue; 

      clip_marking = marking(j/n, j%n);
      res = faces.back().clip(surfaces[j], clip_marking, false);
      if (res == -1) { 
	cerr << "Was clipping face: " << (faces.size()-1) <<" by "<<clip_marking <<endl;
	return -1; 
      }
      if (faces.back().is_empty()) break;
    }
    if (report) faces.back().print_edge_markings();
    if (!faces.back().is_empty()) new_stuff = true;
  }

  if (!new_stuff) {
    surfaces.resize(ns0); // Get rid of them again. 
    faces.resize(faces.size()-2);
    return 0; 
  }

  ortholines.push_back(ol);

  // Clip existing faces. 
  for (i=0; i<2*n_ort; i++) {

    if (faces[i].is_empty()) continue; 

    for (j=0; j<2*n; j++) {

      clip_marking = marking(2*n_ort, j);
      res = faces[i].clip(surfaces[j+ns0], clip_marking, false); 
      if (res == -1) {
	cerr << "Was clipping face: " << i << " by " << clip_marking << endl;
	return -1; 
      }
      if (res > 0 && report) {
	cout << "Face: " << i << " clipped by (" << clip_marking << ")\n";
      }
    }
  }


  return 1;
}

void clipping_cache::print_clippings(Ortholine const& o, int i0) const
{
  int i, j, n = holo.size(); 

  // Translate the ortholine ends into the torus fundamental domain.
  Complex a = o.position(0); torus_reduce(H,a);
  Complex b = o.position(1); torus_reduce(H,b); 
  Ortholine ol = o;
  ol.set_to(o.distance(), a, b); 

  // Print all clipping surface coordinates for this ortholine. 
  for (i=0; i<2; i++) {
    for (j=0; j<n; j++) {
      cout << ol.distance() << ' ' << (ol.position(i) + holo[j].holo(H));
      cout << marking(2*i0 + i, j) << endl; 
    }
  }
}

void clipping_cache::add_unclipped_face(vector<tube_face>& faces, 
					vector<Ortholine>& ortholines, 
					Ortholine const& ol)
{
  int ns0 = surfaces.size(); 
  add_surfaces(ol); 
  ortholines.push_back(ol);

  // Create two new faces. 
  int i, n = holo.size();
  for (i=0; i<2; i++) {
    faces.push_back(tube_face(surfaces[ns0 + n*i], marking(faces.size(), 0), false));
  }
}

Complex tube::holonomy(Mark const& mark) const
{
  Complex Ho[2];
  Ho[0]=H[0];
  Ho[1]=H[1];

  return mark.holo(Ho);
}

// Compute tube adds a number of faces to the tube at the same time,
// printing a log of faces as they are created and checking face
// pairings until a match is found, or no more ortholines are
// available, or the computation takes too long. The first version
// uses ortholines in os, but not in used, putting ortholines in os
// into used as it goes. The second version just uses the ortholines
// in ort and checks face pairings at the end.

// Return values:
// 0 -OK
// 1 hit max_rad or max_n_radii, and aborted; can continue computation later.
// 2 integrity of tube was lost 
// 3 bad input or problem computing the core geodesic.  
// 4 match no longer possible, unmatched vertex radius less than half current orthodist. 

int tube::compute_tube(Ortholine_set& os, Ortholine_set& used, double rigorous_to_rad, bool report)
{
  // report = 1;
  Ortholine ol; 
  Ortholine_set::const_iterator osi; 
  int res; 
  double urad = 1000.; 
  bool match = false, check = true; 
  time_t start_time = time(0); 
  for (osi = os.begin(); osi != os.end(); osi++) {
    ol = *osi; 

    if (used.find(ol)!=used.end()) continue; 

    if (report) cout << "Ortholine: " << ol << endl;
    res = add_face(ol, report);
    used.insert(ol); 

    if (res==-1) return 2; 
    if (res==1) check = true; 

    if (check) {
      match = check_face_pairings(urad); 
      if (match) break;
      check = false; 
    }

    // Check if tube has any stranded unmatched vertices. 
    if (urad < ol.distance().real/2.0 - 1e-4) { 
      if (urad < rigorous_to_rad/2.0) { 
	cerr << "Unmatched vertex at radius: " << urad << endl; 
	return 4; 
      }
      break; // Will need more ortholines. 
    }
    if (time(0) - start_time > 60.) return 5; 
  }
 
  if (!match) return 1; // Insufficient ortholines. 

  compute_connected_faces();
  return 0; 
}

// result 2 = error
//        1 = incomplete
//        0 = OK

int tube::compute_tube(vector<Ortholine> const& ort, bool report)
{
  int i, n = ort.size();

  for (i = 0; i < n; i++) {
    if (report) cout << "Ortholine: " << ort[i] << endl;
    if (add_face(ort[i], report) == -1) return 2; 
  }

  double urad; // dummy, not used. 
  if (!check_face_pairings(urad)) return 1; 

  if (report) cout << "Tube has " << num_faces() << " faces.\n"; 

  compute_connected_faces();
  return 0; 
}

int tube::compute_tube(Triangulation* T, GroupPresentation* G, bool report)
{
  OrthoAngles OA(T,G);
  if (!OA.valid()) return 3; 

  OA.get_holonomies(H); 
  CC->set_holonomies(H[0], H[1]); 

  EdgeEnd ee;
  int i, n = 2*get_num_tetrahedra(T);
  faces.resize(n); 

  da_spec face;
  vector<da_spec> nbrs; 
  bool match = true; 
  int res; 
  Complex a[2];

  for (ee.first(T); !ee.done(T); ++ee) {
    i = ee.index(); 

    if (i < 0 || i >= n) 
      { cerr << "Index out of range in compute tube\n"; return 2; }

    OA.get_da_spec(ee, face); 
    OA.get_da_spec(ee, nbrs); 

    if (report) cout << "Face " << i << ": "; 
    res = faces[i].make_face(face, nbrs, report);
    if (res < 0) return 2; 
    if (res > 0) match = false; 

    // get the ortholines
    a[i%2] = face.angle(); 
    if (i%2) { // i odd
      ortholines.push_back(Ortholine(face.distance(), CGRay(a[0],-2), CGRay(a[1],-2)));
    }
  }

  if (match) compute_connected_faces();

  return match ? 0:1; 
}

void tube::print_face_matrices(int n) const
{
  if (n < 0 || n > faces.size()) {
    cout << "There is no face numbered " << n << endl;
    return; 
  }
  faces[n].print_clippings(cout); 
}


bool tube::check_face_pairings(double& urad, int report) const
{
  bool match = true, edge_at_infinity = false;

  int i, n_ort = ortholines.size();

  urad = 1000.; 
  double this_urad; 
  for (i=0; i<n_ort; i++) {

    if (report)
      cout << "Faces: " << 2*i << ',' << 2*i+1 << ' ';

    if (!faces_match(faces[2*i], faces[2*i+1], this_urad, match_epsilon, report)) {
      match = false; 
      if (this_urad < urad) urad = this_urad; 
    } else if (report) {
      if (faces[2*i].is_empty()) cout << "are empty.\n";
      else cout << "match.\n";
    }

    if (faces[2*i].has_edge_at_infinity() || faces[2*i+1].has_edge_at_infinity()) {
      edge_at_infinity = true; 
    }
  }

  if (report && match && edge_at_infinity) 
    cout << "Faces match but tube still has edges at infinity.\n"; 

  return match && !edge_at_infinity; 
}

// Find and remove any unmatched short edges. 

bool tube::tidy_edges(bool report)
{
  bool ok = true; 
  Mark mk; 
  eqs_iter e;
  double eps = 1e-5;
  int i, n = faces.size(), j;
  for (i=0; i<n; i++) {
    if (faces[i].is_empty()) continue; 
    mk = faces[i].mark(); 
    for (e = faces[i].begin(); e != faces[i].end();) {
      j = e->mark().end_index(); 
      if (j < 0 || j >= n) {
	if (report) cout << "Edge refers to non-existent face " << j << endl;
	++e; 
	continue; 
      }

      if (!faces[j].has_matching(e, mk, H)) {
	if (ok) {
	  if (report) cout << "Unmatched edges\n";
	  ok = false; 
	}
	if (report) cout << "face " << mk.end_index() << " " << (*e) << endl; 
      }

      if (complex_modulus(e->end(0) - e->end(1)) > eps) {
	++e; 
	continue; 
      }

      faces[i].erase(e);
    }
  }
  return ok; 
}


int tube::num_faces() const
{
  int i, n, n_fac = 2*ortholines.size();

  for (i=0, n=0; i<n_fac; i++) {
    if (!faces[i].is_empty()) n++; 
  }
  return n; 
}

void tube::print_ortholines() const
{
  int i, n_ort = ortholines.size();

  for (i=0; i<n_ort; i++) {
    if (faces[2*i].is_empty()) continue;
    cout << "Faces " << setw(2) << (2*i) << ','
	 << setw(2) << (2*i+1) << ": " << ortholines[i] << '\n';
  }
}

void tube::sort_edges()
{
  // cout << "Sorting tube face edges.\n";
  int i, n_fac = faces.size();
  for (i=0; i<n_fac; i++) {
    if (faces[i].is_empty()) continue; // Skip empty faces. 
    faces[i].sort_edges(); 
  }
}

void tube::print_faces(bool natural) const
{
  int i, n_fac = faces.size();
  for (i=0; i<n_fac; i++) {
    
    if (faces[i].is_empty()) continue; // Skip empty faces. 

    cout << "Face " << i << endl;
    if (natural) cout << faces[i]; 
    else faces[i].uc_print(cout, H); 
  }
}

static void compute_banana_projection_matrix(Complex dist, double radius, R_matrix<3>& M)
{
  // Some intermediate values
  double R2 = cosh(radius)*cosh(radius);
  double SC = sinh(dist.real)*cosh(dist.real); 
  double S2 = sinh(dist.real)*sinh(dist.real); 

  // Center and principal radii for normalized ellipsoid in Klein model
  double x0 = SC/(R2+S2); 
  double ry2 = (R2-1.0)/(R2+S2);
  double rz2 = R2/(R2+S2);
  double rx2 = ry2*rz2;

  // More coeffs for the Klein model ellipsoid equation:
  // (x - x0)^2/rx2 + (c y - s z)^2/ry2 + (s y + c z)^2/rz2 = 1.
  // s2 = s^2 etc. 
  double s2 = sin(dist.imag)*sin(dist.imag); 
  double sc = sin(dist.imag)*cos(dist.imag); 
  double c2 = cos(dist.imag)*cos(dist.imag); 

  // Coefficients for restriction of the ellipsoid equation to the 
  // projection line t -> (t, a t, z). Equation takes the form: 
  // (a0 + a1 a^2) t^2 + (b0 + b1 a z) t + (c0 + c1 z^2) = 0.
  double a0 = 1.0; 
  double a1 = c2*rz2 + s2*ry2;
  double b0 = -2.0*x0;
  double b1 = 2.0*sc*(-rz2 + ry2);
  double c0 = x0*x0 - rx2; 
  double c1 = s2*rz2 + c2*ry2; 

  // Outline of the projection is the set of (a, z) such that 
  // the above equation has repeated zeros, namely the zero set
  // of the discriminant. Coefficients of the discriminant as
  // polynomial in a and z. 
  M(0,0) /* = b0*b0-4.0*a0*c0 */ = 4.0*rx2;  // Constant term. 
  M(1,1) = 2.0*b0*b1;                       // Coeff of a z. 
  M(2,2) /* = b1*b1-4.0*a1*c1 */ = -4.0*rx2; // Coeff of a^2 z^2. 
  M(0,2) = -4.0*a1*c0;                      // Coeff of a^2. 
  M(2,0) = -4.0*a0*c1;                      // Coeff of z^2. 

  M(0,1) = 0.; M(2,1) = 0.; M(1,0) = 0.; M(1,2) = 0.; 
}

static Complex proj_to_cyl(Complex const& z)
{ return Complex(atanh(z.imag), -atan(z.real)); }

static void get_cylinder_projection(list<Complex>& curve, Complex const& dist, double radius)
{
  curve.clear(); 

  R_matrix<3> mx;
  compute_banana_projection_matrix(dist, radius, mx);
  bq_curve M(mx,Mark(),0.); 

  interval_list L;
  L.make_intervals(M);
  if (!L.L.size()) {
    printf("Empty cylinder projection computed!\n");
    return; 
  }

  list<Complex> segment; 
  list<eqs_interval>::const_iterator k;

  list<Complex>::iterator it;
  Complex loop, vx; 
  k=L.L.begin();

  loop = k->end(0);
  while (true) {
    segment = k->get_polyline(); 
    curve.splice(curve.end(), segment); 
    
    vx = k->end(1); 
    k++;
    if (k==L.L.end() || !complex_close(vx, k->end(0), eqs_interval::epsilon)) break; 
  }

  for (it = curve.begin(); it!=curve.end(); it++) 
    *it = proj_to_cyl(*it); 

  if (k!=L.L.end()) {
    printf("Cylinder projection has multiple components!\n"); 
  }
}

static bool contains(const string& s, const string& sub)
{
  return s.find(sub)!=string::npos; 
}

void tube::save_picture(const string& print_options, picfile& pic) const
{
  double diam = max(complex_modulus(H[0]+H[1]), complex_modulus(H[0]-H[1]));

  // Page size is 8.5 by 5.5 inches. One inch is diam/4. 
  // We then add a margin of diam, so that every tile which
  // intersects with the page will be represented. 
  vector<half_plane> page(4); 
  page[0] = half_plane(One, -1.5 * diam * One);
  page[1] = half_plane(-One, 1.5 * diam * One);
  page[2] = half_plane(I, -2.1 * diam * I);
  page[3] = half_plane(-I, 2.1 * diam * I);

  // Find a set of tiles (translations of the fundamental
  // parallelogram) which cover the page. 
  tiler a_tiler = tiler(H[0], H[1]);
  list<Complex> translations = a_tiler.tile(Zero, 3.5 * diam); 
  translations = filter(translations, page); 

  list<Complex>::const_iterator lci; 

  int i, n = 2*ortholines.size(); 

  // Count non-empty faces. We use this count in coloring. 
  int n_fac = num_faces(), fnum = 0; 
  
  if (contains(print_options,"-f")) { // tube domain faces
    color col; 

    for (i=0; i<n; i++) {
      if (faces[i].is_empty()) continue; 

      col = hsbcolor(double(fnum/2)/double(n_fac/2));
      
      if (contains(print_options,"-fbp")) // filled
	faces[i].picfile_print(pic, col);
      if (contains(print_options,"-fbo")) // outline
	faces[i].picfile_print(pic, black); 
	
      if (contains(print_options,"-ftp")) // filled
	faces[i].picfile_print(pic, col, translations);
      if (contains(print_options,"-fto")) // outline
	faces[i].picfile_print(pic, black, translations); 

      ++fnum;
    }
  }
  
  list<Complex> box, big_box; 
  if (contains(print_options,"-b")) { 
    box.push_back(0.5 * (-H[0] - H[1]));
    box.push_back(0.5 * ( H[0] - H[1]));
    box.push_back(0.5 * ( H[0] + H[1]));
    box.push_back(0.5 * (-H[0] + H[1]));
    box.push_back(0.5 * (-H[0] - H[1]));
    pic.print_line(box); 
  }

  if (contains(print_options,"-B")) { 
    big_box.push_back(Complex(-4.0, -PI));
    big_box.push_back(Complex(4.0, -PI));
    big_box.push_back(Complex(4.0, PI));
    big_box.push_back(Complex(-4.0, PI));
    big_box.push_back(Complex(-4.0, -PI));
    pic.print_line(big_box); 
  }
  
  if (contains(print_options,"-o")) {   
    for (i=0; i<n; i++)
      pic.print_point(ortholines[i/2].position(i%2));
  }

  if (contains(print_options,"-t")) {

    // Compute the injectivity radius. 
    int n_ort = ortholines.size();
    double i_rad = 100., r; 
    for (i=0; i<n_ort; i++) {
      r = ortholines[i].distance().real/2.;
      if (i_rad > r) i_rad = r; 
    }

    // Print the cylinder projections. 
    Complex a, b; 
    color col; 
    for (i=0; i<n_ort; i++) {
      list<Complex> outline;
      get_cylinder_projection(outline, ortholines[i].distance(), i_rad); 

      col = hsbcolor(double(i)/double(n_ort),.3,1.); 
      a = ortholines[i].position(0);
      b = ortholines[i].position(1);

      if (contains(print_options,"-tbo")) {
	pic.print_line(outline, a);
	pic.print_line(outline, b);
      }
      if (contains(print_options,"-tbp")) {
	pic.print_polygon(outline, col, true, a);
	pic.print_polygon(outline, col, true, b);
      }
      if (contains(print_options,"-tt")) {
	for (lci=translations.begin(); lci!=translations.end(); lci++) {
	  if (contains(print_options,"-tto")) {
	    pic.print_line(outline, a + *lci);
	    pic.print_line(outline, b + *lci);
	  }
	  if (contains(print_options,"-ttp")) {
	    pic.print_polygon(outline, col, true, a + *lci);
	    pic.print_polygon(outline, col, true, b + *lci);
	  }
	}
      }
    }
  }
}

// We return an array of Complexes each of which gives a 
// point on the page of shape width by height. The points
// are so arranged as to keep them as far apart as possible
// while keeping them on a square grid. The size of the 
// grid spacing is returned. We keep the number of columns
// even: the natural picture showing face pairings looks 
// better that way. 

double arrange_on_page(double width, double height, int num, vector<Complex>& pos)
{
  pos.resize(num); 
  if (num==0) return 0.0; 

  int h = 1, w = 2; 
  while (h*w < num) {
    if (height/h < width/w) w+=2; 
    else h++;
  }

  double size, vmarg, hmarg;
  
  if (height/h < width/w) {

    size = height/h;
    vmarg = 0.0; 
    hmarg = 0.5 * (width - w * size); 

  } else {

    size = width/w; 
    vmarg = 0.5 * (height - h * size); 
    hmarg = 0.0; 

  }

  // Can fit w by h squares of side size on your page.
  // Horizontal margin size hmarg, vertical margin size vmarg.

  int i, j, n=0;
  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {

      pos[n] = Complex(
	     hmarg + size/2.0 + double(j) * size, 
	     height - vmarg - size/2.0 - double(i) * size);

      if (++n==num) return size;
    }
  }

  return size; // But should never get here. 
}

void tube::save_natural_face_picture(picfile& pic) const
{
  int i, n = faces.size();

  // Figure out some page positioning info. 
  vector<Complex> pg_pos;
  double size = arrange_on_page(5.5, 8.5, n, pg_pos); 
  Complex pg_center(2.25, 4.25); 

  for (i=0; i<n; i++)
    faces[i].picfile_print(pic, pg_pos[i]-pg_center, size); 
}

// Forwarding functions to connected_faces object. 

bool tube::compute_connected_faces()
{ 
  CF->set_pairings(ortholines); 
  // tidy_edges(false); 
  if (!CF->compute_connected_faces(faces)) return false; 
  // CF->orient_connected_faces(); 
  return true; 
}

void tube::print_connected_faces() const
{ CF->print_connected_faces(); }
void tube::print_face_corners() const
{ CF->print_face_corners(); }

Triangulation* tube::drill()
{ return CF->drill(); }
void tube::print_triangulation_data()
{ CF->print_triangulation_data(); }


// 0 - tubal
// 1 - non-tubal
// 2 - error

int triangulation_is_tubal(Triangulation* T, GroupPresentation* G)
{
  OrthoAngles OA(T,G);
  if (!OA.valid()) return 2; 

  EdgeEnd ee;
  int i;

  da_spec face;
  vector<da_spec> nbrs; 
  bool match = true; 
  int res; 
  bool report = true;

  for (ee.first(T); !ee.done(T); ++ee) {
    i = ee.index(); 
    if (!(i%2)) continue; 

    OA.get_da_spec(ee, face); 
    OA.get_da_spec(ee, nbrs); 

    { 
      tube_face F;

      cout << "Edge " << (i/2) << ": "; 
      res = F.make_face(face, nbrs, report);
      if (res < 0) return 2; 
      if (res > 0) match = false; 
    }
  }

  return match ? 0:1; 
}

#if 0
bool triangulation_is_tubal(Triangulation* manifold, GroupPresentation* G)
{
  // Get all the edge orthodistances
  vector<Complex> orth;
  if (!get_edge_orthodistances(manifold, G, orth)) {
    return false; 
  }
  if (orth[0]==Infinity) {
    cerr << "Sorry, can't compute orthoangles for complete structure.\n";
    return false;
  }

  OrthoAngles OA(manifold, G); 
  if (!OA.valid()) return false;

  EdgeClass *e; 
  vector<Complex> da_spec; 

  tube_face* tf; 
  bool ok=true; 
  Complex d0; 

  // Run through the edges
  for (e = manifold->edge_list_begin.next; 
       e!= &manifold->edge_list_end;
       e = e->next) {

    cout << "Edge: " << e->index << ' '; 

    // Get orthodist/ortho-angle spec around this edge. 
    OA.get_da_spec(e, da_spec); 

    // Compute the corresponding face and check it. 
    d0 = complex_acosh(orth[e->index]); 
    if (d0.real < 0.) d0 = -d0; 
    if (!verify_tube_face(d0, da_spec, tf)) ok=false;
    delete tf; 
  }

  return ok; 
}
#endif

