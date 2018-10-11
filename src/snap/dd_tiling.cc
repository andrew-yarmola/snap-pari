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
#include "dd_tiling.hh"
#include "snappea/length_spectrum.h"
#include "kernel_extras.hh"
#include "closed_geodesics.hh"
#include "gv_print_Dirichlet.hh"

#include <cstring>

using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::strcmp;

class odist_cache {
  const WEPolyhedron* domain;
  const vector<GeodesicWord>* GL;
  double odist; 
  int max_rq;
  list<interval> ivls; 
  Ortholine_set os; 
  static double max_rad;
public:
  int refs; 

  odist_cache(const WEPolyhedron* d, const vector<GeodesicWord>* gw, int g);

  // double rigorous() const { return tiling->radius - offset; }
  int max_requested() const { return max_rq; }
  int get_orth(int n, Ortholine& ol); 
};

double odist_cache::max_rad = 6.2; 

odist_cache::odist_cache(const WEPolyhedron* d, const vector<GeodesicWord>* gw, int g) 
  : domain(d), GL(gw), odist(.5), max_rq(-1), refs(1) 
{
  get_crossing_lifts(domain, GSpec(g, &(*gw)[g]), ivls, true);
}

static Ortholine const& get(Ortholine_set const& os, int n)
{
  Ortholine_set::const_iterator it=os.begin(); 
  int i; 
  for (i=0; i<n; i++) ++it; 
  return *it;
}

// result = 0, no more ortholines (null returned)
// result = 1, non-rigorous next ortholine returned (no longer supported).
// result = 2, rigorous ortholine returned. 

int odist_cache::get_orth(int n, Ortholine& ol)
{
  while (os.size() <= n) {
    odist += 0.5;

    if (odist > max_rad) break; 

    new_ortholines(domain, ivls, *GL, odist, os); 
  }


  if (n > max_rq) max_rq = n; 

  if (os.size() < n) {
    return 0;
  }
  
  ol = get(os,n); 
  return 2; 
}

geo_orbit::geo_orbit(const WEPolyhedron* domain, const vector<GeodesicWord> *gw, int gn)
{
  OC = new odist_cache(domain, gw, gn); 
}

geo_orbit::~geo_orbit()
{
  if (OC && !(--OC->refs)) delete OC;
}

geo_orbit::geo_orbit(geo_orbit const& go)
: OC(go.OC), conjugate(go.conjugate), length(go.length), gnum(go.gnum)
{
  if (OC) ++(OC->refs);
}

geo_orbit& geo_orbit::operator = (geo_orbit const& go)
{
  if (this == &go) return *this; 
  if (OC && !(--OC->refs)) delete OC;
  OC = go.OC; 
  if (OC) ++(OC->refs);
  conjugate = go.conjugate;
  length = go.length;
  gnum = go.gnum; 

  return *this; 
}

int geo_orbit::get_orth(int n, Ortholine& o) const
{
  if (!OC) return 0; 
  return OC->get_orth(n, o); 
}

vector<Ortholine> geodesic_ortholines(const WEPolyhedron* domain,
	 const vector<GeodesicWord> *gw, int gn, int n) 
{
  odist_cache OC(domain, gw, gn);

  Ortholine ol; 
  vector<Ortholine> ort; 
  int i, ok=2;
  for (i=0; i<n; i++) {
    if (OC.get_orth(i,ol) < 2) {
      cout << "max radius reached, out of ortholines\n";
      break; 
    }
    ort.push_back(ol);
  }
  return ort; 
}

void print_Dirichlet(const WEPolyhedron* domain, const GroupPresentation* g, 
		     FILE* fp, unsigned int opts)
{
  domain->conjugacy.print(fp); 
  fprintf(fp, "\n"); 
  print_polyhedron(domain, g, opts, fp);
}

// Saves the face pairing transformations for an existing Dirichlet domain
// and a set of symmetries. All are saved in the Dirichlet domain coordinate 
// system. 

void save_domain_and_symmetries(const WEPolyhedron* domain, vector<O31_matrix> const& symms, FILE* out)
{
  int n = domain->num_faces + symms.size()-1; 

  // Write header for output file.
  fprintf(out, "%% Generators\n%d\n\n", n); 

  WEFace *face;
  for (face = domain->face_list_begin.next;
       face != &domain->face_list_end;
       face = face->next)
    {
      snappea_print_o31(*face->group_element, out); 
    }
  int i; 
  n = symms.size(); 
  for (i=1; i<n; i++) {
    snappea_print_o31(symms[i], out);
  }
}


void convert_Dirichlet(const GroupPresentation* g, FILE* in, FILE* out)
{
  char buf[1024];
  int c; 

  // Skip comments at start. 
  while ((c = getc(in)) == '#' || c == '\n') {
    if (c == '#') fgets(buf, 1023, in);
  }
  ungetc(c, in); 

  // Read the conjugacy matrix. 
  O31_matrix conj; 
  if (!conj.read(in)) {
    printf("Unable to read conjugacy matrix from dirichlet file.\n"); 
    return; 
  }

  // Read words. 
  list<FGWord> words; 
  while (fgets(buf, 1023, in)) {
    if (buf[0]=='#') continue; 
    if (strcmp(buf, ".\n")==0) break;
    for (c=0; buf[c]!='\0'; c++) {
      if (buf[c]=='\n') buf[c]='\0'; // kill '\n'. 
    }
    words.push_back(FGWord(buf));
  }

  // Write header for output file.
  fprintf(out, "%% Generators\n%u\n\n", words.size()); 

  // Write matrices.
  list<FGWord>::const_iterator i; 
  O31_matrix mx; 
  for (i=words.begin(); i!=words.end(); i++) {
    word_to_matrix(g, *i, mx, conj, 1);
    snappea_print_o31(mx, out); 
  }
}

bool group_contains(const WEPolyhedron* H, O31_matrix g, double eps, FGWord& wd)
{
  O31_vector pt(g * O31_origin); 
  wd = FGWord(); // clear wd.
  dirichlet_normalize(H, pt, g, wd, 0); 
  return g.distance_from_identity() < eps; 
}

bool group_contains(const WEPolyhedron* H, O31_matrix g, double eps)
{
  O31_vector pt(g * O31_origin); 
  FGWord wd; // dummy.
  dirichlet_normalize(H, pt, g, wd, 0); 
  return g.distance_from_identity() < eps; 
}

// Check if G is a subgroup of H after conjugation by c. 
//
// i.e. c.G.c^-1 < H, or equivalently, G < c^-1.H.c.

bool is_subgroup(const GroupPresentation* G, const WEPolyhedron* H, 
		 O31_matrix const& c, double eps, bool report=false)
{
  O31_matrix c_inv(inverse(c)); 
  FGWord w;
  int i, n = fg_get_num_generators(G); 
  for (i=0; i<n; i++) {
    if (!group_contains(H, c * inverse(fg_gen(G,i)) * c_inv, eps, w)) break; 
    if (report) cout << tochar(i+1) << " -> " << w << endl; 
  }
  return i==n; 
}

bool is_symmetry(const WEPolyhedron* domain, const GroupPresentation* G, 
		 O31_matrix const& mx, double eps, bool report)
{
  return is_subgroup(G, domain, mx * inverse(domain->conjugacy), eps, report);
}

// Tiling code.

dd_tiling::~dd_tiling()
{
  if (root) free_tiling(root);
}

void dd_tiling::tile(const WEPolyhedron* domain, double tile_rad, int report)
{
  if (tile_rad <= radius) return; 
  if (!root) {
    if (report) cout << "Creating tiling.\n";
    ::tile(domain, tile_rad, &root);
  } else {
    if (report) cout << "Expanding tiling.\n";
    expand_tiling(domain, tile_rad, root);
  }
  if (report) 
    cout << "Tiling has " << count_translates(root) << " tiles.\n"; 
  radius = tile_rad; 
}

void dd_tiling::clear()
{
  if (root) {
    free_tiling(root);
    root = 0; 
  }
  radius = 0.0;
}

int dd_tiling::num_tiles() const
{
  return count_translates((Tile*)root);
}

bool compute_geodesics(const WEPolyhedron* domain, dd_tiling& tiling, 
		       double max_length, vector<GeodesicWord>& geodesics)
{
  /* Tile out sufficiently far to compute all geodesics up to max_length. */
  tiling.tile(domain, 2.0 * acosh(cosh(domain->spine_radius) * cosh(max_length/2.0))); 

  int i, num_geodesics;
  GeodesicWord* the_geodesics = 0;
  list_geodesics(tiling.root, domain->spine_radius, max_length, &the_geodesics, 
		 &num_geodesics, TRUE); 

  // Copy them into the vector. 
  geodesics.resize(0); 
  geodesics.reserve(num_geodesics); 
  for (i=0; i<num_geodesics; i++) 
    geodesics.push_back(the_geodesics[i]); 

  my_free_array(the_geodesics); 

  sort(geodesics.begin(), geodesics.end()); 
  return num_geodesics!=0; 
}

bool compute_geodesics(const WEPolyhedron* domain, dd_tiling& tiling, 
		       int min_n, vector<GeodesicWord>& geodesics)
{
  int prev_num = geodesics.size();
  if (prev_num > min_n) return true; 

  double current_cutoff; 
  if (prev_num) current_cutoff = geodesics.back().length.real + .5;
  else current_cutoff = 2.0; 

  int guard = 5; 
  GeodesicWord* the_geodesics = 0;
  int num_geodesics;
  while (--guard) {

    /* Tile out sufficiently far to compute all geodesics up to current_cutoff. */
    tiling.tile(domain, 2.0 * acosh(cosh(domain->spine_radius) * cosh(current_cutoff/2.0)));

    list_geodesics(tiling.root, domain->spine_radius, 
		   current_cutoff, &the_geodesics, &num_geodesics, TRUE); 

    if (num_geodesics >= min_n) break;
    if (num_geodesics > prev_num) guard = 5; // We're getting somewhere. 
    prev_num = num_geodesics; 
    current_cutoff += 0.3; // Tile a bit further next time. 
  }

  int i; 
  if (num_geodesics >= min_n) {
      // Copy them into the vector. 
      geodesics.resize(0); 
      geodesics.reserve(num_geodesics); 
      for (i=0; i<num_geodesics; i++) 
	geodesics.push_back(the_geodesics[i]); 
  }

  my_free_array(the_geodesics); 
  sort(geodesics.begin(), geodesics.end()); 

  return num_geodesics >= min_n; 
}

static void same_length_range(vector<GeodesicWord> const& geodesics, int gn, int& b, int& e)
{
  Complex uhp_len, this_len; 
  uhp_len = geodesics[gn].length; 
  uhp_len.imag = fabs(uhp_len.imag);

  b = gn;
  while (true) {
    --b;
    if (b < 0) break; 
    this_len = geodesics[b].length; 
    if (fabs(uhp_len.real-this_len.real) > 1.0e-6) break; 
    if (fabs(uhp_len.imag-fabs(this_len.imag)) > 1.0e-6) break; 
  }
  ++b; 
  e = gn; 
  while (true) {
    ++e;
    if (e == geodesics.size()) break; 
    this_len = geodesics[e].length; 
    if (fabs(uhp_len.real-this_len.real) > 1.0e-6) break; 
    if (fabs(uhp_len.imag-fabs(this_len.imag)) > 1.0e-6) break; 
  }
}


// MoebiusTransformation mt = word_to_Moebius(group, wd);
// ie. mt is expected in fundamental group coordinates rather than
// Dirichlet domain coordinates.

void print_element_info(const WEPolyhedron* domain, 
		     vector<GeodesicWord> const& geodesics, 
			MoebiusTransformation const& mt, bool dd_coords) 
{
  // cout << mt << endl; 
  line l;
  int nfix = fixed_points(mt.matrix, l);
  if (nfix==2) {
    Complex clen = complex_length_mt(mt); 
    cout << "Complex length: " << clen << endl; 
    cout << "Axis: " << l << endl; 

    // Get line in Dirichlet domain coordinates. 
    MoebiusTransformation M = inverse(domain->conjugacy); 
    O31_line L = O31_line(dd_coords ? l : M*l); 

    // Print distance from Dirichlet domain basepoint. 
    cout << "Distance from Dirichlet basepoint: " << distance_to_org(L) << endl; 

    FGWord conj; 
    int ori, gn; 
    if (!find_geodesic(domain, geodesics, L, 0, gn, ori, conj)) {
      cout << "No matching geodesic was found.\n"; 
      return; 
    } 

    // Word matches geodesic gn, conj = conjugacy, ori = orientation.
    // Print out all the relevant details. 
    double mult = clen.real/geodesics[gn].length.real; 
    double imult = floor(mult + 0.5); 
    if (fabs(mult-imult)>.0001)
      cout << "Bad multiplicity: " << mult << endl; 
    else ori *= int(imult); 
    cout << "Geodesic: " << gn << " Conjugacy: ";
    if (conj.length()) cout << conj;
    else cout << "."; 
    cout << " Multiplicity: " << ori << endl; 
  } else if (nfix==1) {
    cout << "Parabolic fixed point: " << l.end[0] << endl; 
  } else if (nfix==-1) {
    cout << "Identity\n"; 
  }
}

void get_crossing_lifts(const WEPolyhedron* domain,
			GSpec const& g, 
			list<interval>& ivls, bool full_set) 
{
  line l; 
  if (fixed_points(g.g->mt.matrix, l) != 2) {
    cout << "Geodesic " << g.index << " in get_crossing_lifts does not correspond to a hyperbolic element!\n";
    return;
  }
  
  all_crossing_lifts(domain, l, g, ivls, full_set);
}

// ivls here should not be a full set: want a set of intervals covering 
// each closed geodesic once and once only. 

void compute_ortholines(const WEPolyhedron* domain, dd_tiling& tiling, 
			list<interval> const& ivls, double radius, 
			Ortholine_set& os, bool check_range) 
{
  if (!ivls.size()) return; 

  // Ensure tiling is big enough. 
  tiling.tile(domain, 2.0 * outradius(ivls) + radius); 

  // Get the ortholines. 
  get_ortholines(tiling.root, ivls, 0., tiling.radius, os, radius, check_range);
}

void compute_ortholines(const WEPolyhedron* domain,
			dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
			vector<int> const& vgn, double radius, Ortholine_set& os)
{
  // Get the intervals corresponding to this set of geodesics.
  int i, gn; 
  list<interval> ivls;
  for (i=0; i<vgn.size(); i++) {
    gn = vgn[i];
    get_crossing_lifts(domain, GSpec(gn, &geodesics[gn]), ivls, false);
  }
  compute_ortholines(domain, tiling, ivls, radius, os);
}

vector<Ortholine> copy_ortholines(Ortholine_set const& OS, double cutoff)
{
  int i, n = 0; 
  Ortholine_set::const_iterator osi; 

  // Count them. 
  for (osi = OS.begin(); osi != OS.end(); osi++) {
    if ((*osi).distance().real > cutoff) break; 
    n++; 
  }

  // Copy them. 
  vector<Ortholine> ort(n);
  for (i = 0, osi = OS.begin(); i < n && osi != OS.end(); osi++, i++) ort[i] = *osi;

  return ort; 
}

#if 0
// don't want a full set of intervals here. 

static bool get_first_ortholine(const WEPolyhedron* domain, dd_tiling& tiling, 
				list<interval> const& ivls, 
				double& orad, double max_orad, Ortholine_set& os, 
				Ortholine_set::const_iterator& it)
{
  double offset = 2.0 * outradius(ivls);
  double prev_orad;

  tiling.tile(domain, orad+offset); 
  get_ortholines(tiling.root, ivls, 0.0, orad+offset, os, max_orad); 
  if (os.size() > 0 && (*os.begin()).distance().real < orad) {
    it = os.begin(); 
    return true; 
  }

  // Make sure there is at least one ortholine. 
  while (orad < max_orad) {
    prev_orad = orad;
    orad += 0.5; 
    if (orad > max_orad) orad = max_orad + 0.01; 

    tiling.tile(domain, orad+offset); 
    get_ortholines(tiling.root, ivls, prev_orad+offset, orad+offset, os, max_orad); 

    if (os.size() > 0 && (*os.begin()).distance().real < orad) {
      it = os.begin(); 
      return true; 
    }
  }

  return false;
}

static bool get_next_ortholine(const WEPolyhedron* domain, dd_tiling& tiling, 
			       list<interval> const& ivls, 
			       double& orad, double max_orad, Ortholine_set& os, 
			       Ortholine_set::const_iterator& it)
{
  Ortholine prev = *it; 

  // Check if there is a next ortholine to return. 
  it++; 
  if (it != os.end() && (*it).distance().real < orad) return true; 

  double prev_orad, offset = 2.0 * outradius(ivls); 
  int old_size; 

  while (orad < max_orad) {

    // Expand our collection of ortholines. 
    old_size = os.size();
    prev_orad = orad;
    orad += 0.5; 
    if (orad > max_orad) orad = max_orad + 0.01; 
    
    tiling.tile(domain, orad+offset); 
    get_ortholines(tiling.root, ivls, prev_orad+offset, orad+offset, os, max_orad); 

    // Nothing new this time around. 
    if (old_size == os.size()) continue; 

    // Find our position again. 
    it = os.begin();
    while (it != os.end() && !((*it) == prev)) it++; 
    if (it == os.end()) {
      cout << "lost our place in get_next_ortholine\n";
      return false; 
    }

    ++it; // Go past our previous position. 
    if (it == os.end()) { // Didn't get any new ones. 
      cout << "ortholine has appeared below presumed rigorous radius!\n";
      return false; 
    }

    if ((*it).distance().real < orad) return true; 
  }

  return false;
}

#endif 

/* This function returns a new set of equal real orthodistance
   ortholines each time it is called, incrementally expanding its
   ortholine set whenever necessary. */

static bool next_ortholine_group(const WEPolyhedron* domain, dd_tiling& tiling, 
			     list<interval> const& ivls, 
			     double max_radius, Ortholine_set& os, 
			     Ortholine_set::const_iterator& it)
{
  double radius; 
  int prev_size = os.size(); 

  if (prev_size) {
    radius = it->distance().real; 
    // Skip to the end of the current group of ortholines of equal length. 
    while (it != os.end() && it->distance().real < radius + 1e-4) it++;
    if (it != os.end()) return true; 
  } else {
    radius = 0.5; 
  }

  // Get some new ortholines. 
  while (1) {
    radius += 0.5; 
    if (radius > max_radius) radius = max_radius;
    compute_ortholines(domain, tiling, ivls, radius, os);
    if (os.size() > prev_size || radius == max_radius) break; 
  }

  if (os.size() <= prev_size) {
    it = os.end();
    return false; // Failed to find any new ortholines. 
  } 

  it = os.begin(); 
  int i; 
  for (i=0; i<prev_size; i++) it++; // Skip to end of what we had already.
  return true; 
}

static bool next_geodesic_group(const WEPolyhedron* domain, 
       			      dd_tiling& tiling, vector<GeodesicWord>& geodesics, int& g_num,
				list<interval>& ivls, double max_length)
{
  double length;
  int i, next_gn; 

  // Deal with all cases in which we already have enough geodesics. 
  if (g_num >= 0) {
    length = geodesics[g_num].length.real;
    same_length_range(geodesics, g_num, i, next_gn);
    g_num = next_gn;
    if (g_num < geodesics.size())
      goto all_OK; 
  } else {
    g_num = 0; 
    if (geodesics.size()) {
      goto all_OK; 
    }
    length = 0.5; 
  }

  // At this point g_num = geodesics.size(). 
  // Get some new geodesics. 
  while (1) {
    length += 0.5; 
    if (length > max_length) length = max_length;
    compute_geodesics(domain, tiling, length, geodesics);
    if (geodesics.size() > g_num || length == max_length) break; 
  }

  if (geodesics.size() <= g_num) return false;

 all_OK:

  // Get intervals for this set of geodesics. 
  ivls.clear(); 
  same_length_range(geodesics, g_num, i, next_gn);
  for (; i<next_gn; i++) 
    get_crossing_lifts(domain, GSpec(i, &geodesics[i]), ivls, false);
  if (!ivls.size()) {
    cout << "Get_crossing_lifts failed for geodesics " << g_num;
    cout << " to " <<(next_gn-1)<< endl;
    return false; 
  }
  return true;
}

/* A go pair is a pair consisting of a set of equal length geodesics
   and a set of equal real orthodistance ortholines between them. The
   following function returns go pairs in lexicographic order of
   (geodesic_length, real_orthodistance) within a given (max_length,
   max_distance) bound. We skip geodesics of length < 0.3 and
   ortholines of orthodistance 0. */

bool next_go_pair(const WEPolyhedron* domain, 
			 dd_tiling& tiling, vector<GeodesicWord>& geodesics, int& g_num,
			 list<interval>& ivls, double max_len, 
			 Ortholine_set& os, Ortholine_set::const_iterator& osi,
		  double max_dist)
{
  if (g_num < 0) {
    double length = 0.0;
    while (length < 0.3) {
      if (!next_geodesic_group(domain, tiling, geodesics, g_num, ivls, max_len))
	return false; 
      length = geodesics[g_num].length.real; 
    }
  }

  double dist;
  while (true) {
    dist = 0.0;
    while (dist < 1e-4) {
      if (!next_ortholine_group(domain, tiling, ivls, max_dist, os, osi)) break; 
      dist = osi->distance().real;
    }
    if (dist >= 1e-4) break; 
    
    if (!next_geodesic_group(domain, tiling, geodesics, g_num, ivls, max_len)) 
      return false; 
    os.clear(); 
  }

#if 0
  // debugging
  cout << '(' << g_num;
  Ortholine_set::const_iterator it;
  for (it = osi; it!=os.end() && it->distance().real < dist+1e-4; ++it)
    cout << ',' << it->distance();
  cout << ')' << endl;
#endif

  return true; 
}

void print_go_pairs(const WEPolyhedron* domain, double max_len, double max_dist)
{
  dd_tiling tiling; 
  vector<GeodesicWord> geodesics;
  int g_num = -1; 
  list<interval> ivls; 
  Ortholine_set os; 
  Ortholine_set::const_iterator osi;

  while (next_go_pair(domain, tiling, geodesics, g_num, ivls, max_len, os, osi, max_dist)) {
    cout << '(' << geodesics[g_num].length.real << ',' << (*osi).distance().real << ')' << endl;
  }
}




void negate_row(O31_matrix& m, int r)
{
  int c;
  for (c = 0; c < 4; c++) {
    m(r,c) = -m(r,c); 
  }
}

// A T-matrix maps (x-axis, y+-axis) onto (ivl, ortholine-seg).
// The list of intervals is assumed to contain one at end i of ol. 

O31_matrix get_T_matrix(Ortholine const& ol, int i, list<interval> const& ivls)
{
  Complex olp = ol.position(i); 
  int gnum = ol.geodesic_num(i); 

  // Find the interval with the right geodesic number which contains
  // end i of the ortholine. 

  list<interval>::const_iterator j; 
  for (j=ivls.begin(); j!=ivls.end(); j++) {
    if (j->gnum() != gnum) continue; 
    if (j->in_range(olp.real)) break; 
  }

  if (j==ivls.end()) {
    cout << "Unable to find an interval containing end " << i << " of\n";
    cout << "ortholine: " << ol; 
    return O31_matrix(1); // return identity. 
  }

  return j->the_line().T() * O31_x_trans(olp);
}

static double mod_pi(double x)
{
  while (x >= PI) x-= PI;
  while (x < 0.) x+= PI; 
  return x; 
}

static double quarter_fold(double x)
{
  while (x >= PI) x-= PI;
  while (x < 0.) x+= PI; 
  if (x > PI_OVER_2) x = PI - x; 
  return x; 
}

// Goes through ortholines of equal or conjugate orthodistance to osi. 
// For each it constructs a set of candidate symmetries. 
// It then tests each candidate symmetry to see if it is an actual symmetry.
// Symmetry of what? 
// If group is non-null, check if it is a symmetry of the quotient space.
// If group is null, check if it is a symmetry of the lifts of the geodesics
// indicated by ivls. 

static int symmetries_from_ort(const WEPolyhedron* domain, 
			       const GroupPresentation* group, 
		       list<interval> const& ivls, 
		       Ortholine_set const& os, 
			       Ortholine_set::const_iterator osi, 
		       vector<O31_matrix>& matrices, int report)
{
  matrices.resize(0);
  Complex od = osi->distance(); 
  double qf_arg = quarter_fold(od.imag); 

  bool symmetric_od = fabs(od.imag) < 1e-6 || 
    fabs(fabs(od.imag) - PI_OVER_2) < 1e-6 ||
    fabs(fabs(od.imag) - PI) < 1e-6; 

  if (report) cout << "Using ortholines with orthodistance: " << od << endl; 
  if (symmetric_od && report) {
    cout << "These ortholines have extra reflection symmetries.\n"; 
  }

  O31_matrix T, S; 
  O31_matrix T0[2][2]; 

  T0[0][0] = inverse(get_T_matrix(*osi, 0, ivls)); 
  T0[0][1] = T0[0][0]; negate_row(T0[0][1], 3); // T0[0][1] reverses orientation
  T0[1][0] = T0[0][1]; negate_row(T0[1][0], 1); // T0[1][0] reverses geodesic
  T0[1][1] = T0[0][0]; negate_row(T0[1][1], 1); // T0[1][1] reverses geodesic and ori

  int ori, gdir, oln = 0; 

  vector<O31_vector> T_points; 
  
  if (report) cout << "Symmetries:\n"; 

  int i, j; 
  bool sym;

  // Iterate over ortholines of equal length. 
  for (; osi != os.end() && osi->distance().real < od.real + 1e-4; 
       osi++, oln++) {

    // Skip any with complex orthodistance neither equal nor conjugate mod pi. 
    if (fabs(quarter_fold(osi->distance().imag) - qf_arg) > 1e-4) {
      if (report) cout << "Skipping " << osi->distance() << endl; 
      continue; 
    }

    // Check if we are mapping to an equivalent orthodistance or a conjugate one. 
    ori = (fabs(mod_pi(od.imag - osi->distance().imag + 1.) - 1.) < 1e-4) ? 0 : 1; 

    for (i=0; i<2; i++) {
      T = get_T_matrix(*osi, i, ivls);

      T_points.push_back(T * O31_origin);

      for (gdir = 0; gdir < 2; gdir++) {
	if (!symmetric_od) {
	  S = T * T0[gdir][ori];
	  sym = group ? is_symmetry(domain, group, S) : 
	    is_symmetry_of_geodesic_lifts(domain, ivls, S, report); 
	  if (sym) {
	    matrices.push_back(S); 
	    if (group && report)
	      cout << oln << ':' << i << ':' << gdir << ':' << ori << endl; 
	  }
	} else {
	  for (ori=0; ori < 2; ori++) {
	    S = T * T0[gdir][ori];
	    sym = group ? is_symmetry(domain, group, S) : 
	      is_symmetry_of_geodesic_lifts(domain, ivls, S, report); 
	    if (sym) {
	      matrices.push_back(S);
	      if (group && report) 
		cout << oln << ':' << i << ':' << gdir << ':' << ori << endl; 
	    }
	  }
	}
      }
    }
  }

  // Check for coincident T-points: at one time I thought these might
  // result in symmetries being missed.  Actually I don't think so now
  // but I still want to check for them just in case.

  for (i=0; i<T_points.size(); i++) {
    for (j=i+1; j<T_points.size(); j++) {
      if (-(T_points[i] * T_points[j]) < 1.001) {
	return 5;
      }
    }
  }
  if (report) cout << "Symmetry finding successful.\n"; 
  return 0; 
}

static inline double safe_acosh(double x)
{
  return (x < 1.+1e-8) ? 0. : acosh(x); 
}

static Tile** make_tile_array(vector<O31_matrix> elts, Tile** array)
{
  int i, n=elts.size();
  Tile* a     = new Tile  [n];
  Tile** ptrs = new Tile* [n];
  
  *array = a; 

  for (i=0; i<n; i++) {
    ptrs[i] = &a[i]; 
    a[i].g = elts[i]; 
    a[i].length = complex_length(elts[i]);
    a[i].parity = (determinant(elts[i]) > 0.) ?
      orientation_preserving : orientation_reversing; 
  }
  return ptrs; 
}

static void free_tile_array(Tile* array, Tile** ptrs)
{
  delete [] array;
  delete [] ptrs; 
}

/* 
   Find the "axes" of a symmetry. The set of lifts of a symmetry S
   is the set of elements GS for G in the fundamental group (or 
   equivalently SG since S normalizes G). If S and all G are orientation
   preserving then (for S non-trivial) each GS has an axis. 
   The axes of the GS project to geodesics of the quotient manifold. 
   We call these the axes of S in the manifold. If S has fixed points
   then it will have a lift GS which is elliptic. If S is free we can
   still try to find the axes along which S has minimum translation 
   distance. These are the elements which this function returns. 

   Let D denote the Dirichlet domain. To find all elements GS with 
   translation distance < t it will suffice to find all G such that
   GS(D) and D can be joined by a geodesic of length < t. It follows
   that the distance between 0 and GS(0) is at most

       A = 2 * acosh( cosh(spine_radius) * cosh(t/2) ),

   (see discussion of spine radius in winged_edge.h and Dirichlet.cc).
   If we first normalize S so that S(0) is in D, it follows that G(0)
   has distance at most A + d(0,S(0)) from 0. We find all such
   elements G by tiling to the appropriate radius.
 */

void get_symmetry_axes(const WEPolyhedron* domain, O31_matrix const& Sy, 
		       dd_tiling& tiling, vector<O31_matrix>& axes)
{
  double eps = 1e-4;

  O31_matrix S  = Sy;
  O31_vector S0 = Sy * O31_origin; 
  FGWord w; // dummy. 

  // move S0 into domain, and adjust S by same element.
  dirichlet_normalize(domain, S0, S, w, 0);

  bool S_or_rev = (determinant(S) < 0.), or_rev;
  Complex length; 

  double S0_dist = safe_acosh(S0[0]);
  double t = 0.5; // find all elements with translation dist < t.
  double t_dist, min_t_dist = complex_length(S).real;
  double rqd_radius, cosh_rad;
  TilingIterator i;
  O31_matrix GS;

  axes.resize(0); 

  while (true) {

    rqd_radius = 2 * acosh(cosh(domain->spine_radius) * cosh(t/2)) + S0_dist; 

    tiling.tile(domain, rqd_radius);

    cosh_rad = cosh(rqd_radius);
    for (i=tiling.root; i; i++) {
      if (i->g(0,0) > cosh_rad) continue; // save a few cycles
      GS = i->g * S;
      length = complex_length(GS); 

      t_dist = length.real; 
      if (t_dist < min_t_dist - eps) { // found a shorter one
	axes.resize(0);
	min_t_dist = t_dist; 
      }

      // only want axes within spine_radius of the origin.
      or_rev = (i->parity==orientation_preserving)
	? S_or_rev : (!S_or_rev);
      if (distance_to_origin(GS, length, or_rev) > domain->spine_radius)
	continue; 

      if (t_dist < min_t_dist + eps) { // found same length one
	axes.push_back(GS);
      }
    }

    if (min_t_dist < t) break; 
    t += 0.5;
  }

  // eliminate conjugates using a rather awkward function. 
  Tile *array, **ptrs; 
  int j, n = axes.size(), n_tiles = tiling.num_tiles(); 
  ptrs = make_tile_array(axes, &array);
  eliminate_conjugates(ptrs, &n, tiling.root, n_tiles, domain->spine_radius);
  axes.resize(n);
  for (j=0; j<n; j++) {
    axes[i] = ptrs[i]->g;
  }
  free_tile_array(array,ptrs);
}

void print_symmetry_info(const WEPolyhedron* domain, 
			 dd_tiling& tiling, 
			 vector<GeodesicWord> const& geodesics, 
			 vector<O31_matrix> const& symms)
{
  const double eps = 1e-5; 

  int i, j, k, ord;
  O31_matrix M0, M; 
  MoebiusTransformation mt; 

  vector<O31_matrix> axes; 
  bool is_elliptic; 

  for (i=0; i<symms.size(); i++) {
    if (i>0) cout << endl; // Format nicely. 

    M0 = symms[i];
    cout << 'S' << i; 

    // compute the order of the element
    M = M0; 
    for (j=1; j<10; j++) {
      if (group_contains(domain, M, eps)) {
	cout << " has order " << j << endl; 
	break; 
      }
      M *= M0; 
    }
    if (j==10) {
      cout << " has order > 9" << endl; 
      continue; // do the next symmetry. 
    }
    ord = j; 
    
    // compute shortest translation axes
    get_symmetry_axes(domain, M0, tiling, axes);
    if (!axes.size()) continue; 

    is_elliptic = complex_length(axes[0]).real < eps;
    if (is_elliptic) {
      cout << 'S' << i << ":\n";
    } else {
      cout << 'S' << i << '^' << ord << ":\n";
    }

    for (k=0; k<axes.size(); k++) {
      M0 = axes[k];
      M = M0; 

      if (!is_elliptic) {
	// get M = M0^j which will be in the group. 
	for (j=1; j<ord; j++) M *= M0; 
      }

      mt = MoebiusTransformation(M); 
      print_element_info(domain, geodesics, mt, true);
    }
  }

}


int geodesic_symmetries(const WEPolyhedron* domain, list<interval>& ivls, 
		       vector<O31_matrix>& matrices, 
		       double max_dist, int report)
{
  dd_tiling tiling; 

  Ortholine_set os; 
  Ortholine_set::const_iterator osi;

  if (!next_ortholine_group(domain, tiling, ivls, max_dist, os, osi) ||
      osi->distance().real < 1e-4 &&
      !next_ortholine_group(domain, tiling, ivls, max_dist, os, osi)) {

    if (report) { 
      cout << "Unable to find nonzero orthodistances < " << max_dist << endl;
    }
    return 1;
  }

  symmetries_from_ort(domain, 0, ivls, os, osi, matrices, 0);

  int n; 
  for (n=0; n<matrices.size(); n++) {
    is_symmetry_of_geodesic_lifts(domain, ivls, matrices[n], 1); 
    cout << endl; 
  }

  return 0; 
}


int compute_symmetries(const WEPolyhedron* domain, const GroupPresentation* group, 
		       vector<O31_matrix>& matrices, 
		       double max_len, double max_dist, int report)
{
  dd_tiling tiling; 
  vector<GeodesicWord> geodesics;
  int g_num = -1; 
  list<interval> ivls; 
  Ortholine_set os; 
  Ortholine_set::const_iterator osi;

  if (!next_go_pair(domain, tiling, geodesics, g_num, ivls, max_len, 
		    os, osi, max_dist)) {
    if (report) { 
      cout << "Unable to find geodesics and ortholines with length in\n";
      cout << "the range 0.3-" << max_len << " and positive distance less than ";
      cout << max_dist << endl; 
    }
    return 1;
  }

  // Debugging version: I wasn't sure if ortholines with 
  // coincident T-points might not confuse the symmetry 
  // computation. Since it seems they don't, but just in 
  // case there are still problems, I'll let it recheck once
  // only in the case of coincident T-points. 

  int res, prev_num=0, tries = 2;
  
  while (tries > 0) {
    --tries;

    if (report) cout << '(' << g_num << ',' << osi->distance().real << ')' << endl;
    res = symmetries_from_ort(domain, group, ivls, os, osi, matrices, report);

    // If symmetry calculation found coincident T-points, recheck it
    // with a different set of ortholines/geodesics. See note at end of
    // compute symmetries. May remove this later.

    if (prev_num) {
      if (prev_num!=matrices.size()) {
	cerr << "Recheck failed! Symmetry group computation may be wrong!\n"; 
	cerr << "With    coincident T-points: " << prev_num << endl; 
	cerr << "Without coincident T-points: " << matrices.size() << endl; 
	break;
      } else {
	if (report) cout << "Ok.\n";
      }
    }

    if (res != 5) break;

    if (report) cout << "Ort. set had coincident T-points, rechecking..\n";
    if (!next_go_pair(domain, tiling, geodesics, g_num, ivls, max_len, os, osi, max_dist)) {
      if (report) cout << "No further ort. set available for checking.\n"; 
      break;
    }

    prev_num = matrices.size(); // save result with coincident T's. 
  }

  return 0; 
}

// We construct candidate maps from d1 to d2 and then determine if d1
// covers d2, ie. if every element of g1 maps to one in g2 (implicitly
// determined by d2).

static bool covering_from_ort(
  const WEPolyhedron* d1, const GroupPresentation* g1, 
  list<interval> const& ivls1, Ortholine const& ol1,
  const WEPolyhedron* d2, 
  list<interval> const& ivls2, Ortholine_set const& os2, Ortholine_set::const_iterator osi2, 
  O31_matrix& S)
{

  Complex od = ol1.distance(); 
  double qf_arg = quarter_fold(od.imag); 

  bool symmetric_od = fabs(od.imag) < 1e-6 || 
    fabs(fabs(od.imag) - PI_OVER_2) < 1e-6 ||
    fabs(fabs(od.imag) - PI) < 1e-6; 

  O31_matrix T; 
  O31_matrix T0[2][2]; 
  O31_matrix d1c = inverse(d1->conjugacy);

  // Get isometries from d1 coordinates to a standard position. 

  T0[0][0] = inverse(get_T_matrix(ol1, 0, ivls1)); 
  T0[0][1] = T0[0][0]; negate_row(T0[0][1], 3); // T0[0][1] reverses orientation
  T0[1][0] = T0[0][1]; negate_row(T0[1][0], 1); // T0[1][0] reverses geodesic
  T0[1][1] = T0[0][0]; negate_row(T0[1][1], 1); // T0[1][1] reverses geodesic and ori

  int ori, gdir, oln = 0; 

  vector<O31_vector> T_points; 
  
  int i; 

  // Iterate over ortholines of equal length. 
  for (; osi2 != os2.end() && (*osi2).distance().real < od.real + 1e-4; 
       osi2++, oln++) {

    // Skip any with complex orthodistance neither equal nor conjugate mod pi. 
    if (fabs(quarter_fold((*osi2).distance().imag) - qf_arg) > 1e-4) continue; 

    // Check if we are mapping to an equivalent orthodistance or a conjugate one. 
    ori = (fabs(mod_pi(od.imag - (*osi2).distance().imag + 1.) - 1.) < 1e-4) ? 0 : 1; 

    for (i=0; i<2; i++) {
      T = get_T_matrix(*osi2, i, ivls2);

      for (gdir = 0; gdir < 2; gdir++) {
	if (!symmetric_od) {
	  S = T * T0[gdir][ori];
	  if (is_subgroup(g1, d2, S*d1c, 1e-4)) return true;
	} else {
	  for (ori=0; ori < 2; ori++) {
	    S = T * T0[gdir][ori];
	    if (is_subgroup(g1, d2, S*d1c, 1e-4)) return true; 
	  }
	}
      }
    }
  }

  return false; 
}

// Strictly speaking this only checks if d1 is a covering of d2. We
// shall assume however that both manifolds have the same
// volume.  If d1 is a covering of d2 the program will not necessarily
// report this correctly either since more care is then needed in
// trying to find comparable geodesic and ortholine sets in the two
// manifolds.

bool compute_isometry(const WEPolyhedron* d1, const GroupPresentation* g1, 
		      const WEPolyhedron* d2, 
		      O31_matrix& matrix, 
		      double max_len, double max_dist, int report)
{
  dd_tiling t1, t2; 
  vector<GeodesicWord> gl1, gl2;
  int gn1 = -1, gn2 = -1; 
  list<interval> ivls1, ivls2; 
  Ortholine_set os1, os2; 
  Ortholine_set::const_iterator osi1, osi2;

  if (!next_go_pair(d1, t1, gl1, gn1, ivls1, max_len, os1, osi1, max_dist) ||
      !next_go_pair(d2, t2, gl2, gn2, ivls2, max_len, os2, osi2, max_dist)) {
    if (report) { 
      cout << "Unable to find geodesics and ortholines with length in\n";
      cout << "the range 0.3-" << max_len << " and positive distance less than ";
      cout << max_dist << endl; 
    }
    return false;
  }

  // Check we found comparable geodesic/ortholine sets in each manifold. 

  if (fabs(gl1[gn1].length.real  - gl2[gn2].length.real) > 1e-4) return false;
  if (fabs((*osi1).distance().real - (*osi2).distance().real) > 1e-4) return false; 

  return covering_from_ort(d1, g1, ivls1, *osi1, d2, ivls2, os2, osi2, matrix);
}

static void get_commens_length_geodesics(vector<GeodesicWord> const& geodesics, int gn, vector<int>& results)
{
  double length = geodesics[gn].length.real, ratio;
  const double eps = 1e-9;
  const double conf = 1e-5;
  int i, n = geodesics.size();
  long num, den;
  for (i=0; i<n; i++) {
    ratio = geodesics[i].length.real/length;
    if (appears_rational(ratio-eps, ratio+eps, conf, &num, &den))
      results.push_back(i); 
  }
}

void print_symmetry_action(const WEPolyhedron* domain, vector<O31_matrix> const& mats, vector<GeodesicWord> const& geodesics, int gn, bool all_commens)
{
  int n; 
  list<interval> ivls; 

  if (all_commens) {
    vector<int> commens;
    get_commens_length_geodesics(geodesics, gn, commens);
    int N=commens.size();
    for (n=0; n<N; n++) {
      get_crossing_lifts(domain, GSpec(commens[n], &geodesics[commens[n]]), ivls,true);
    }
  } else {
    int fg, lg; 
    same_length_range(geodesics, gn, fg, lg); 

    // Get the set of intervals corresponding to these geodesics. 
    for (n=fg; n<lg; n++) {
      get_crossing_lifts(domain, GSpec(n, &geodesics[n]), ivls, true);
    }
  }

  for (n=0; n<mats.size(); n++) {
    is_symmetry_of_geodesic_lifts(domain, ivls, mats[n], 1); 
    cout << endl; 
  }
  return; 
}

#if 0
void print_symmetry_action(const WEPolyhedron* domain, const GroupPresentation* group, 
			   vector<GeodesicWord> const& geodesics, int gn)
{
  vector<O31_matrix> mats;
  if (compute_symmetries(domain, group, mats)!=0) {
    cout << "Problem computing symmetries.\n";
    return;
  }

  int fg, lg; 
  same_length_range(geodesics, gn, fg, lg); 

  int g_num, im_num, mult; 
  line l; 
  O31_line L, lim; 
  FGWord wd;
  GeodesicWord gw; 

  // Print a line of symmetry numbers. 
  int j;
  cout << '\t'; 
  for (j=0; j<mats.size(); j++)
    cout << 'S' << j << '\t'; 
  cout << endl;

  // For each geodesic print image geodesic numbers.

  for (g_num=fg; g_num < lg; g_num++) {
    cout << g_num << '\t';
    gw = geodesics[g_num];
    fixed_points(gw.mt.matrix, l); 
    L = O31_line(l); 
    
    for (j=0; j<mats.size(); j++) {
      
      lim = mats[j] * L; 
      if (find_geodesic(domain,geodesics,lim,fg,im_num,mult,wd)) {
	if (mult < 0) cout << '-';
	cout << im_num << '\t';
      } else {
	cout << "?\t";
      }
    }
    cout << endl;
  }
}
#endif

static int num_odist_to_distinguish(geo_orbit const& a, geo_orbit const& b)
{
  if (!complex_small(a.length-b.length)) return 0;
  int i; 
  for (i=0; i < 5; i++) {
    if (!complex_small(a.odist(i) - b.odist(i))) return i+1; 
  }
  cout << "Warning: indistinguishable geo_orbits found.\n";
  return i; 
}

static void print_orbits(ostream& out, vector<geo_orbit>::const_iterator i, 
			 vector<geo_orbit>::const_iterator e, int ind) 
{
  int nod_for_prev, nod_for_next = 1, nod; 
  
  if (i==e) return; 
  for (; i != e-1; i++) {
    nod_for_prev = nod_for_next;
    nod_for_next = num_odist_to_distinguish(*i, *(i+1));
    if (nod_for_next < 1) nod_for_next = 1; // Always show at least one. 
    nod = max(nod_for_prev, nod_for_next); 
    out << '<' << ind << '>';
    if (ind < 10) out << ' '; // Line things up nicely. 
    i->print(out, nod); 
    ind++; 
  }

  out << '<' << ind << '>';
  if (ind < 10) out << ' '; // Line things up nicely. 
  i->print(out, nod_for_next); 
}

// info = 0, don't bother sorting orbits within same-length groups
// info = 1, compute only enough ortholengths to sort the orbits. 
// info = 2, compute at least 3 ortholengths per orbit. 

bool symmetry_orbits(const WEPolyhedron* domain, 
		     dd_tiling& tiling, vector<GeodesicWord>& geodesics, 
		     vector<O31_matrix> const& mats,
		     vector<geo_orbit>& orbits, 
		     int n_orbits, int n_geod, double max_length, 
		     int info, int report)
{
  // Check if manifold is amphicheiral. (Should be a better way than this!)
  int i, n_symms = mats.size();
  bool amphicheiral = false; 
  MoebiusTransformation mt; 
  for (i=0; i<n_symms; i++) {
    mt = MoebiusTransformation(mats[i]);
    if (mt.parity == orientation_reversing) {
      amphicheiral = true; 
      break; 
    }
  }

  int gn=0, next_gn, first_gn, j, n, im_num, mult, prev_num_orbits;
  FGWord wd; 
  GeodesicWord gw;
  line l; 
  O31_line L, lim; 
  bool first_orbit, only_orbit, conjugate;  
  Complex dist; 
  vector<Ortholine> ort; 
  if (n_orbits <= 0) n_orbits = 10000; 
  if (n_geod <= 0) n_geod = 10000; // ie. infinity.

  // If the cutoff is by max_length, precompute all required geodesics. 
  // unless > 4, in which case they will be computed incrementally below.
  if (n_orbits==10000 && n_geod==10000) {
    if (max_length <= 4.) {
      if (!compute_geodesics(domain, tiling, max_length, geodesics)) {
	cout << "Problem computing geodesics in symmetry_orbits.\n";
	return false; 
      }
    }
  }

  while (orbits.size() < n_orbits && gn < n_geod) {

    if (gn==geodesics.size()) { // Get some more geodesics to work with. 
      if (!compute_geodesics(domain, tiling, gn+10, geodesics)) {
	cout << "Problem computing geodesics in symmetry_orbits.\n";
	return false; 
      }
    }

    if (geodesics[gn].length.real > max_length) break; 

    same_length_range(geodesics, gn, first_gn, next_gn);
    n = next_gn - first_gn;

    vector<bool> geodesic_done(n, false);

    first_orbit = true; 
    prev_num_orbits = orbits.size();

    while (true) { // break when geodesics in [first_gn, next_gn-1] are done.

      // Look for a geodesic not yet done. 
      for (i=0; i<n; i++) 
	if (!geodesic_done[i]) break;
      if (i==n) break; 
      gn = first_gn+i; 

      // Compute its end points. 
      gw = geodesics[gn];
      fixed_points(gw.mt.matrix, l); 
      L = O31_line(l); 

      // Start a new orbit. 
      orbits.push_back(geo_orbit(domain, &geodesics, gn));
      conjugate = (amphicheiral && gw.length.imag < -1e-10); 
      orbits.back().conjugate = conjugate; 
      orbits.back().length = conjugate ? complex_conjugate(gw.length) : gw.length;

      // Compute geodesics in the orbit of gn. 
      for (j=0; j<n_symms; j++) {
	lim = mats[j] * L; 
	if (!find_geodesic(domain,geodesics,lim,first_gn,im_num,mult,wd)) {
	  cout << "Unable to identify symmetric image " << j << " of geodesic " << gn << endl;
	  continue; 
	}
	if (im_num < first_gn || im_num >= next_gn) {
	  return false; 
	}
	if (!geodesic_done[im_num - first_gn]) 
	  orbits.back().gnum.push_back(im_num); 
	geodesic_done[im_num - first_gn] = true; 
      }

      // Check if this is the only orbit in same length range
      if (first_orbit) {
	only_orbit = (orbits.back().gnum.size() == n); 
	first_orbit = false; 
      }

#if 0
      // Get orthodistances if required. 
      if (info==2 || (info==1 && !only_orbit)) {
	ort = geodesic_ortholines(domain, tiling, geodesics[gn], gn, 3);
	for (j=0; j<ort.size(); j++) {
	  dist = ort[j].distance();
	  orbits.back().odist.push_back(conjugate ? complex_conjugate(dist) : dist);
	}
      }
#endif
    }

    // Now sort the orbits we've just added. 
    if (!only_orbit && info>0)
      sort(orbits.begin()+prev_num_orbits, orbits.end());

    if (report) {
      print_orbits(cout, orbits.begin()+prev_num_orbits, orbits.end(), prev_num_orbits);
    }

    gn = next_gn;
  }
  return true;
}

#if 0
static void sort_same_length_group(const WEPolyhedron* domain, 
				   dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
				   vector<geo_orbit>& orbits, int first_on, int next_on)
{
  Complex dist; 
  bool conjugate; 
  vector<Ortholine> ort; 
  int i, j, gn; 
  for (i=first_on; i<next_on; i++) {
    gn = orbits[i].gnum[0];
    ort = geodesic_ortholines(domain, tiling, geodesics[gn], gn, 3);
    conjugate = !(orbits[i].length == geodesics[gn].length); 
    orbits[i].odist.clear(); 
    for (j=0; j<ort.size(); j++) {
      dist = ort[j].distance();
      orbits[i].odist.push_back(conjugate ? complex_conjugate(dist) : dist);
    }
  }
  sort(orbits.begin()+first_on, orbits.begin()+next_on); 
}
#endif


static bool orbit_same_length_range(vector<GeodesicWord> const& geodesics, 
			     vector<geo_orbit> const& orbits, 
			     int gn, int& first_on, int& next_on)
{
  int first_gn, next_gn; 
  same_length_range(geodesics, gn, first_gn, next_gn); 

  int i; 
  for (i=0; i<orbits.size(); i++) {
    if (orbits[i].gnum[0] >= first_gn) {
      first_on = i; 
      break; 
    }
  }
  if (i==orbits.size()) return false; 
  for (;i<orbits.size(); i++) {
    if (orbits[i].gnum[0] >= next_gn) break; 
  }
  next_on = i; 
  return true; 
}

int find_orbit_num(const WEPolyhedron* domain, 
		   dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
		   vector<geo_orbit>& orbits, int gn)
{
  int first_on, next_on; 
  if (!orbit_same_length_range(geodesics, orbits, gn, first_on, next_on)) return -1; 

  // Orbit is alone so no need to sort it. 
  if (next_on == first_on + 1) return first_on; 

  // If any orbit in this same length group already has
  // orthodistance information it is (hopefully) safe to 
  // assume that this information is present for all orbits
  // in the same length group and that they have already
  // been sorted. 
  if (!orbits[first_on].has_od()) {
    sort(orbits.begin()+first_on, orbits.begin()+next_on); 
  }

  int i, j;
  for (i=first_on; i<next_on; i++) {
    for (j=0; j<orbits[i].gnum.size(); j++) {
      if (orbits[i].gnum[j] == gn) return i; 
    }
  }
  return -1; 
}

geo_orbit find_orbit(const WEPolyhedron* domain, 
		     dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
		     vector<geo_orbit>& orbits, int on)
{
  // If number out of range, return empty answer. 
  if (on < 0 || on >= orbits.size()) return geo_orbit(); 

  // If any orbit in this same length group already has
  // orthodistance information it is (hopefully) safe to 
  // assume that this information is present for all orbits
  // in the same length group and that they have already
  // been sorted. 
  if (orbits[on].has_od()) return orbits[on]; 

  // Since sorting preserves the same-length groups, orbit[on] before sorting 
  // has the same length as orbit[on] after sorting. Get the limits of the 
  // required same-length group. 
  int gn = orbits[on].gnum[0]; 
  int first_on, next_on; 
  if (!orbit_same_length_range(geodesics, orbits, gn, first_on, next_on)) 
    return geo_orbit(); 

  // Orbit is alone so no need to sort it. 
  if (next_on == first_on + 1) return orbits[on]; 

  sort(orbits.begin()+first_on, orbits.begin()+next_on); 

  return orbits[on]; 
}

bool geo_orbit::has_od() const
{
  if (!OC) return false; 
  return OC->max_requested() >= 0; 
}

Complex geo_orbit::odist(int n) const
{
  if (!OC) return Zero; 
  Ortholine ol; 
  if (OC->get_orth(n, ol) < 2) {
    cout << "max tiling radius reached, out of ortholines\n";
    return Zero; 
  }

  return conjugate ? 
    complex_conjugate(ol.distance()) : ol.distance();
}

bool operator < (geo_orbit const& a, geo_orbit const& b)
{
  int sgn = uhs_compare(a.length,b.length,1e-6);
  if (sgn) return sgn < 0; 
  int i; 
  for (i=0; i < 8; i++) { // compare at most 8 orthodistances
    sgn = uhs_compare(a.odist(i), b.odist(i), 1e-6);
    if (sgn) return sgn < 0;
  }
  return false; 
}


void geo_orbit::print(ostream& out, int num_odist) const
{
  out << length;
  int i; 

  for (i=0; i<num_odist; i++) {
    out << odist(i);
  }

  out << ' '; 
  for (i=0; i<gnum.size(); i++) {
    out << gnum[i]; 
    if (i < gnum.size()-1) out << ',';
  }
  out << endl; 
}

ostream& operator << (ostream& out, geo_orbit const& a)
{ a.print(out, 1); return out; }

ostream& operator << (ostream& out, vector<geo_orbit> const& orbits)
{ print_orbits(out, orbits.begin(), orbits.end(), 0); return out; }



#if 0

static int compute_symmetries(const WEPolyhedron* domain, const GroupPresentation* group, 
		       dd_tiling& tiling, list<interval> const& ivls, double radius, 
		       vector<O31_matrix>& matrices, int report);
static int compute_symmetries(const WEPolyhedron* domain, const GroupPresentation* group, 
		       list<interval> const& ivls, 
		       Ortholine_set const& os, Ortholine_set::const_iterator osi, 
		       vector<O31_matrix>& matrices, int report);

// Our strategy for finding symmetries of a manifold is to compute a
// set of equal length geodesics and then compute a set of equal
// length ortholine segments. Any symmetry of the manifold must
// permute these geodesics and also permute the ortholine
// segements. Call the point where an ortholine segment meets a
// geodesic a T-point. Then these points are permuted. If they are all
// distinct then for each pair of T-points there are most four
// symmetries mapping one to the other. These can be characterized by
// whether they preserve or reverse the direction of the geodesic and
// whether they preserve or reverse the orientation of the
// manifold. Once we have a finite list of candidate symmetries we
// simply use is_symmetry() to filter down to a list of actual
// symmetries.

// At the top level we try each group of equal length geodesics in turn. 

// 0 success
// 1 lower level computation failed for every group of geodesics. 
// 2 get_crossing_lifts failed for some group of geodesics. 

int compute_symmetries(const WEPolyhedron* domain, const GroupPresentation* group, 
		       dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
		       vector<O31_matrix>& matrices, 
		       double ort_radius, int report) 
{
  // Repeat code in compute_ortholines. Have to do this rather than
  // calling it because we also want the list of intervals. 

  int gn = 0, i, next_gn, N = geodesics.size(); 

  list<interval> ivls;
  while (gn < N) {
    same_length_range(geodesics, gn, i, next_gn);

    // Skip short geodesics, they cause problems in ortholine
    // computation.
    if (geodesics[gn].length.real < 0.3) {
      gn = next_gn; 
      continue;
    }

    for (; i<next_gn; i++) 
      get_crossing_lifts(domain, group, geodesics, i, ivls);

    if (!ivls.size()) {
      cout << "Get_crossing_lifts failed for geodesics " << gn << " to " <<(next_gn-1)<< endl;
      return 2; 
    }

    if (compute_symmetries(domain, group, tiling, ivls, ort_radius, matrices, report)==0)
      return 0; 

    ivls.clear();
    gn = next_gn; 
  }

  if (report) cout << "Compute_symmetries failed for all geodesics.\n";

  return 1; 
}


static int compute_symmetries(const WEPolyhedron* domain, const GroupPresentation* group, 
		       dd_tiling& tiling, list<interval> const& ivls, double radius, 
		       vector<O31_matrix>& matrices, int report) 
{
  // Ensure tiling is big enough. 
  tiling.tile(domain, 2.0 * outradius(ivls) + radius); 

  // Get the ortholines. 
  Ortholine_set os; 
  get_ortholines(tiling.root, ivls, 0., tiling.radius, os, radius);

  Ortholine_set::const_iterator osi = os.begin(); 

  // Skip any ortholines with real orthodistance zero. 
  while (osi != os.end() && osi->distance().real < 1e-4) osi++; 

  double dist;

  while (1) {
    if (osi == os.end()) {
      if (report) cout << "No ortholines found under radius " << radius << endl;
      return 1; 
    }

    if (symmetries_from_ort(domain, group, ivls, os, osi, matrices, report)==0) return 0; 

    dist = osi->distance().real;

    // Skip to the end of the current group of ortholines of equal length. 
    while (osi != os.end() && osi->distance().real < dist + 1e-4) osi++;
  }
}

#endif
