#ifndef _dd_tiling_
#define _dd_tiling_
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


#include "ortholine.hh"
#include <vector>
#include <iostream>

using std::vector;
using std::ostream;

class Tile;
class odist_cache;
class dd_tiling;

class geo_orbit {
  odist_cache *OC;
public:
  bool conjugate; 
  Complex length;
  vector<int> gnum;

  geo_orbit() : OC(0) {}
  geo_orbit(const WEPolyhedron* domain, const vector<GeodesicWord> *gw, int gn);
  geo_orbit(geo_orbit const& go); 

  ~geo_orbit();

  geo_orbit& operator = (geo_orbit const& go); 

  bool has_od() const;
  Complex odist(int n) const; 
  int get_orth(int n, Ortholine& o) const; 

  void print(ostream& out, int num_odist=1) const; 
};

bool operator < (geo_orbit const& a, geo_orbit const& b);
ostream& operator << (ostream& out, geo_orbit const& o);
ostream& operator << (ostream& out, vector<geo_orbit> const& orbits);

// Dirichlet domain. 

// To compute a Dirichlet domain do:
// WEPolyhedron* domain = Dirichlet(manifold, dirichlet_epsilon, FALSE, FALSE);

// Except where specifically noted otherwise, all matrices, Moebius transformations 
// and geodesics are in Dirichlet domain coordinates. 

void print_Dirichlet(const WEPolyhedron* domain, const GroupPresentation* g, 
		     FILE* fp=stdout, unsigned int opts=0x39); 

void convert_Dirichlet(const GroupPresentation* g, FILE* in, FILE* out);

void save_domain_and_symmetries(const WEPolyhedron* domain, vector<O31_matrix> const& symms, FILE* out);

bool group_contains(const WEPolyhedron* H, O31_matrix g, double eps = 1e-3);
bool group_contains(const WEPolyhedron* H, O31_matrix g, double eps, FGWord& wd);

bool is_symmetry(const WEPolyhedron* domain, const GroupPresentation* G, 
		 O31_matrix const& mx, double eps = 1e-3, bool report = false);

// Tiling. 

class dd_tiling {
public:
  Tile* root;
  double radius;

  dd_tiling() : root(0), radius(0.0) {}
  ~dd_tiling();

  void tile(const WEPolyhedron* domain, double tile_rad, int report=0);
  void clear();
  int num_tiles() const; 
};

// Geodesics. 
bool compute_geodesics(const WEPolyhedron* domain, dd_tiling& tiling, 
		       double max_length, vector<GeodesicWord>& geodesics);
bool compute_geodesics(const WEPolyhedron* domain, dd_tiling& tiling, 
		       int min_n, vector<GeodesicWord>& geodesics);

void print_element_info(const WEPolyhedron* domain, 
			vector<GeodesicWord> const& geodesics, 
			MoebiusTransformation const& mt, bool dd_coords);

void get_crossing_lifts(const WEPolyhedron* domain,
			GSpec const& g, 
			list<interval>& ivls, bool full_set);

// Print_element_info expects mt in fundamental group coordinates rather than 
// Dirichlet domain coordinates. 

// Ortholines. 

// compute_ortholines(...), finds a set of ortholines, possibly between
//   multiple geodesics, up to a given radius. 
// geodesic_ortholines(...), finds at least n ortholines for one geodesic. 
// copy_ortholines(...), is just a utility to extract a vector of ortholines
//   from a set. 

void compute_ortholines(const WEPolyhedron* domain, 
			dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
			vector<int> const& vgn, double radius, Ortholine_set& os);

void compute_ortholines(const WEPolyhedron* domain, dd_tiling& tiling, 
			list<interval> const& ivls, double radius, 
			Ortholine_set& os, bool check_range=true);

vector<Ortholine> geodesic_ortholines(const WEPolyhedron* domain,
      const vector<GeodesicWord> *gw, int gn, int n);

vector<Ortholine> copy_ortholines(Ortholine_set const& OS, double cutoff);


// Symmetries.

void print_go_pairs(const WEPolyhedron* domain, 
		    double max_len, double max_dist); // Diagnostic. 

int compute_symmetries(const WEPolyhedron* domain, const GroupPresentation* group, 
		       vector<O31_matrix>& matrices, 
		       double max_len=3.0, double max_dist=2.0, int report=0);

void print_symmetry_action(const WEPolyhedron* domain, vector<O31_matrix> const& mats, vector<GeodesicWord> const& geodesics, int gn, bool all_commens=false);

void print_symmetry_info(const WEPolyhedron* domain, 
			 dd_tiling& tiling,
			 vector<GeodesicWord> const& geodesics, 
			 vector<O31_matrix> const& symms);

bool symmetry_orbits(const WEPolyhedron* domain, 
		     dd_tiling& tiling, vector<GeodesicWord>& geodesics, 
		     vector<O31_matrix> const& mats,
		     vector<geo_orbit>& orbits, int n_orbits, int n_geod, double max_length, 
		     int info, int report);

int find_orbit_num(const WEPolyhedron* domain,
		   dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
		   vector<geo_orbit>& orbits, int gn);

geo_orbit find_orbit(const WEPolyhedron* domain,
		     dd_tiling& tiling, vector<GeodesicWord> const& geodesics, 
		     vector<geo_orbit>& orbits, int on);

int geodesic_symmetries(const WEPolyhedron* domain, list<interval>& ivls, 
		       vector<O31_matrix>& matrices, 
			double max_dist, int report);

// Isometry checker. 

bool compute_isometry(const WEPolyhedron* d1, const GroupPresentation* g1, 
		      const WEPolyhedron* d2,
		      O31_matrix& matrix, 
		      double max_len=3.0, double max_dist=2.0, int report=0);



//    void gv_print_Dirichlet(FILE* fp, vector<int> const& geodesics) const;
//    void gv_print_Dirichlet(FILE* fp, int gnum, int n_radii) const;

//  void it_with_dd_tiling::gv_print_Dirichlet(FILE* fp, vector<int> const& show) const
//  {
//    gv_print_w_geodesics(fp, domain, u_group(), _geodesics, show); 
//  }

//  void it_with_dd_tiling::gv_print_Dirichlet(FILE* fp, int gn, int n_radii) const
//  {
//    if (gn < 0 || gn > num_geodesics()) {
//      cout << "Invalid geodesic number: " << gn << endl; return; }

//    ((it_with_dd_tiling*)this)->cache_ortholines(gn, tiling_radius);
//    gv_print_geodesic_w_ortholines(fp, domain, u_group(), _geodesics[gn], os, n_radii);
//  }


#endif
