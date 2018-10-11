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
#include <cstdio>
#include <iomanip>
#include "env_D.hh"
#include "snappea/length_spectrum.h"
#include "snappea/kernel.h"
#include "kernel_extras.hh"
#include "snappea/fundamental_group.h"
#include "closed_geodesics.hh"
#include "gv_print_Dirichlet.hh"
#include "snappea/Dirichlet2D.h"
#include "find_commens.hh"
#include "commens.hh"
#include "helpers.hh"

using std::sscanf;
using std::cout;
using std::cerr;
using std::endl;

#define DUPLICATE_LENGTH_EPSILON    1e-5

double env_D::dirichlet_epsilon = 1e-7;
unsigned int env_D::Dirichlet_print_options = 0x30;

enum {
  SetDirichletEpsilon =250,
  SetDirichletPrint,

  ConvertDirichlet, // D-D not needed

  ComputeDirichlet, // Dirichlet domain needed.
  Perturb,
  PrintDirichlet,
  SaveDirichlet,
  TileOut, 

  PrintTileRadius,
  PrintNumTiles,

  PrintLengthSpectrum,
  PrintGeodesics,
  PrintInjectivityRadii,
  PrintOrtholines,
  PrintOrthoTraces,
  PrintK,
  PrintGOPairs,
  PrintIntervals,
  
  NewOrtholines, 

  PrintSymmetries,
  PrintSymmetryMatrices,
  PrintSymmetryInfo,
  PrintSymmetryOrbits,
  PrintSymmetryAction, 
  PrintSymmetryGroupAction,
  SaveSymmetryMatrices,
  PrintSymmetryGroupTable,
  OrbitNumber,
  PrintOrbit,

  GeodesicSymmetries,

  Evaluate,
  EvaluateOrthodistance,
  FindIsometry, 

  GvSaveGeodesics,
  GvSaveOrtholines,

  TestNonArithmetic,

  FindCore, // Triangulation needed
  HiddenSymmetryAction,

  FlatDirichlet,
  SaveFlatDirichlet
};

void env_D::setup_menu()
{
  env_T::setup_menu();

  mm.add_item("set dirichlet epsilon", SetDirichletEpsilon); 
  mm.add_item("set dirichlet print_options", SetDirichletPrint); 
  mm.add_item("convert dirichlet", ConvertDirichlet);
  mm.add_item("flat dirichlet", FlatDirichlet); 
  mm.add_item("save flat dirichlet", SaveFlatDirichlet); 
  mm.add_item("compute dirichlet", ComputeDirichlet); 
  // mm.add_item("perturb_dirichlet", Perturb);
  mm.add_item("print dirichlet", PrintDirichlet); 
  mm.add_item("save dirichlet", SaveDirichlet); 
  mm.add_item("tile", TileOut); 
  mm.add_item("print tile radius", PrintTileRadius); 
  mm.add_item("print number tiles", PrintNumTiles);
  mm.add_item("print length_spectrum", PrintLengthSpectrum);
  mm.add_item("print geodesics", PrintGeodesics);
  mm.add_item("print injectivity_radii", PrintInjectivityRadii); 
  mm.add_item("print ortholines", PrintOrtholines);
  mm.add_item("print traces", PrintOrthoTraces); 
  mm.add_item("print k", PrintK);
  mm.add_item("print go_pairs", PrintGOPairs); 
  mm.add_item("print intervals", PrintIntervals);  

  mm.add_item("ortholines", NewOrtholines);

  mm.add_item("print number symmetries", PrintSymmetries); 
  mm.add_item("print symmetry matrices", PrintSymmetryMatrices); 
  mm.add_item("print symmetry info", PrintSymmetryInfo); 
  mm.add_item("print symmetry orbits", PrintSymmetryOrbits); 
  mm.add_item("print symmetry action", PrintSymmetryAction); 
  mm.add_item("symmetry", PrintSymmetryGroupAction);
  mm.add_item("save symmetry_group", SaveSymmetryMatrices); 
  mm.add_item("print symmetry group", PrintSymmetryGroupTable); 
  mm.add_item("orbit number", OrbitNumber); 
  mm.add_item("print orbit", PrintOrbit); 

  mm.add_item("geodesic symmetries", GeodesicSymmetries); 
  mm.add_item("print -h symmetry action", HiddenSymmetryAction); 

  mm.add_item("evaluate word", Evaluate); 
  mm.add_item("evaluate orthodistance", EvaluateOrthodistance); 
  mm.add_item("isometry", FindIsometry);
  mm.add_item("find core", FindCore); 

  mm.add_item("gv_save geodesics", GvSaveGeodesics); 
  mm.add_item("gv_save ortholines", GvSaveOrtholines); 

  mm.add_item("test -l arithmetic", TestNonArithmetic); 
}

void env_D::validate_event(int& what)
{
  env_T::validate_event(what); 

  if (what > SaveFlatDirichlet) return; 

  if (what >= FindCore && !T) { 
    cout << "This function requires a triangulation.\n";
    what = Nothing; 
    return; 
  }

  if (what >= PrintDirichlet) {
    if (!get_group()) {
      cout << "This function requires a manifold.\n";
      what = Nothing;
    } else if (!get_Dirichlet_domain()) {
      what = Nothing; 
    }
  }
}

static void chop(Complex& z, double eps)
{
  if (fabs(z.real) < eps) z.real = 0.0; 
  if (fabs(z.imag) < eps) z.imag = 0.0; 
}

// User may input a word, either as a word in the generators of the 
// fundamental group, or as a geodesic number. 

bool get_word(FGWord& w, vector<GeodesicWord> const& geodesics)
{
  string s; 
  if (!get_input(s, "word in fundamental group ? "))
    return false; 
  int n; 
  if (isdigit(s[0])) {
    if (sscanf(s.c_str(), "%d", &n)!=1) return false; 
    if (n < 0 || n >= geodesics.size()) {
      cout << "Geodesic number " << n << " out of range.\n";
      return false; 
    } 
    w = geodesics[n].word;
    return true; 
  }
  w = FGWord(s.c_str()); 
  return true; 
}

#if 0
static bool select_a_geodesic(vector<GeodesicWord> const& geodesics, int& gnum, int& n_radii)
{
  string s; 
  if (geodesics.size()==0) {
    printf("Use \"print geodesics\" to compute geodesics first\n"); 
    return false;
  }
  return get_input(s,"which geodesic ? ") &&
    (sscanf(s.c_str(), "%d", &gnum) == 1) && 
      get_input(s,"how many radii ? ") &&
	(sscanf(s.c_str(), "%d", &n_radii) == 1); 
}
#endif

static bool select_a_radius(vector<GeodesicWord> const& geodesics, int& gnum, double& radius)
{
  string s; 
  if (geodesics.size()==0) {
    printf("Use \"print geodesics <radius>\" to compute geodesics first\n"); 
    return false;
  }
  return get_input(s,"which geodesic ? ") &&
    (sscanf(s.c_str(), "%d", &gnum) == 1) && 
      get_input(s,"cutoff radius ? ") &&
	(sscanf(s.c_str(), "%lf", &radius) == 1); 
}

// The idea is to modify G by adding another generator S
// which is supposed to be a symmetry (in the normalizer)
// of G. 

static bool symmetry_quotient(GroupPresentation* G, const WEPolyhedron* H, O31_matrix const& S)
{
  double eps = 1e-4; 
  O31_matrix c(S*inverse(H->conjugacy));
  O31_matrix c_inv(inverse(c)); 
  FGWord w;
  int i, n = fg_get_num_generators(G); 
  for (i=0; i<n; i++) {
    if (!group_contains(H, c * fg_gen(G,i) * c_inv, eps, w)) 
      return false; 
  }

  // The new generator represents S.
  int ng = add_generator(G, H->conjugacy*c);

  FGWord nrel; 
  for (i=0; i<n; i++) {
    group_contains(H, c * fg_gen(G,i) * c_inv, eps, w);

    // The new relations are of the form: S * gen[i] * S^-1 * w.
    nrel = FGWord();
    nrel *= ng;
    nrel *= i+1;
    nrel *= -ng;
    nrel *= w; 

    add_relation(G, to_cyclic_word(nrel));
  }
  return true; 
}


env_D::~env_D() 
{
  if (domain) my_free(domain); 
}

bool env_D::get_Dirichlet_domain()
{
  if (domain) return true; 
  if (!get_group()) return false; 

  int n; 
  O31_matrix* o31_generators; 
  MoebiusTransformation* Moebius_generators; 
  FGWord* words=0; 
  if (T && simplify) { 
    int i; 
    n = get_num_generators(T->M()); 
    Moebius_generators = NEW_ARRAY(n, MoebiusTransformation);
    o31_generators = NEW_ARRAY(n, O31_matrix);
    words = NEW_ARRAY(n, FGWord); 
    matrix_generators(T->M(), Moebius_generators, FALSE, TRUE);
    for (i=0; i<n; ++i) {
      o31_generators[i] = Moebius_generators[i]; 
      words[i] = fg_original_generator(G, i); 
    }
  } else { 
    o31_generators = fg_get_generators(G); 
    n = fg_get_num_generators(G);
  }

  domain = Dirichlet_from_generators(o31_generators, n, 
				     dirichlet_epsilon, TRUE, FALSE, words);
  if (!domain) // try again without moving basepoint.
    domain = Dirichlet_from_generators(o31_generators, n, 
				       dirichlet_epsilon, FALSE, FALSE, words);

  if (words) {
    my_free_array(words); 
    my_free_array(o31_generators);
    my_free_array(Moebius_generators); 
  }

  if (!domain) {
    cout << "Problem computing a Dirichlet domain for this group.\n";
    cout << "Try changing Dirichlet epsilon.\n"; 
    return false;
  }
  return true;
}

void env_D::clear()
{
  if (symmetries.size()) { symmetries.resize(0); }
  if (num_geodesics()) { geodesics.resize(0); }
  tiling.clear();
  if (domain) { my_free(domain); domain = 0; }
}

void env_D::conjugate(O31_matrix const& c)
{
  if (!domain) return; 
  update_conjugacy(domain, c); 
#if 0
  O31_matrix ci = inverse(c); 
  int i, n=symmetries.size(); 
  for (i=0; i<n; ++i)
    symmetries[i] = c * symmetries[i] * ci; 
#endif
}

void env_D::T_changed()
{
  clear(); 
  env_T::T_changed(); 
}

void env_D::T_shape_changed()
{
  clear();
  env_T::T_shape_changed(); 
}

bool env_D::compute_geodesics(double max_length)
{
  if (geodesics.size() > 0 && geodesics.back().length.real >= max_length) return true; 
  return ::compute_geodesics(domain, tiling, max_length, geodesics);
}

bool env_D::compute_geodesics(int min_n)
{
  return ::compute_geodesics(domain, tiling, min_n, geodesics);
}

void env_D::print_lengths(double max_len) const
{
  int i=0, n = num_geodesics(); 
  if (n==0) return;
  Complex len;
  double eps = 1e-6; 
  int j = 0; 
  if (max_len==0.) max_len = 1000.; // no cutoff specified.

  while (i < n) {
    len = geodesic(i).length; 
    if (len.real > max_len) break; 
    cout << '(' << j << ')';
    if (j < 10) cout << ' ';
    cout << len << ' ' << i; 
    i++; j++; 
    while (i < n) {
      if (fabs(geodesic(i).length.real - len.real) > eps) break;
      if (fabs(geodesic(i).length.imag - len.imag) < eps) {
	cout << ',' << i; 
      } else if (fabs(geodesic(i).length.imag + len.imag) < eps) {
	cout << ',' << i << '*'; 
      } else break; 
      i++; 
    }
    cout << endl; 
  }
}

bool env_D::find_core_geodesic(int index, Complex& clen, int& gn, int& ori, FGWord& conj) const
{
  int nc = get_num_cusps(T->M()), sing_index;
  if (index < 0 || index >= nc) {
    cerr << "Cusp/core index must be in range 0-" << (nc-1) << '.' << endl; 
    return false; 
  }

  // Get complex length of this core geodesic. 
  core_geodesic(T->M(), index, &sing_index, &clen, NULL);
  if (clen==Zero) {
    cerr << "Cusp must be filled (with integer coefficients)." << endl; 
    return false; 
  }

  // Find a lift of the same core geodesic.
  line core = core_geodesic(T->M(), G, index); 
  MoebiusTransformation C(O31_matrix(domain->conjugacy));
  core = inverse(C) * core; 
      
  // Locate it in list of geodesics. 
  return find_geodesic(domain, geodesics, core, 0, gn, ori, conj);
}

vector<Ortholine> env_D::ortholines(vector<int> const& vgn, double radius)
{
  Ortholine_set O; 
  ::compute_ortholines(domain, tiling, geodesics, vgn, radius, O);
  return copy_ortholines(O, radius); 
}

int env_D::compute_symmetries(double lmx, double dmax, int report)
{
  if (symmetries.size()) return 0; 
  return ::compute_symmetries(domain, G, symmetries, lmx, dmax, report);
}

void env_D::print_symmetry_action(int gn, const vector<O31_matrix>* mats, bool all_commens) 
{
  if (compute_symmetries()!=0) {
    cout << "Problem computing symmetries.\n"; return; }

  if (!mats) mats = &symmetries;
  ::print_symmetry_action(domain, *mats, geodesics, gn, all_commens); 
}

void env_D::print_geodesic_symmetries(vector<int> const& vgn) const
{
  // Get the intervals corresponding to this set of geodesics.
  int i, gn; 
  list<interval> ivls;
  for (i=0; i<vgn.size(); i++) {
    gn = vgn[i];
    if (gn < 0 || gn >= geodesics.size()) {
      cout << "Geodesic number " << vgn[i] << " out of range.\n";
      if (gn >= 0) cout << "Please compute more geodesics and try again.\n"; 
      return;
    }
    get_crossing_lifts(domain, GSpec(gn, &geodesics[gn]), ivls, true);
  }

  vector<O31_matrix> matrices; 
  geodesic_symmetries(domain, ivls, matrices, 4.0, 1); 
}

// info = 0, don't bother sorting orbits within same-length groups
// info = 1, compute only enough ortholengths to sort the orbits. 
// info = 2, compute at least 3 ortholengths per orbit. 

bool env_D::symmetry_orbits(vector<geo_orbit>& orbits, int n_orbits, int n_geod, double max_length, int info, int report)
{
  // Find the symmetries. 
  if (compute_symmetries()!=0) {
    cout << "Failed to find symmetries in symmetry_orbits.\n";
    return false; 
  }

  return ::symmetry_orbits(domain, tiling, geodesics, symmetries, orbits, n_orbits, n_geod, max_length, info, report); 
}

void env_D::gv_print_Dirichlet(FILE* fp, vector<int> const& show) const
{
  gv_print_w_geodesics(fp, domain, geodesics, show); 
}

void env_D::gv_print_Dirichlet(FILE* fp, int gn, int n_ort)
{
  if (gn < 0 || !compute_geodesics(gn+1)) return; 
  vector<Ortholine> ort = geodesic_ortholines(domain, &geodesics, gn, n_ort); 
  gv_print_geodesic_w_ortholines(fp, domain, geodesics[gn], ort); 
}

// This doesn't quite belong here but there is no other sensible 
// place to put it. 

Complex nearby_complex(Complex z, Complex const& z0)
{
  while (z.imag - z0.imag > PI) z.imag -= TWO_PI;
  while (z.imag - z0.imag < -PI) z.imag += TWO_PI;
  return z;
}

bool words_to_ortholines(const GroupPresentation* group, 
			 Triangulation* manifold,
			 vector<Ortholine>& ort, bool ac)
{
  int i, n = ort.size(); 
  if (n == 0) return true; 

  // We want to apply the words to a lift of the core geodesic. 
  // The right lift to use is the axis of the longitude. 

  FGWord longitude = fg_longitude(group, 0); 
  line core; 
  fixed_points(word_to_Moebius(group, longitude).matrix, core);

  // We also want the ends of core to stay in a fixed order as
  // we vary the dehn surgery. We do this by having end 0 always 
  // coincide with the base vertex of the ideal triangulation. 

  // We're assuming a 1-cusped manifold. 
  Complex basepoint = cusp_vertex(manifold, 0);
  if (same_point(basepoint, core.end[1])) {
    basepoint = core.end[1]; // OK, semi redundant, but avoids inaccuracies. 
    core.end[1] = core.end[0];
    core.end[0] = basepoint;
  } else if (!same_point(basepoint, core.end[0])) {
    return false; 
  }

  // Ortholines are computed from MoebiusTransformations. To get an
  // ortholine we need the transformation which carries [0,infinity] to
  // another geodesic. Therefore we want a normalizing transformation 
  // which takes core to [0, infinity]. We also normalize such that
  // the first_image has position 0+0*i. 

  line first_image = word_to_Moebius(group, ort[0].word) * core;
  line first_ort = ortholine(core, first_image);
  Complex a[3], b[3]; 
  b[0] = Zero; b[1] = Infinity; b[2] = One;
  a[0] = core.end[0]; a[1] = core.end[1]; a[2] = first_ort.end[1];
  MoebiusTransformation M, N(a,b); 
  MoebiusTransformation N_inv = inverse(N); 

  Complex e0, e1; 
  for (i=0; i<n; i++) {
    if (ac) { 
      e0 = ort[i].position(0); 
      e1 = ort[i].position(1); 
    }
    M = N * word_to_Moebius(group,ort[i].word) * N_inv; 
    ort[i].set_from_Moebius(M); 
    if (ac) {
      // Avoid jumping by multiples of 2*PI*I. 
      ort[i].set_to(ort[i].distance(), 
		    nearby_complex(ort[i].position(0), e0),
		    nearby_complex(ort[i].position(1), e1)); 
    } // Initial computation should give positions close to 0, not +/- PI*I. 
  }
  return true; 
}

void env_D::process_event(int what)
{
  string s;
  int j;

  switch (what) {

  case SetDirichletEpsilon:
    {
      double new_epsilon;
      if (!get_input(s,"Dirichlet epsilon ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &new_epsilon) != 1) break; 
      dirichlet_epsilon = new_epsilon; 
    }
    break; 
  case SetDirichletPrint:
    {
      int new_opts;
      if (!get_input(s,"Dirichlet print options ? ")) break; 
      if (sscanf(s.c_str(), "%i", &new_opts) != 1) break; 

      Dirichlet_print_options = new_opts; 
    }
    break;
  case ComputeDirichlet:
    verbose_Dirichlet = TRUE; 
    get_Dirichlet_domain();
    verbose_Dirichlet = FALSE; 
    break; // All done above. 

  case PrintDirichlet:
    printf("Volume: %.16g\n", domain->approximate_volume); 
    printf("Conjugacy\n"); 
    print_Dirichlet(domain, G, stdout, Dirichlet_print_options); 
    break;

  case SaveDirichlet:
    {
      FILE* fp = fopen("dirichlet.words", "w");
      if (!fp) {
	printf("Unable to open file dirichlet.words for writing.\n");
	break; 
      }
      printf("Writing file dirichlet.words.\n");
      print_Dirichlet(domain, G, fp, 0x20); 
      fprintf(fp, ".\n");
      fclose(fp); 
    }
    break;

  case ConvertDirichlet:
    {
      if (!get_group()) break; 
      FILE* in = fopen("dirichlet.words", "r");
      if (!in) {
	printf("Unable to open file dirichlet.words for reading.\n");
	break; 
      }
      FILE* out = fopen("dirichlet.mats", "w");
      if (!out) {
	printf("Unable to open file dirichlet.mats for writing.\n");
	fclose(in); 
	break;
      }
      printf("Writing file dirichlet.mats.\n"); 
      convert_Dirichlet(G, in, out); 
      fclose(out);
      fclose(in); 
    }
    break;
  case FlatDirichlet:
    {
      WE2DPolygon* polygon = Dirichlet2D(T->M(), dirichlet_epsilon, FALSE, FALSE); 
      if (!polygon) {
	printf("Problem computing 2D Dirichlet domain.\n"); 
	break;
      }

      print_polygon(polygon); 

      if (Dirichlet2D_extras(polygon)!=func_OK) {
	printf("Problem computing extra information for 2D Dirichlet domain.\n"); 
	break;
      }
      print_Dirichlet2D_extras(polygon); 

      printf("Group representation:\n");
      print_group(polygon);
      printf("Edge identifications:\n"); 
      print_label(polygon); 

      free_Dirichlet2D_domain(polygon); 
    }
    break;
  case SaveFlatDirichlet:
    {
      WE2DPolygon* polygon = Dirichlet2D(T->M(), dirichlet_epsilon, FALSE, FALSE); 
      if (!polygon || Dirichlet2D_extras(polygon)!=func_OK) {
	printf("Problem computing 2D Dirichlet domain.\n"); 
	break;
      }

      FILE* fp = fopen("dirichlet.words", "w");
      if (!fp) {
	printf("Unable to open file dirichlet.words for writing.\n");
	break; 
      }
      printf("Writing file dirichlet.words.\n");

      fprintf(fp, "# %s\n", name().c_str()); 
	
      fprintf(fp, "# "); print_polygon(polygon, fp); 
      print_Dirichlet2D_extras(polygon, fp, "# ");

      fprintf(fp, "# Edge identifications:\n# "); 
      print_label(polygon, fp); 

      fprintf(fp, "# \n[1,0,0,0;\n 0,1,0,0;\n 0,0,1,0;\n 0,0,0,1]\n\n");

      print_group(polygon, fp);
      fprintf(fp, ".\n"); 

      fclose(fp); 
      free_Dirichlet2D_domain(polygon); 
    }
    break;


  case TileOut:
    {
      double radius;
      if (!get_input(s,"radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break;
      tiling.tile(domain, radius, 1); 
    }
    break;
    
  case PrintTileRadius:
    cout << "Current tiling radius: " << current_tiling_radius() << endl;
    break;

  case PrintNumTiles:
    cout << tiling.num_tiles() << endl; 
    break;

  case PrintGeodesics:
    {
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      
      if (!compute_geodesics(radius)) break;
      int i,n = num_geodesics();
      int old_cout_prec = cout.precision(16);
      for (i=0; i<n; i++) {
	if (geodesic(i).length.real > radius) break; 
	cout << '[' << i << ']';
	if (i < 10) cout << ' ';
	cout << geodesic(i) << endl;
      }
    cout.precision(old_cout_prec);
    }
    break; 
    
  case PrintLengthSpectrum:
    {
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      
      if (!compute_geodesics(radius)) break;
      print_lengths(radius); 
    }
    break; 
    
  case Evaluate:
    {
      if (!get_group()) break; 
      FGWord wd;
      if (!get_word(wd, geodesics)) break; 
      MoebiusTransformation mt = word_to_Moebius(G, wd);
      cout << mt << endl; 
      print_element_info(domain, geodesics, mt, false);
    }
    break; 
    
  case EvaluateOrthodistance:
    {
      if (!get_group()) break; 
      FGWord wd, conj;
      if (!get_word(wd, geodesics)) break; 
      if (!get_word(conj, geodesics)) break; 
      
      MoebiusTransformation g = word_to_Moebius(G, wd);
      line l, cl;
      if (fixed_points(g.matrix, l)!=2) {
	cout << "element is not hyperbolic\n";
	break; 
      }
      MoebiusTransformation c = word_to_Moebius(G, conj);
      cl = c*l; 
      
      cout << "Axis : " << l << endl; 
      cout << "Image: " << cl << endl; 
      cout << "Orthodistance:" << complex_log(orthodist(l, cl),0.0) << endl;
    }
    break; 
    
  case FindCore:
    {
      int core_index;
      if (!get_input(s,"index of core geodesic ? ")) break; 
      if (sscanf(s.c_str(), "%d", &core_index) != 1) break; 
      
      int nc = get_num_cusps(T->M()), sing_index;
      if (core_index < 0 || core_index >= nc) {
	cout << "Cusp/core index must be in range 0-" << (nc-1) << '.' << endl; 
	break;
      }
      // Get complex length of this core geodesic. 
      Complex clen; 
      core_geodesic(T->M(), core_index, &sing_index, &clen, NULL);
      if (clen==Zero) {
	cout << "Cusp must be filled (with integer coefficients)." << endl; 
	break;
      }
      // Make sure we have enough geodesics in our list. 
      cout << "Core geodesic has length " << clen << endl;
      if (!compute_geodesics(clen.real)) {
	cout << "Sorry, problem computing geodesics." << endl; 
	break;
      }
      
      // Find a lift of the same core geodesic.
      line core = core_geodesic(T->M(), G, core_index); 
      MoebiusTransformation C(O31_matrix(domain->conjugacy));
      core = inverse(C) * core; 
      
      // Locate it in list of geodesics. 
      int gnum, ori;
      FGWord conj;
      if (!find_geodesic(domain, geodesics, core, 0, gnum, ori, conj)) {
	cout << "Sorry, unable to locate core geodesic in the current list of geodesics.\n";
	break;
      }
      cout << "Geodesic: " << gnum << '.' << endl; 
    }
    break;
    
  case PrintSymmetryOrbits:
    {
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      
      vector<geo_orbit> orbits; 
      if (!symmetry_orbits(orbits, 0, 0, radius, 2, true)) break;
    }
    break; 
    
  case PrintOrbit:
    {
      int n; 
      if (!get_input(s,"orbit number ? ")) break;
      if (sscanf(s.c_str(), "%d", &n)!=1) break;
      
      vector<geo_orbit> orbits; 
      if (!symmetry_orbits(orbits, n+1, 0, 10.0, 0, false)) {
	cout << "Problem computing symmetry orbits.\n";
	break;
      }
      cout << find_orbit(domain, tiling, geodesics, orbits, n);
    }
    break; 
    
  case OrbitNumber:
    {
      int n; 
      if (!get_input(s,"geodesic number ? ")) break;
      if (sscanf(s.c_str(), "%d", &n)!=1) break;
      
      vector<geo_orbit> orbits; 
      if (!symmetry_orbits(orbits, n+1, n+1, 10.0, 0, false)) {
	cout << "Problem computing symmetry orbits.\n";
	break;
      }
      cout << "Orbit number: " << 
	find_orbit_num(domain, tiling, geodesics, orbits, n) << endl;
    }
    break; 
    
  case PrintSymmetryAction:
    {
      int gn;
      if (!get_input(s,"geodesic number ? ")) break; 
      if (sscanf(s.c_str(), "%d", &gn)!=1) break;
      if (!compute_geodesics(gn+1)) break; 
      
      print_symmetry_action(gn);
    }
    break;
    
  case PrintInjectivityRadii:
    // Prints a list of geodesics and the tube radius of each.
    {
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 

      if (!compute_geodesics(radius)) break;
      int i,n = num_geodesics();

      vector<Ortholine> S;

      for (i=0; i<n; i++) {
	if (geodesic(i).length.real > radius) break; 

	S = geodesic_ortholines(domain, &geodesics, i, 1);
	cout << "[" << i << "]"; if (i<10) cout << ' ';
	cout << geodesic(i).length << " ";
	if (S.begin() == S.end()) {
	  cout << "> " << radius/2.0;
	} else {
	  fwprint(cout, S.front().distance().real/2.0, 0, '\0');
	}
	cout << " " << geodesic(i).word << endl; 
      }
    }
    break; 

#if 0
    // I am not sure what this function was supposed to be for. 
  case Perturb:
    {
      if (!get_Dirichlet_domain()) break; 
      double disp[3] = {.003,-.002,.0011};
      change_basepoint(&domain, 0, 0, 0, disp, dirichlet_epsilon, FALSE, TRUE, FALSE); 
      tiling.clear();
      geodesics.resize(0);
      symmetries.resize(0);
    }
    break;
#endif
  case PrintIntervals:
    { 
      int n; 
      if (!get_input(s,"geodesic number ? ")) break;
      if (sscanf(s.c_str(), "%d", &n) != 1) break; 
      if (!num_geodesics()) {
	cout << "Please compute some geodesics first\n";
	break; 
      }
      if (n >= num_geodesics()) {
	cout << "Please compute more geodesics first\n";
	break; 
      }
      if (n < 0) {
	cout << "Invalid geodesic number\n";
	break; 
      }

      list<interval> ivls; 
      list<interval>::const_iterator it; 
      O31_matrix M = domain->conjugacy; 
      get_crossing_lifts(domain, GSpec(n, &geodesic(n)), ivls, true);
      for (it=ivls.begin(); it!=ivls.end(); it++) {
	cout << line(M * it->the_line()) << ' ' << *it << endl;
      }
      break; 
    }

  case GeodesicSymmetries:
    {
      vector<int> geodesics; 
      get_input(s, "geodesics ? ", -1); // Read the rest of the line. 
      if (s.length() > 0) geodesics = read_numbers(s.c_str()); 
      if (!geodesics.size()) break; 
      print_geodesic_symmetries(geodesics); 
      break; 
    }

  case PrintOrtholines:
    {
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      
      vector<int> geodesics; 
      get_input(s, "geodesics to display ? ", -1); // Read the rest of the line. 
      if (s.length() > 0) geodesics = read_numbers(s.c_str()); 
      if (!geodesics.size()) break; 
      
      vector<Ortholine> S = ortholines(geodesics, radius);
      int old_cout_prec = cout.precision(16);
      vector<Ortholine>::const_iterator it; 
      for (it = S.begin(); it != S.end(); it++)
	cout << *it << endl;
      cout << S.size() << " ortholines\n";
      cout.precision(old_cout_prec);
    }
    break;

  case NewOrtholines:
    {
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      
      vector<int> gl; 
      get_input(s, "geodesics to display ? ", -1); // Read the rest of the line. 
      if (s.length() > 0) gl = read_numbers(s.c_str()); 
      
      list<interval> ivls; 
      int i, n = gl.size(); 
      for (i=0; i<n; ++i) 
	get_crossing_lifts(domain, GSpec(gl[i], &geodesic(gl[i])), ivls, true);
      
      Ortholine_set OS;
      new_ortholines(domain, ivls, geodesics, radius, OS);

      int old_cout_prec = cout.precision(16);
      Ortholine_set::const_iterator it; 
      for (it = OS.begin(); it != OS.end(); it++)
	cout << *it << endl;
      cout << OS.size() << " ortholines\n";

      // cout << "There were " << count_lifts(domain, ivls, radius) << " lifts\n"; 
      cout.precision(old_cout_prec);
    }
    break;
    
  case PrintOrthoTraces:
    {
      double radius;
      int gnum;
      if (!select_a_radius(geodesics, gnum, radius)) break;
      Complex trh2 = 2.0 * complex_cosh(geodesic(gnum).length) + Two; 
      chop(trh2, 1e-10);
      cout << "tr^2(geodesic(" << gnum << ")): " << trh2 << endl;

      vector<int> vgn(1, gnum); 
      vector<Ortholine> S = ortholines(vgn, radius);
      vector<Ortholine>::const_iterator it; 

      Complex trcomm;
      for (it = S.begin(); it != S.end(); it++) {
	trcomm = (trh2 - complex_cosh(it->distance()) * (trh2 - Four))/2.0;
	chop(trcomm, 1e-10); 
	cout << trcomm << endl;
      }
    }
    break;

  case PrintK:
    {
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      
      vector<int> geodesics; 
      get_input(s, "geodesics to display ? ", -1); // Read the rest of the line. 
      if (s.length() > 0) geodesics = read_numbers(s.c_str()); 
      if (!geodesics.size()) break; 
      
      vector<Ortholine> S = ortholines(geodesics, radius);
      int i, n = S.size(); 
      for (i=0; i<n; i++) {
	cout << S[i].distance() << ' '; 
	fwprint(cout, sin(S[i].distance().imag)/sinh(S[i].distance().real), 0);
	cout << ' ' << S[i].word << endl; 
      }
    }
    break;
    
  case SaveSymmetryMatrices:
    {
      if (compute_symmetries()!=0) break; 
      string fname = name(1) + ".symms";
      cout << "Writing file " << fname << endl; 
      FILE* fp = fopen(fname.c_str(), "w");
      if (!fp) {
	printf("Unable to open file for writing.\n");
	break; 
      }
      save_domain_and_symmetries(domain, symmetries, fp); 
      fclose(fp); 
    }
    break;

  case PrintSymmetryGroupTable:
    {
      SymmetryGroup *G = NEW_STRUCT(SymmetryGroup); 
      if (compute_symmetries()!=0) break; 
      int i, j, k, n = symmetries.size(); 
      G->order = n; 
      G->symmetry_list = 0;
      G->product = NEW_ARRAY(n, int *);
      for (i=0; i<n; i++)
	G->product[i] = NEW_ARRAY(n, int);

      O31_matrix mx; 
      FGWord wd; 
      bool amphicheiral = false; 
      for (i=0; i<n; i++) {
	if (!amphicheiral && determinant(symmetries[i]) < 0.) amphicheiral = true; 
	for (j=0; j<n; j++) { 
	    
	  for (k=0; k<n; k++) {
	    mx = symmetries[i]*symmetries[j]*inverse(symmetries[k]);
	    if (group_contains(domain, mx)) {
	      printf("%4d", k);
	      G->product[i][j] = k; 
	      break; 
	    }
	  }
	  if (k==n)
	    printf("   ?"); 
	}
	printf("\n"); 
      }
      fill_in_symmetry_group(G); 
      print_symmetry_group(G);
      cout << "Amphicheiral: ";
      if (amphicheiral) 
	cout << "yes" << endl; 
      else 
	cout << "no" << endl;
 
      break; 
    }

  case PrintSymmetries:
    {
      if (compute_symmetries()==0) {
	if (symmetries.size()==1) 
	  cout << "Manifold has no non-trivial symmetries.\n"; 
	else 
	  cout << "Manifold has " << symmetries.size() << " symmetries.\n"; 
      }
    }
    break;
    
  case PrintSymmetryMatrices:
    {
      int i; 
      if (compute_symmetries()==0) {
	for (i=0; i<symmetries.size(); i++) 
	  cout << MoebiusTransformation(symmetries[i]) << endl; 
      }
    }
    break;

  case PrintSymmetryGroupAction:
    {
      if (compute_symmetries()!=0) break;
      int s_num;
      if (!get_input(s, "symmetry number ? ")) break;
      if (sscanf(s.c_str(), "%d", &s_num)!=1) break;
      if (s_num < 0 || s_num >= symmetries.size()) {
	cout << "invalid symmetry number " << s_num << endl;
	break;
      }
      if (!symmetry_quotient(G, domain, symmetries[s_num])) break;
      print(G, true);

      // Print distance from identity for each relation. 
      cout << "Relation accuracies: "; 
      int nr = fg_get_num_relations(G);
      int i;
      for (i=0; i<nr; ++i) {
	cout << id_distance(word_to_Moebius(G, fg_relation(G,i))) << ' ';
      }
      cout << endl; 

      // Get a name for the quotient.
      char buf[20];
      string quot_name = name();
      sprintf(buf, ",S%d", s_num);
      quot_name += buf;

      // Discard the triangulation now, and set things
      // up as if we'd just read the group from a file. 
      GroupPresentation* tmp = G;
      G = 0;    // Avoid deleting the saved G.
      set_T(0); // Clears everything. 
      set_G(tmp);
      gp_name = quot_name;
    }
    break;

  case PrintSymmetryInfo:
    {
      ios_fmtflags old_fmt = cout.setf(ios::fixed);
      if (compute_symmetries()==0) {
	print_symmetry_info(domain, tiling, geodesics, symmetries);
      }
      cout.flags(old_fmt); 
    }
    break;

  case FindIsometry:
    {
      j = choose_manifold("Find isometry with "); 
      if (!get_T(j)->has_manifold()) {
	cout << "There is no manifold " << (j+1) << endl; 
	break; 
      }
      
      if (T && get_TT(j) && fabs(T->volume() - get_TT(j)->volume()) > 1e-5) {
	cout << "Volumes differ.\n";
	break; 
      }
      
      if (!get_D(j)->get_Dirichlet_domain()) {
	cout << "Problem computing domain for manifold " <<(j+1)<< '.' << endl;
	break; 
      }
      
      O31_matrix mx;
      if (compute_isometry(domain, G, get_D(j)->domain, mx))
	cout << "Manifolds are isometric.\n";
      else
	cout << "Manifolds are not isometric.\n";
      
    }
    break;
    
  case PrintGOPairs:
    {
      double len;
      if (!get_input(s,"geodesic length ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &len) != 1) break; 
      
      double dist;
      if (!get_input(s,"orthodistance ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &dist) != 1) break; 
      
      print_go_pairs(domain, len, dist);
    }
    break;
    
  case GvSaveGeodesics:
    {
      vector<int> geodesics; 

      get_input(s,"geodesics to display ? ",-1); /* get a whole line of input */ 
      if (s.length()>0) geodesics = read_numbers(s.c_str()); 
      
      FILE* fp = fopen("snap.gv", "w");
      if (!fp) break; 
      
      gv_print_Dirichlet(fp, geodesics); 
      
      fclose(fp); 
    }
    break;
  case GvSaveOrtholines:
    {
      if (num_geodesics()==0) {
	cout << "Use \"print geodesics\" to compute some geodesics first.\n"; 
	break; 
      }

      int gnum; 
      if (!get_input(s, "which geodesic ? ")) break; 
      if (sscanf(s.c_str(), "%d", &gnum) != 1) break; 
      int n_radii; 
      if (!get_input(s, "how many radii ? ")) break; 
      if (sscanf(s.c_str(), "%d", &n_radii) != 1) break; 
	
      FILE* fp = fopen("snap.gv", "w");
      if (!fp) break; 
      gv_print_Dirichlet(fp, gnum, n_radii); 

      fclose(fp); 
    }
    break; 
    
  case TestNonArithmetic:
    {
      if (T && get_num_cusps(T->M())==num_filled_cusps(T->M())) {
	cout << "This test only applies to cusped manifolds.\n"; 
	break; 
      }
      compute_geodesics(0.98); 
      if (num_geodesics() > 0 && 
	  length_is_non_arithmetic(geodesic(0).length.real))
	cout << "Manifold is non-arithmetic.\n";
      else 
	cout << "Manifold may be arithmetic.\n"; 
    }
    break; 
  case HiddenSymmetryAction:
    hsymm_geodesic_action(); 
    break;
    
  default:
    env_T::process_event(what);
  }
}

void env_D::print_settings() const
{
  env_T::print_settings(); 

  // Use scientific notation for all epsilons. 
  ios_fmtflags old_fmt = cout.flags();
  cout.unsetf(ios::fixed);
  cout << "Dirichlet epsilon: " << dirichlet_epsilon << endl; 
  cout.flags(old_fmt); 

  printf("Options for print dirichlet: 0x%x\n", Dirichlet_print_options);
}

void env_D::hsymm_geodesic_action()
{
  Triangulation* m = T->M(); 

  if (num_filled_cusps(m) > 0) {
    cout << "Manifolds may not have any filled cusps\n"; 
    return; 
  }

  n_vector cusp_sizes = get_cusp_sizes(*T, true); 
  int gnum; 
  string s; 
  if (!get_input(s,"geodesic number ? ")) return;
  if (sscanf(s.c_str(), "%d", &gnum)!=1) return;

  bool copied; 
  if (!my_canonize(&m, cusp_sizes, copied)) {
    cout << "Canonize failed\n"; 
    return; 
  }

  vector<MoebiusTransformation> TR;
  bool or_pres;
  GraphNode* mg; 
  commensurabilities(m, m, TR, or_pres, &mg); 
  O31_matrix conj = domain->conjugacy; 
  O31_matrix c_inv = inverse(conj);
  O31_matrix M; 
  vector<O31_matrix> mats(TR.size()); 
  int i, n = TR.size(); 
  for (i=0; i<n; i++) {
    M = TR[i]; // Convert to O31. 
    mats[i] = c_inv*M*conj; 
  }
  print_symmetry_action(gnum, &mats, true);
  if (copied) free_triangulation(m); 
}

#if 0
static void print_length_commens(const i_triangulation* m, double radius)
{
  int i, n=m->num_geodesics();

  // First just count them. 
  for (i=0; i<n; i++) {
    if (m->geodesic(i).length.real > radius) break; 
  }
  n = i; 
}
#endif



