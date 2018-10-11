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
#include <ctype.h>
#include <fstream>
#include <cstdio>
#include "printable.hh"
#include "env_T.hh"
#include "snappea/unix_io.h"
#include "snappea/fundamental_group.h"
#include "snap_io.hh"
#include "terse_io.hh"
#include "linktable.hh"
#include "snappea/kernel.h"
#include "kernel_extras.hh"
#include "orthoangles.hh"
#include "helpers.hh"
#include "commens.hh"
#include "find_commens.hh"
#include "i_triangulation.hh"
#include "tilt_polytope.hh"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::sscanf;

#ifndef DATADIR
#define DATADIR "/usr/local/share/snap_data"
#endif

extern int show_holonomy_contribs;

bool env_T::simplify = true; 

static string default_dir = DATADIR; 

static const char* snl_text[] = {"none", "oriented", "full"};

enum {
  Comment=3,

  SaveMenu,
  SaveTest,
  ShowMenu,
  ShowSortedMenu,
  ShowDocumented,
  ShowUndocumented,
  ShowHelp, 
  ShowHistory,

  ReadFile, 
  ReadCensus, 
  ReadClosed,
  ReadManifold,
  ReadKnot,
  ReadDowker,
  ReadLink,
  ReadGroup, 

  ReadTerse,
  // ListTerse,
  CommitPendingKeys,

  ClearManifold, 

  SetManifold,
  ListManifolds,

  SetPath,
  SetShapeNormalization,
  SetSimplify,
  PrintSettings,

  PrintGroup, // group/triangulation required
  SaveGroup,
  SaveGenerators,

  Save, // triangulation required from here on down
  SaveTerse,
  IdentifyInCensus,
  CheckTwoBridge, 
  Surgery, 
  ConeSurgery,
  HolonomySolve,
  PrintSolutionType, 
  ClearSurgery, 

  Canonize,
  Randomize,
  Simplify,
  ReOrient,
  AcyclicEdges,
  FillCusps, 
  ChangeCurves,
  StandardizeCurves, 

  CheckCurves, 
  SplitEdge,

  DrillToStandard,

  PrintSnappeaSymmetries, 
  SnappeaIsometry,

  PrintHiddenSymmetries,
  PrintTiltMatrix,
  ComputeTiltPolytope,
  PrintTiltPolytope,
  SaveTiltPolytope,
  FindCommensurator, 
  BruteCommensurator, 
  Commensurable, 

  PrintHomology, 

  PrintVolume,
  PrintChernSimons,
  BarycentricSubdivision, 

  PrintShapes,
  PrintFullShapes, 
  LogShapes, 
  ShapeHistories,
  PrintTetrahedra, 
  PrintTetDomain, // Prints fund. dom. of ideal tetrahedra and images under face pairings.
  PrintOrthodistances, 
  PrintEdgeK,

  PrintFullEquations,
  PrintEquations,
  PrintPolyEquations,
  PrintFillingEquations, 
  PrintPolyFillingEquations, 
  VerifyEquationsApprox, 

  PrintCuspShapes,
  CuspDensity,
  CuspVolumes,
  PrintHolonomies,
  PrintCoreLength
};

void env_T::setup_menu()
{
  mm.add_item("?", Help); 
  mm.add_item("help", Help); 
  mm.add_item("show help", ShowHelp);
  mm.add_item("%", Comment); 
  mm.add_item("manifold", SetManifold); 
  mm.add_item("menu", ShowMenu);
  mm.add_item("show menu", ShowSortedMenu);
  mm.add_item("show documented", ShowDocumented);
  mm.add_item("show undocumented", ShowUndocumented);
  mm.add_item("history", ShowHistory);
  // mm.add_item("save test", SaveTest);
  // mm.add_item("save pending_keys", CommitPendingKeys);
  // mm.add_item("show terse", ListTerse); 
  mm.add_item("list", ListManifolds); 

  mm.add_item("read census", ReadCensus);
  mm.add_item("read closed", ReadClosed);
  mm.add_item("read file", ReadFile);
  mm.add_item("read manifold", ReadManifold);
  mm.add_item("read link", ReadLink);
  mm.add_item("read knot", ReadKnot);
  mm.add_item("read dowker", ReadDowker);
  mm.add_item("read terse", ReadTerse);
  mm.add_item("read group", ReadGroup);
  mm.add_item("clear manifold", ClearManifold); 
  mm.add_item("quit", Quit);

  mm.add_item("print settings", PrintSettings); 
  mm.add_item("set path", SetPath);
  mm.add_item("set shape_normalization", SetShapeNormalization); 
  mm.add_item("set simplify", SetSimplify); 

  mm.add_item("save menu", SaveMenu); 
  mm.add_item("save manifold", Save); 
  mm.add_item("save terse", SaveTerse); 

  mm.add_item("surgery", Surgery); 
  mm.add_item("solve", HolonomySolve);
  mm.add_item("clear surgery", ClearSurgery);
  mm.add_item("canonize", Canonize); 
  mm.add_item("randomize", Randomize); 
  mm.add_item("simplify", Simplify); 
  mm.add_item("reorient", ReOrient);
  mm.add_item("acyclic_edges", AcyclicEdges); 
  mm.add_item("fill", FillCusps); 
  mm.add_item("drill standard", DrillToStandard); 
  mm.add_item("check curves", CheckCurves); 
  mm.add_item("split_edge", SplitEdge);
  mm.add_item("change curves", ChangeCurves); 
  mm.add_item("standardize curves", StandardizeCurves); 

  mm.add_item("show tilt_matrix", PrintTiltMatrix); 
  mm.add_item("compute tilt_polytope", ComputeTiltPolytope); 
  mm.add_item("print tilt_polytope", PrintTiltPolytope); 
  mm.add_item("print hidden_symmetries", PrintHiddenSymmetries); 
  mm.add_item("find commensurator", FindCommensurator); 
  mm.add_item("save tilt_polytope", SaveTiltPolytope); 
  mm.add_item("brute_commensurator", BruteCommensurator); 

  mm.add_item("cone_surgery", ConeSurgery); // want "co" to mean compute.


  mm.add_item("print homology", PrintHomology); 

  mm.add_item("print volume", PrintVolume);
  mm.add_item("print chern_simons", PrintChernSimons); 
  mm.add_item("snappea volume", PrintVolume);
  mm.add_item("snappea chern_simons", PrintChernSimons); 

  // mm.add_item("bcs", BarycentricSubdivision);

  mm.add_item("print shapes", PrintShapes);
  mm.add_item("print all shapes", PrintFullShapes); 
  mm.add_item("print log_shapes", LogShapes); 
  mm.add_item("print histories", ShapeHistories); 
  mm.add_item("print tetrahedra", PrintTetrahedra);
  mm.add_item("print domain", PrintTetDomain); 

  mm.add_item("print gluing_equations", PrintEquations); 
  mm.add_item("print -p gluing_equations", PrintPolyEquations); 
  mm.add_item("print filling_equations", PrintFillingEquations); 
  mm.add_item("print -p filling_equations", PrintPolyFillingEquations); 
  mm.add_item("print full_equations", PrintFullEquations); 
  mm.add_item("verify approximate", VerifyEquationsApprox); 

  mm.add_item("print cusp shapes", PrintCuspShapes);
  mm.add_item("print cusp density", CuspDensity); 
  mm.add_item("print cusp volumes", CuspVolumes); 
  mm.add_item("print holonomies", PrintHolonomies);

  mm.add_item("print edge orthodistances", PrintOrthodistances); 
  mm.add_item("print edge k", PrintEdgeK); 

  mm.add_item("snappea symmetry_group", PrintSnappeaSymmetries); 

  mm.add_item("snappea isometry", SnappeaIsometry);
  mm.add_item("commensurable", Commensurable); 

  mm.add_item("print solution_type", PrintSolutionType);
  mm.add_item("print core length", PrintCoreLength);
  mm.add_item("print group", PrintGroup); 
  mm.add_item("save group", SaveGroup); 
  mm.add_item("save generators", SaveGenerators); 
  mm.add_item("identify", IdentifyInCensus); 
  mm.add_item("check two_bridge", CheckTwoBridge); 
}

void env_T::set_T(i_triangulation* m)
{
  delete T;
  T = m;
  T_changed();
}

void env_T::set_G(GroupPresentation* g)
{
  if (G) free_group_presentation(G); 
  G = g; 
}

bool env_T::get_group(bool required)
{
  if (T && !G) G = T->group(simplify); 
  if (!G && required) cout << "Function requires a group\n";
  return G!=0; 
}

void env_T::clear_group()
{
  set_G(0); 
}

i_triangulation* env_T::new_manifold(Triangulation* tri) const
{
  return new i_triangulation(tri); 
} 

void env_T::T_changed()
{
  set_G(0); 
  closed_num = 0; 
}

void env_T::T_shape_changed()
{
  set_G(0); 
  closed_num = 0; 
}

string env_T::name(int style) const
{
  if (!T) return gp_name; 
  return get_filling_name(T->M(), style);
}

int env_T::num_cusps() const
{
  if (T) {
    int nc = get_num_cusps(T->M()); 
    int nf = num_filled_cusps(T->M()); 
    return nc-nf;
  }
  if (!G) return 0;
  return fg_get_num_cusps(G); 
}

int env_T::choose_manifold(string const& prefix)
{
  int j; 
  string s; 
  string pr = prefix + "which manifold (1-10) ? ";
  if (!get_input(s, pr)) return -1; 
  sscanf(s.c_str(), "%d", &j); 
  j--; 
  if (j < 0) j = 0; 
  if (j > 9) j = 9; 
  return j;
}

int env_T::choose_save_manifold()
{
  string s; 
  int j=-1; 
  if (!get_input(s,"save result as which manifold (1-10) ? ")) return j; 
  if (sscanf(s.c_str(), "%d", &j)!=1) return -1; 
  j--; 
  if (j < 0) j = 0; 
  if (j > 9) j = 9; 
  env_T* alt = get_T(j); 
  if (alt->has_manifold()) {
    if (!ask("delete existing manifold ?", 0)) return -1; 
    alt->set_T(0);
  }
  return j; 
}

string choose_save_file(string save_dir, string def_name)
{
  // give user a chance to save somewhere other than the default location.
  string s;
  string uname = (save_dir == ".") ? def_name : (save_dir + "/" + def_name); 
  string pr = "Save in " + uname + "\n(hit return or type required filename) ? ";
  if (get_input(s, pr)) {
    string dir; 
    get_name_path(s, def_name, dir); 
    if (dir.size()) save_dir = dir; 
  }
  uname = (save_dir == ".") ? def_name : (save_dir + "/" + def_name); 

  // check if file already exists. 
  FILE* f = fopen(uname.c_str(), "r"); 
  if (f) { 
    fclose(f); 
    pr = "file " + uname + " exists! overwrite ? "; 
    if (!ask(pr, 0)) return ""; 
  }

  return uname; 
}

void env_T::set_manifold_from_surgery_description(string const& s)
{
  char census_char; 
  vector<double> sc(2);
  int N, mn;
  int o=0;
  while (s[o] && s[o]==' ') o++; 
  int narg = sscanf(s.c_str()+o, "%c%d( %lf , %lf )", 
		    &census_char, &mn, &sc[0], &sc[1]);
  if (narg!=2 && narg!=4) return; 
  
  /* get which census */ 
  if (census_char=='m') N = 5; 
  else if (census_char=='s') N = 6;
  else if (census_char=='v') N = 7;
  else return; 
  
  /* Read the manifold. */ 
  Triangulation* tri = read_census_manifold(path.c_str(), N, mn);
  if (!tri) {
    printf("Unable to read census manifold %d %d\n", N, mn);
    return; 
  }
  set_T(new_manifold(tri));
  
  /* Fill it if necessary. */
  if (narg==4) {
    T->do_Dehn_surgery(sc); 
  } 
}

void env_T::show_menu()
{
  string name = find_file_name(path, prog+"_menu"); 
  if (name.length()) { // If it has been saved in a file, "more" it. 
    ifstream file(name.c_str());
    string cmd = "less " + name; 
    system(cmd.c_str());
  } else {
    cout << mm;
  }
}

void env_T::show_help()
{
  string name = find_file_name(path, prog+"_help"); 
  if (!name.length()) { 
    cout << "Couldn't locate help file: " << prog << "_help.\n";
    return; 
  }
  ifstream file(name.c_str());
  string cmd = "less " + name; 
  system(cmd.c_str());
}

const char* solution_types[7] = {
  "not attempted", "geometric", "nongeometric", 
  "flat", "degenerate", "other", "no solution" 
};

static void print_abelian_group(AbelianGroup* hg)
{
  int i, n = hg->num_torsion_coefficients;
  for (i=0; i<n; i++) {
    if (hg->torsion_coefficients[i])
      cout << "Z/" << hg->torsion_coefficients[i];
    else 
      cout << "Z"; 
    if (i < n-1) cout << " + "; 
  }
}

void print_symmetry_group(SymmetryGroup* symm_gp)
{
  cout << "Order: " << symmetry_group_order(symm_gp) << endl; 
  AbelianGroup *abel=NULL;
  if (symmetry_group_is_abelian(symm_gp, &abel)) {
    cout << "Abelian group: ";
    print_abelian_group(abel); 
    cout << endl; 
  } else {
    cout << "Group is non-abelian.\n";
  }
  if (symmetry_group_is_dihedral(symm_gp)) {
    cout << "Group is dihedral.\n"; 
  }
}

static bool do_dehn_surgery(i_triangulation& mfld)
{
  char pr[100];
  string s; 
  int n = get_num_cusps(mfld.M());
  vector<double> coeffs(2*n);

  sprintf(pr, "Please enter %d pair(s) of surgery coefficients\n? ", n);
  get_input(s, pr, 0); 

  int cusp;
  double m, l; 

  for (cusp=0; cusp < n; ++cusp) {

    if (!get_input(s, "? (surgery coefficients) ", 2)) 
      return false;
    if (sscanf(s.c_str(), "%lf %lf", &m, &l) != 2)
      return false; 

    coeffs[2 * cusp] = m; 
    coeffs[2 * cusp + 1] = l; 
  }

  SolutionType sol = mfld.do_Dehn_surgery(coeffs);

  cout << "solution type: " << solution_types[sol] << endl; 

  return true; 
}

static bool do_cone_surgery(i_triangulation& mfld)
{
  char pr[100];
  string s; 
  int n = get_num_cusps(mfld.M());
  vector<double> coeffs(2*n);

  sprintf(pr, "Please enter %d rational \"slope order\" pair(s)\n? ", n);
  get_input(s, pr, 0); 

  int cusp;
  double m, l=1.0, o, d=1.0;

  for(cusp=0; cusp < n; ++cusp) {

    if (!get_input(s, "? slope", 1)) 
      return false;
    if (sscanf(s.c_str(), " %lf / %lf ", &m, &l) == 0)
      return false; 

    if (!get_input(s, "? order", 1)) 
      return false;
    if (sscanf(s.c_str(), " %lf / %lf ", &o, &d) == 0)
      return false; 

    if (d==0) {
      coeffs[2 * cusp] = 0.; 
      coeffs[2 * cusp + 1] = 0.; 
    } else {
      coeffs[2 * cusp] = m*o/d; 
      coeffs[2 * cusp + 1] = l*o/d; 
    }
  }

  SolutionType sol = mfld.do_Dehn_surgery(coeffs);

  cout << "solution type: " << solution_types[sol] << endl; 

  return true; 
}

static Complex normalize_shape(Complex s, int snl)
{
  if (!snl) return s;
  if (fabs(s.real-0.5)<1e-8 && fabs(fabs(s.imag)-ROOT_3_OVER_2)<1e-8)
    return s; 
  if (complex_modulus(s)-1.0 > -1e-8 || 
      complex_modulus(s - One)-1.0 > 1e-8) {
    if (s.real - .5 < 1e-8) 
      s = One/(One-s); 
    else 
      s = One - One/s;
  }
  if (snl<2) return s; 
  if (s.real - .5 > 1e-8)
    s.real = 1.0 - s.real; 
  return s;
}

static void print_shapes(const vector<Complex>& vc, int snl=0)
{
  vector<Complex>::const_iterator vcci;
  cout << "[";
  for (vcci = vc.begin(); vcci != vc.end();) {
    print(normalize_shape(*vcci,snl)); 
    vcci++;
    if (vcci != vc.end()) cout << ", "; 
  }
  cout << "]"; 
}

void print_shape_histories(Triangulation* manifold)
{
  int n = get_num_tetrahedra(manifold);

  Complex* shapes = new Complex [n]; 

  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    shapes[t->index] = t->shape[filled]->cwl[ultimate][0].rect;

  int j, k = get_num_tetrahedra(manifold); 
  ShapeInversion **histories = shape_histories(manifold); 
  ShapeInversion *inv; 
  for (j=0; j < k; j++) {
    printf("shape(%d) = ", j); 
    print(shapes[j]);
    printf(" ("); 
    inv = histories[j];
    while (inv) {
      printf("%d", inv->wide_angle); 
      inv = inv->next; 
      if (inv) printf(" "); 
    }
    printf(")\n"); 
  }
  delete [] shapes; 
  delete [] histories; 
}

static void verify_equations_approx(int_matrix const& eqns, 
				    vector<Complex> const& shapes)
{
  int ns = shapes.size(); 
  if (eqns.cols != ns + 1) {
    printf("Oops, mismatch between equations and number of tetrahedra!\n"); 
    return; 
  }

  int w = field_width(eqns); 
  int r, i; 
  Complex z, PiI(0.,PI); 
  for (r=0; r<eqns.rows; r++) {
    pretty_print(eqns, w, r); 
    printf(" -> ");

    z = Zero; 
    for (i=0; i<ns; i++) {
      z += double(eqns[r][i]) * shapes[i];
    }
    z += double(eqns[r][ns]) * PiI; 
 
    fwprint(stdout, z, 8); 
    printf("\n"); 
  }
}

static void print_holonomies(Triangulation* manifold, vector<Complex> const& tgt)
{
  int n = get_num_cusps(manifold);
  int cusp;
  Complex mh, lh; 

  print_holonomies(manifold); 

  cout << "TH=[";
  for (cusp=0; cusp < n; ++cusp) {
    if (cusp > 0) cout << ", ";
    cout << "[" << tgt[2*cusp] << ", " << tgt[2*cusp+1] << "]";
  }
  cout << "]\n"; 

  cout << "TE=[";
  for (cusp=0; cusp < n; ++cusp) {
    if (cusp > 0) cout << ", ";
    get_holonomy(manifold, cusp, &mh, &lh, 0, 0);
    cout << "[" << 0.5*complex_exp(mh/2.0)*tgt[2*cusp] << ", " << 0.5*complex_exp(lh/2.0)*tgt[2*cusp+1] << "]";
  }
  cout << "]\n"; 
}

static bool solve_for_holonomies(i_triangulation& mfld)
{
  char pr[100];
  string s, t; 
  int n = get_num_cusps(mfld.M());
  vector<Complex> holo(n);
  vector<bool> hism(n);

  if (n==1) 
    sprintf(pr, "Please enter m[e] or l[e] followed by a complex holonomy\n? ");
  else 
    sprintf(pr, "Please enter %d holonomies (m[e] or l[e] followed by complex holonomy)\n? ", n);
  get_input(s, pr, 0); 

  int cusp;
  Complex z, mh, lh; 

  for (cusp=0; cusp < n; ++cusp) {

    if (!get_input(s, "? (m[e] or l[e]) ", 1)) return false; 
    if (s[0] == 'm') hism[cusp] = true; 
    else if (s[0] == 'l') hism[cusp] = false; 
    else return false;

    if (!get_input(t, "? (holonomy) ", 1)) 
      return false;
    if (!complex_from_string(t.c_str(), z))
      return false; 

    // If user specified an eigenvalue, we need to take 2 * its logarithm
    // before calling do_Dehn_surgery. The branch is chosen by assuming 
    // that the new holonomy is close to the previous one. 
    if (s.length() > 1 && s[1]=='e') {
      get_holonomy(mfld.M(), cusp, &mh, &lh, 0, 0); 
      z = 2.0 * complex_log(z, ((s[0]=='m') ? mh.imag : lh.imag)/2.0);
    }
    holo[cusp] = z; 
  }

  vector<Complex> tgt(2*n, Zero);

  SolutionType sol = mfld.do_Dehn_surgery(holo, hism, tgt);

  cout << "solution type: " << solution_types[sol] << endl; 
  print_holonomies(mfld.M(), tgt); 
  return true; 
}

static int intmatrix22_from_string(const char* str, MatrixInt22 mx)
{
  char delim;
  const char* s = str;
  int i, j; 

  while (*s && *s!='[') s++;
  if (!*s) return 0; 
  s++;

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      if (sscanf(s, " %d ", &mx[i][j]) != 1) return 0; 
      delim = (j==0) ? ',' : ((i==0) ? ';' : ']');
      while (*s && *s!=delim) s++;
      if (!*s) return 0; 
      s++;
    }
  }
  return s - str; 
}

static bool get_basis_change_matrices(int n, const char* s, MatrixInt22 mats[])
{
  char delim;
  int i, nc; 
  
  while (*s && *s!='[') s++;
  if (!*s) return false; 
  s++;
  
  for (i=0; i<n; i++) {
    nc = intmatrix22_from_string(s, mats[i]);
    if (!nc) return false; 
    s += nc; 

    delim = (i<n-1) ? ',' : ']';
    while (*s && *s!=delim) s++;
    if (!*s) return false; 
    s++;
  }
  return true; 
}

static vector<int> get_cusp_ties(int nc)
{
  string s; 
  vector<int> ties(nc); 
  if (!get_input(s, "cusp tie classes [0,...] ? ") || 
      !scan_vector(s.c_str(), ties) || 
      ties.size()!=nc ||
      !validate_class_spec(ties)) {
    int i; 
    ties.resize(nc); 
    for (i=0; i<nc; i++) ties[i]=i; 
  } 
  return ties; 
}

static vector<int> get_area_vector(int nc)
{
  vector<int> areas(nc, 1); 

  string s; 
  if (!get_input(s, "area ratios (a1:a2:...) ? ")) 
    return areas; // equal area one. 

  // Read area multipliers. 
  int i=0; 
  const char* cs = s.c_str(); 
  while (*cs && i<nc) { 
    if (sscanf(cs, "%d", &areas[i])!=1) break;
    i++;
    while (*cs) { // Skip past next ':'
      if (*cs==':') { cs++; break; }
      cs++;
    }
  }
  return areas; 
}

// Get a : separated list of c doubles from the user,
// normalized such that the largest is 1.0. 

n_vector get_cusp_sizes(int c)
{
  string s; 
  if (!get_input(s, "area ratios (a1:a2:...) ? ")) 
    return n_vector(c,1.0); // default is equal areas. 

  n_vector z;
  int i; 

  // check if user wants to give sizes instead. 
  bool areas=true; 
  if (s==string("-s")) {
    areas=false; 
    if (!get_input(s, "cusp cross section sizes (s1:s2:...) ? ")) return z;
  } 

  n_vector size_vector(c); 

  // Read area ratios. 
  i=scan_ratio_vector(s.c_str(),size_vector); 
  // Set any unassigned ratios to 1.0. 
  for (;i<c;i++) size_vector[i]=1.0; 
  
  // Normalize so that maximum multiplier is 1.0; 
  double max_mul=1e-3, mul;
  for (i=0; i<c; i++) {
    mul = size_vector[i];
    if (mul <= 0.) mul = size_vector[i] = 1e-3;
    if (mul > max_mul) max_mul=mul; 
  }
  for (i=0; i<c; i++) size_vector[i] /= max_mul; 

  if (areas) for (i=0; i<c; i++) size_vector[i] = sqrt(size_vector[i]); 

  return size_vector; 
}


n_vector get_cusp_sizes(i_triangulation const& m, bool always)
{
  int c = get_num_cusps(m.M()); 
  tilt_polytope *TP = m.tp(); 

  string s; 
  if (c<2 && !always) return n_vector(c,1.0); 
  if (!get_input(s, "area ratios (a1:a2:...) ? ")) 
    return n_vector(c,1.0); // default is equal areas. 

  n_vector z;
  int i; 

  // check if user wants to give sizes instead. 
  bool areas=true; 
  if (s==string("-s")) {
    areas=false; 
    if (!get_input(s, "cusp cross section sizes (s1:s2:...) ? ")) return z;
  } else if (s==string("-t")) {
    if (!TP) { cout << "please compute the tilt polytope first\n"; return z; }
    if (!get_input(s, "tilt polytope vertices (v1,v2,...) ? ")) return z;    
    
    // get size vector for a tilt polytope face
    n_vector size_vector(c,0.);
    const char* cs = s.c_str(); 
    int vnum, n=0; 
    polytope<VS>::vertex_cp vp; 
    if (*cs=='(') cs++; //')'
    while (*cs) { 
      if (sscanf(cs, "%d", &vnum)!=1) break;
      for (vp = TP->P.VL.begin(); vp!=TP->P.VL.end(); ++vp) {
	if (vp->index==vnum) {
	  size_vector += vp->V; 
	  ++n;
	}
      }
      while (*cs) { // Skip past next ','
	if (*cs==',') { cs++; break; }
	cs++;
      }
    }
    if (on_boundary(size_vector)) {
      cout << "bad cusp size vector " << size_vector << endl; 
      return z; 
    }
    size_vector /= n; 

    return size_vector; 
  }

  n_vector size_vector(c); 

  // Read area ratios. 
  i=scan_ratio_vector(s.c_str(),size_vector); 
  // Set any unassigned ratios to 1.0. 
  for (;i<c;i++) size_vector[i]=1.0; 
  
  // Normalize so that maximum multiplier is 1.0; 
  double max_mul=1e-3, mul;
  for (i=0; i<c; i++) {
    mul = size_vector[i];
    if (mul <= 0.) mul = size_vector[i] = 1e-3;
    if (mul > max_mul) max_mul=mul; 
  }
  for (i=0; i<c; i++) size_vector[i] /= max_mul; 

  if (areas) for (i=0; i<c; i++) size_vector[i] = sqrt(size_vector[i]); 

  return size_vector; 
}



Triangulation* get_DT_knot()
{
  string s; 
  if (!get_input(s, "Alphabetic/numeric Dowker notation ? ")) return 0; 

  Triangulation* tri; 
  int j; 
  string name; 

  if (sscanf(s.c_str(), "%d", &j)==0) { // Alphabetic
    tri = DT_alpha_to_triangulation((char*)s.c_str()); 
    if (tri) set_triangulation_name(tri, (char*)s.c_str()); 
    return tri; 

  } else {

    // Numeric code supplied. j = number of crossings. 
    if (j < 3) {
      cout << "Invalid number of crossings: " << j << endl; 
      return 0; 
    } 
    
    int c, i; 
    int *code = new int [j]; 
    for (i = 0; i < j; i++) {
      if (!get_input(s, "next crossing ? ")) break; 
      if (!sscanf(s.c_str(), "%d", &c)) break; 
      code[i] = c; 
      name += char((c >= 0) ? 'a' + (c/2) - 1 : 'A' - (c/2) - 1); 
    }
    if (i < j) {
      delete [] code; 
      cout << "Incomplete code.\n";
      return 0; 
    }
    tri = DT_int_to_triangulation(j, code); 
    delete [] code; 
  }

  if (tri) set_triangulation_name(tri, (char*)name.c_str()); 
  else 
    cout << "Sorry, unable to compute a triangulation for that code.\n"; 

  return tri; 
}

static string dbdir(string path)
{
  char* home = getenv("HOME");
  if (home) path = string(home)+" "+path; 
  vector<string> tpath;
  split(path, tpath); 

  // Look in the path first. 
  int i, n=tpath.size(); 
  for (i=0; i<n; ++i) {
    if (file_exists(tpath[i]+"/terse/terse_by_volume"))
      return tpath[i]+"/terse";
  }

  // If not found, assume $HOME/terse
  string dir("terse");
  if (home) dir = string(home) + "/terse";
  return dir;
}

static void print_poly(const int* coeff, int n)
{
  string rhs;
  char buf[256];
  bool is_one=true;
  int i, c;
  for (i=0; i<n; ++i) {
    c = coeff[i];
    if (c==0) continue;
    if (c < 0) {
      if (c==-1) sprintf(buf, "z%d", i);
      else sprintf(buf, "z%d^%d", i, -c);
      if (rhs.size()) rhs += "*";
      rhs += buf;
    } else {
      if (!is_one) cout << "*";
      is_one = false;
      if (c==1) sprintf(buf, "z%d", i);
      else sprintf(buf, "z%d^%d", i, c);
      cout << buf;
    }
  }
  for (i=0; i<n; ++i) {
    c = coeff[i+n]; 
    if (c==0) continue;
    if (c < 0) {
      if (c==-1) sprintf(buf, "(1-z%d)", i);
      else sprintf(buf, "(1-z%d)^%d", i, -c);
      if (rhs.size()) rhs += "*";
      rhs += buf;
    } else {
      if (!is_one) cout << "*";
      is_one = false;
      if (c==1) sprintf(buf, "(1-z%d)", i);
      else sprintf(buf, "(1-z%d)^%d", i, c);
      cout << buf;
    }
  }
  if (is_one) cout << "1"; 
  if (coeff[2*n]%2) cout << " + ";
  else cout << " - ";

  if (rhs.size()==0) cout << "1\n"; 
  else cout << rhs << endl;
}

text_menu env_T::mm;
string env_T::path;
string env_T::prog; 
int env_T::cm = 0; 
bool env_T::batchmode = false; 
bool env_T::echomode  = false; 
int env_T::snl = 0; 
int env_T::tp_max = 2000; 

void env_T::init(string prg)
{
  prog = prg; 
  const char *str = getenv("SNAPPATH");
  if (!str) path = string(". "+default_dir);
  else path = string(str); 
}

env_T::~env_T()
{
  delete T;
  if (G) free_group_presentation(G); 
}

void env_T::set_options(int argc, char* argv[])
{
  int j;

  for (j=1; j<argc; j++) {
    if (strcmp(argv[j], "-be")==0) { batchmode = true; echomode = true; }
    if (strcmp(argv[j], "-b")==0) batchmode = true; 
  }
}

int env_T::get_event()
{
  int what; 
  char* ln; 
  char prompt[256];

  if (T || G) 
    sprintf(prompt, "%d. %s: ", cm+1, name().c_str()); 
  else 
    sprintf(prompt, "%d. : ", cm+1);
  
  what = mm.read_input(batchmode ? "\0" : prompt, &ln); 

  // this assumes one instruction per line (in echomode)
  // .. if not the whole input line will be repeated. 
  if (echomode) {
    if (what!=Nothing)
      printf("%s%s\n", prompt, ln); 
    else printf("\n"); 
  }

  return what;
}

void env_T::validate_event(int& what)
{
  if (what > PrintCoreLength) return; 

  // Check they have a manifold to do it on. 
  if (what >= Save && !T) {
    cout << "This function requires a triangulation.\n";
    what = Nothing; 
  }

  if (what >= PrintGroup && !(T||G)) {
    cout << "This function requires a manifold.\n";
    what = Nothing; 
  }
}

void env_T::print_settings() const
{
  cout << "Path: " << path << endl; 
  cout << "Shape normalization: " << snl_text[snl] << endl;
  cout << "Tilt polytope face limit: " << tp_max << endl; 
}

bool env_T::read_group()
{
  string s; 
  if (!get_input(s, "Group presentation file ? ")) return false;
  FILE* fp = fopen(s.c_str(), "r");
  if (!fp) {
    cout << "Unable to open file " << s << " for reading\n"; 
    return false; 
  }
  GroupPresentation* G0=read_presentation(fp);
  fclose(fp); 
  if (!G0) return false; 
  set_T(0); 
  set_G(G0);

  string file, path, base, ext;
  get_name_path(s, file, path);
  get_base_ext(file, base, ext); 
  if (ext=="gens" || ext=="grp") gp_name=base;
  else gp_name = file; 

  return true; 
}


void env_T::process_event(int what)
{
  string s; 
  int j; 

  // Do it!
  switch (what) {
  case Nothing:
    break;
  case Comment:
    if (!get_input(s, "?", -1)) break; 
    if (echomode) printf("%% %s\n", s.c_str()); 
    break;
  case Help:
    {
      int t;
      string help_topic = mm.expand_input("<" + prog + "-command> ? ", t); 
      if (!help_topic.length()) {
	t = 1; help_topic = "help"; 
      }
      if (t==-1) {
	cout << '"' << help_topic << '"' << " is not a valid command.\n"; 
	cout << "Type \"menu\" for a list of valid " << prog << " commands.\n"; 
      } else if (t > 0) {
	print_help(help_topic, prog+"_help", path); 
      }
    }
    break;

  case SaveMenu:
    {
      string menu_file = prog + "_menu"; 
      string filename = choose_save_file(".", menu_file); 
      if (!filename.length()) break; 
      ofstream file(menu_file.c_str()); 
      if (!file) { 
	cout << "Sorry, couldn't open " << filename << " for writing.\n"; break; }
      file << mm; 
    }
    break; 
  case SaveTest:
    {
      string test_file = find_unused_name("test", prog+"_", ".in");
      test_file = choose_save_file("test", test_file); 
      if (!test_file.length()) break; 
      ofstream tf(test_file.c_str());
      if (!tf) {
	cout << "Couldn't open " << test_file << " for writing.\n"; 
	break; 
      }
      tf << "% PROG: " << prog << " -be\n\n"; 
      print_history(tf);
    }
    break;
  case ShowMenu:
    show_menu(); 
    break; 
  case ShowHelp:
    show_help(); 
    break; 
  case ShowSortedMenu:
    {
      vector<string> menu_items; 
      mm.get_menu(menu_items); 
      int i, n=menu_items.size();
      for (i=0; i<n; i++) cout << menu_items[i] << endl; 
    }
    break; 
  case ShowDocumented:
    {
      vector<string> doc_items; 
      get_documented_topics(prog+"_help", path, doc_items); 
      int i, n=doc_items.size();
      for (i=0; i<n; i++) cout << doc_items[i] << endl; 
    }
    break; 
  case ShowUndocumented:
    {
      string item; 
      vector<string> menu_items, doc_items; 
      get_documented_topics(prog+"_help", path, doc_items); 
      mm.get_menu(menu_items); 
      int i, j=0, n=menu_items.size();
      for (i=0; i<n; i++) {
	item = menu_items[i]; 
	while (j < doc_items.size() && doc_items[j] < item) j++; 
	if (j==doc_items.size()) {
	  cout << item << endl; 
	  continue; 
	}
	if (doc_items[j].substr(0,item.length())==item) continue; 
	cout << item << endl; 
      }
    }
    break; 
  case ShowHistory:
    print_history(); 
    break;

  case SetManifold:
    if (!get_input(s,"which manifold (1-10) ? ")) break; 
    sscanf(s.c_str(), "%d", &cm); 
    cm--; 
    if (cm < 0) cm = 0; 
    if (cm > 9) cm = 9; 
    break;
  case ListManifolds:
    for (j=0; j<10; j++)
      if (get_T(j)->has_manifold()) 
	printf("%d = %s\n", j+1, get_T(j)->name().c_str());
    break; 
    
  case SetPath:
    if (!get_input(s, "blank separated list of directories ? ", -1)) break;
    path = s; 
    break;
  case SetShapeNormalization:
    if (!get_input(s, "Shape normalization 0=none, 1=oriented, 2=full ? ")) 
      break; 
    sscanf(s.c_str(), "%d", &snl);
    if (snl < 0) snl = 0; 
    if (snl > 2) snl = 2; 
    break; 
  case SetSimplify:
    simplify = ask("simplify fundamental group ?", 1);
    break;
  case PrintSettings:
    print_settings(); 
    break;
    
  case ReadCensus: 
    {
      int N, i;
      if (!get_input(s, "Which census, manifold number? ", 2)) break; 
      if (sscanf(s.c_str(), "%d %d", &N, &i) != 2) break; 
      
      Triangulation* tri = read_census_manifold(path.c_str(), N, i);
      if (!tri) {
	cout << "Unable to read " << N << "-census manifold " << i << endl;
	break;
      }
      set_T(new_manifold(tri));
    }
    break;
  case ReadFile:
    {
      if (!get_input(s, "Which manifold? ")) break;
      string file = find_file_name(path, s, "manifolds"); 
      Triangulation* tri = read_manifold_file(path, s); 
      if (!tri) break; 

      set_T(new_manifold(tri));
    }
    break;

  case ReadKnot:
    {
      Triangulation* tri = get_DT_knot();
      if (!tri) break;
      set_T(new_manifold(tri));
    }
    break;

  case ReadDowker:
    if (!get_input(s, "Dowker/Thistlethwaite code ? ")) break; 
    {
      Triangulation* tri=DT2Triangulation(s); 
      if (!tri) break; 
      set_T(new_manifold(tri));
    }
    break;

  case ReadLink:
    {
      int ncross; 
      if (!get_input(s, "number of crossings ? ")) break; 
      if (sscanf(s.c_str(), "%d", &ncross)!=1) break; 
      
      char alt; 
      if (!get_input(s, "a or n (alternating or non-alternating) ? ")) break; 
      if (sscanf(s.c_str(), "%c", &alt)!=1) break; 
      
      int index; 
      if (!get_input(s, "index of link in table ? ")) break; 
      if (sscanf(s.c_str(), "%d", &index)!=1) break; 
      
      string code = lookup_link_dtcode(path, ncross, alt, index); 
      if (!code.length()) break; // message already given. 
      
      Triangulation* tri=DT2Triangulation(code); 
      if (!tri) break; 
      cout << "Dowker/Thistlethwaite: " << code << endl; 

      char buf[50];
      sprintf(buf, "%d%c%d", ncross, alt, index); 
      set_triangulation_name(tri, buf);

      set_T(new_manifold(tri));
    }
    break;

  case ReadManifold:
    {
      if (!get_input(s, "surgery description? ")) break; 
      set_manifold_from_surgery_description(s);
    }
    break;
  case ReadClosed:
    {
      int index; 
      if (!get_input(s, "index in closed census ? ")) break; 
      if (sscanf(s.c_str(), "%d", &index)!=1) break;

      FILE* fp = open_data_file(path, "ClosedManifolds"); 
      if (!fp) break; 
      char* line = get_line(fp, index, 31); // data starts after line 31
      if (!line) {
	cout << "index must be in range 1-11031\n";
	break;
      }
      // surgery description is in cols 72-83.
      line[84] = '\0';
      set_manifold_from_surgery_description(line+72);
      closed_num = index; 
    }
    break;
  case ReadTerse:
    {
      string dir = dbdir(path);

      int n_tet; 
      if (!get_input(s, "number of tetrahedra ? ")) break; 
      if (sscanf(s.c_str(), "%d", &n_tet)!=1) break; 
      
      int index; 
      if (!get_input(s, "index of triangulation ? ")) break; 
      if (sscanf(s.c_str(), "%d", &index)!=1) break; 

      int maxind = num_terse(dir.c_str(), n_tet);
      if (maxind < 1) {
	cout << "You do not appear to have any " << n_tet
	     << "-tet triangulations saved.\n";
	break; 
      }
      if (index >= maxind || index < 0) {
	cout << "Index must be in range 0-" << maxind-1 << endl; 
	break; 
      }

      Triangulation* t = read_terse(dir.c_str(), n_tet, index); 
      if (!t) {
	cout << "There was some problem reading the triangulation\n"; 
	break; 
      }

      set_T(new_manifold(t)); 
    }
    break;
  case ReadGroup:
    if (read_group()) 
      print(G, true); 
    closed_num = 0; 
    break;

  case ClearManifold:
    set_T(0);
    break; 
    
  case Save: // manifold required from here on down
    {
      string name;
      if (num_filled_cusps(T->M()) == 0) 
	name = get_triangulation_name(T->M()); 
      else 
	name = get_filling_name(T->M(),1); 
      string file = choose_save_file(".", name); 
      if (!file.length()) break; 
      T->save_manifold_file(file.c_str()); 
    }
    break; 
  case SaveTerse:
    {
      int index, nt;
      manifold_database DB(dbdir(path).c_str());
      if (!DB.insert(T->M(), index, nt)) {
	if (index < 0) 
	  cout << "There was a problem saving the triangulation.\n";
	else 
	  cout << "Manifold already exists with index " <<
	    nt << '.' << index << endl;
	break; 
      }
      cout << "Saved triangulation has index " << 
	nt << '.' << index << endl; 
    }
    break; 
    
  case IdentifyInCensus:
    {
      int C, N;
      if (find_census_manifold(path.c_str(), T->M(), C, N)) {
	cout << "Census: " << C << " manifold: " << N << endl; 
      } else {
	cout << "Manifold not found in cusped census\n"; 
      }

      string code; 
      code = find_in_linktable(path, T->M());
      if (code.length()) {
	cout << "Dowker code of link " << code << endl; 
      } else {
	cout << "Manifold not found in link tables\n";
      }

      manifold_database DB(dbdir(path).c_str());
      int index, n_tet; 
      if (DB.locate(T->M(), index, n_tet)) {
	cout << "Manifold in terse database " << n_tet<<'.'<<index<<endl;
      }
    }
    break; 
#if 0
  case ListTerse:
    {
      manifold_database DB(dbdir(path).c_str()); 
      DB.print();
    }
    break;
#endif
  case CommitPendingKeys:
    {
      manifold_database DB(dbdir(path).c_str()); 
      DB.save_pending_keys();
    }
    break;

  case CheckTwoBridge:
    {
      Boolean is_two_bridge;
      long p, q; 
      two_bridge(T->M(), &is_two_bridge, &p, &q); 
      if (is_two_bridge) 
	cout << "Two bridge knot " << p << '/' << q << endl; 
      else 
	cout << "Does not appear to be a 2-bridge knot\n"; 
    }
    break;

  case Surgery: 
    do_dehn_surgery(*T);
    T_shape_changed();
    break;
  case ConeSurgery:
    do_cone_surgery(*T);
    T_shape_changed();
    break;

  case HolonomySolve: 
    solve_for_holonomies(*T);
    T_shape_changed();
    break;
    
  case PrintSolutionType:
    {
      int sol = get_filled_solution_type(T->M());
      cout << "solution type: " << solution_types[sol] << endl; 
    }
    break;
  case ClearSurgery:
    T->find_complete_hyperbolic_structure(); 
    T_shape_changed();
    break;
    
  case ReOrient:
    T->reorient();
    T_changed(); 
    break; 

  case AcyclicEdges:
    {
      vector<int> e_or; 
      if (seek_acyclic_edge_orientations(T->M(), e_or, true)) {
	cout << "Manifold has acyclic edge orientations\n";
	print_vector(cout, e_or); cout << endl; 
      }
#if 0
      else
	cout << "Manifold has no acyclic edge orientations\n";
#endif
    }
    break; 

  case Canonize:
    {
      n_vector size_vector = get_cusp_sizes(*T);
      if (!size_vector.dim) break; 
      if (!T->canonize(size_vector.V)) {
	cout << "Canonize failed\n";
	break; 
      }
      T_changed(); 
      if (!is_canonical_triangulation(T->M())) {
	cout << "Contains subdivided non-tetrahedral cells\n"; 
	choose_generators(T->M(), TRUE, FALSE); // Computes vertices. 
	check_no_generator_transparent(T->M());

	Triangulation* mc;
	copy_triangulation(T->M(), &mc); 
	canonical_retriangulation(mc); 
	print_ideal_cells(mc);
	free_triangulation(mc); 
      }
    }
    break;
  case FillCusps:
    T->fill(); 
    T_changed();
    break;

  case DrillToStandard:
    {
      j = choose_save_manifold();
      if (j<0) break; 

      load_standard_set(path,1); /* ensure that these are present */ 
      if (!standard_set_present) {
	printf("Problem loading standard set.\n"); break; }

      Triangulation* result; 
      int which_std = get_filled_standard(T->M(), &result, 1); 

      if (which_std == -1) {
	printf("Unable to find an ancestor manifold in the standard set.\n");
	break; 
      }
	  
      // Get the eta fudge for new_m with its standard peripheral curves. 
      get_T(j)->set_T(new_manifold(result));
    }
    break;

  case PrintTiltMatrix:
    {
      int nf = 2*get_num_tetrahedra(T->M()); 
      int nc = get_num_cusps(T->M()); 
      int i, j; 
      double* tilt_mx = face_tilt_matrix(T->M()); 
      for (i=0; i<nf; i++) {
	for (j=0; j<nc; j++) {
	  fwprint(cout, tilt_mx[i*nc + j], 8);  
	  if (j<nc-1) cout << ' '; 
	}
	cout << endl; 
      }
      my_free_array(tilt_mx); 
    }
    break; 
  case PrintTiltPolytope:
    {
      if (T->tp()) 
	T->tp()->print(); 
      else cout << "Please call \"compute tilt_polytope\" first!\n"; 
    }
    break;
  case ComputeTiltPolytope:
    {
      vector<int> ties = get_cusp_ties(get_num_cusps(T->M())); 
      n_vector tie_ratios = get_cusp_sizes(*T);
      T->clear_tilt_polytope(); 
      if (!T->compute_tilt_polytope(ties, tie_ratios, tp_max))
	cout << "Polytope has more than " << tp_max << " faces!\n"; 
      else 
	cout << "Tilt polytope has " << (T->tp()->FL.size()) << " faces\n";
    }
    break;
  case SaveTiltPolytope:
    {
      ps_picfile pic = ps_picfile("tilt.ps"); 
      pic.set_scale(72.0);
      T->tp()->save_picture(pic);
      pic.print_text(Complex(-2.,-.4), get_triangulation_name(T->M())); 
      pic.close(); 
    }
    break;
  case PrintHiddenSymmetries:
    {
      n_vector csm = get_cusp_sizes(*T, true); 
      print_hidden_symmetries(*T, csm); 
    }
    break; 
  case FindCommensurator:
    find_commensurator(*T); 
    break; 
  case BruteCommensurator:
    {
      int n=0, N=0; 
      if (!get_input(s, "area vector sum(s) to check ? ")) return; 
      if (sscanf(s.c_str(), "-%d", &N)==1) {}
      else if (sscanf(s.c_str(), "%d-%d", &n, &N)==0) return;

      int c = get_num_cusps(T->M());
      vector<int> ties = get_cusp_ties(c); 
      vector<int> areas = get_area_vector(c);  
      brute_force_hsymms(*T, n, N, ties, areas); 
    }
    break;
  case Commensurable:
    {
      j = choose_manifold(); 
      if (!get_TT(j)) {
	cout << "There is no triangulation " << (j+1) << endl; 
	break; 
      }
      n_vector csm = get_cusp_sizes(*T,true); 
      n_vector csn = get_cusp_sizes(*get_TT(j),true); 
      
      check_commensurability(*T, csm, *get_TT(j), csn);
    }
    break;
    
  case Simplify:
    T->simplify();
    break;
    
  case Randomize:
    T->randomize();
    break;
    
  case PrintShapes:
    {
      vector<Complex> shapes = T->shapes(); 
      if (snl) {
	if (snl<2) cout << "Oriented congruence classes\n";
	else cout << "Congruence classes\n";
      }
      print_shapes(shapes, snl);
      cout << "\n"; 
    }
    break;
  case PrintFullShapes:
    print_full_shapes(T->M()); 
    break; 

  case LogShapes:
    {
      vector<Complex> shapes = get_log_shapes(T->M()); 
      print_shapes(shapes);
      cout << "\n"; 
    }
    break;

  case ShapeHistories:
    print_shape_histories(T->M()); 
    break;

  case PrintTetrahedra:
    print_tetrahedra(T->M()); 
    break;
    
  case PrintTetDomain:
    choose_generators(T->M(), TRUE, FALSE); // Computes vertices. 
    T->print_domain(); 
    break; 
  case PrintEquations:
    {
      int_matrix e_eqns, c_eqns; 
      get_edge_equations(T->M(),e_eqns); 
      get_cusp_equations(T->M(),c_eqns); 

      int w = field_width(c_eqns); 
      if (field_width(e_eqns) > w) w = field_width(e_eqns);

      cout << "Edge equations:\n"; 
      pretty_print(e_eqns, w); 
      cout << "Cusp equations:\n"; 
      pretty_print(c_eqns, w); 
    }
    break;
  case PrintPolyEquations:
    {
      int_matrix e_eqns, c_eqns; 
      get_edge_equations(T->M(),e_eqns); 
      get_cusp_equations(T->M(),c_eqns); 

      int i, r, n = get_num_tetrahedra(T->M());

      cout << "Edge equations:\n"; 
      r = e_eqns.rows;
      for (i=0; i<r; ++i) print_poly(e_eqns[i], n);

      cout << "Cusp equations:\n"; 
      r = c_eqns.rows;
      for (i=0; i<r; ++i) print_poly(c_eqns[i], n);
    }
    break;
  case PrintFillingEquations:
    {
      int_matrix e_eqns, f_eqns; 
      get_edge_equations(T->M(),e_eqns); 
      get_filling_equations(T->M(),f_eqns); 

      int w = field_width(f_eqns); 
      if (field_width(e_eqns) > w) w = field_width(e_eqns);

      cout << "Edge equations:\n"; 
      pretty_print(e_eqns, w); 
      cout << "Filling equations:\n"; 
      pretty_print(f_eqns, w); 
    }
    break;
  case PrintPolyFillingEquations:
    {
      int_matrix e_eqns, f_eqns; 
      get_edge_equations(T->M(),e_eqns); 
      get_filling_equations(T->M(),f_eqns); 

      int i, r, n = get_num_tetrahedra(T->M());

      cout << "Edge equations:\n"; 
      r = e_eqns.rows;
      for (i=0; i<r; ++i) print_poly(e_eqns[i], n);

      cout << "Cusp equations:\n"; 
      r = f_eqns.rows;
      for (i=0; i<r; ++i) print_poly(f_eqns[i], n);
    }
    break;
  case PrintFullEquations:
    {
      int_matrix e_eqns, c_eqns; 
      get_full_edge_equations(T->M(),e_eqns); 
      get_full_cusp_equations(T->M(),c_eqns); 

      int w = field_width(c_eqns); 
      if (field_width(e_eqns) > w) w = field_width(e_eqns);

      cout << "Edge equations:\n"; 
      pretty_print(e_eqns, w); 
      cout << "Cusp equations:\n"; 
      pretty_print(c_eqns, w); 
    }
    break;
  case VerifyEquationsApprox:
    {
      int_matrix eqns; 
      vector<Complex> shapes = get_log_shapes(T->M()); 

      cout << "Edge equations:\n"; 
      get_edge_equations(T->M(), eqns); 
      verify_equations_approx(eqns, shapes);

      cout << "Cusp equations:\n"; 
      get_filling_equations(T->M(), eqns); 
      verify_equations_approx(eqns, shapes);
    }
    break;

  case PrintHomology:
    {
      AbelianGroup* hg = homology(T->M());
      if (!hg) {
	cout << "Unable to compute the homology group of this manifold.\n"; 
	break;
      }

      compress_abelian_group(hg);
      print_abelian_group(hg); 
      cout << "\n";
      free_abelian_group(hg); 
    }
    break; 

  case PrintGroup:
    {
      if (!get_group()) break; 
      print(G, true); 
    }
    break;

  case PrintSnappeaSymmetries:
    {
      SymmetryGroup *symm_gp=NULL;
      SymmetryGroup *symmetry_group_of_link=NULL;
      Triangulation *symmetric_triangulation=NULL;
      Boolean is_full_group;
      FuncResult res = compute_symmetry_group(T->M(), &symm_gp, &symmetry_group_of_link, &symmetric_triangulation, &is_full_group); 
      if (res != func_OK) {
	cout << "Symmetry group computation failed.\n";
	break; 
      }
      if (!is_full_group) {
	cout << "Only part of the symmetry group was found.\n";
      }
      print_symmetry_group(symm_gp); 
      if (get_num_cusps(T->M())==1) {
	cout << "Group of an invertible knot: "; 
	if (symmetry_group_invertible_knot(symm_gp)) 
	  cout << "yes" << endl; 
	else 
	  cout << "no" << endl;
      }
      cout << "Amphicheiral: ";
      if (symmetry_group_is_amphicheiral(symm_gp)) 
	cout << "yes" << endl; 
      else 
	cout << "no" << endl;
    }
    break;

  case SaveGroup:
    {
      if (!get_group()) break; 
      string fname = choose_save_file(".", name(1)+".grp"); 
      if (!fname.size()) break; 
      FILE* fp = fopen(fname.c_str(), "w");
      if (!fp) {
	printf("Unable to open %s for writing.\n", fname.c_str());
	break; 
      }
      cout << "Writing file " << fname << endl; 
      save_presentation(fp, G, num_cusps()==0); 
      fclose(fp); 
    }
    break; 

  case SaveGenerators:
    {
      if (!get_group()) break; 
      string fname = choose_save_file(".", name(1)+".gens");
      if (!fname.size()) break; 
      FILE* fp = fopen(fname.c_str(), "w");
      if (!fp) {
	printf("Unable to open %s for writing.\n", fname.c_str());
	break; 
      }
      cout << "Writing file " << fname << endl; 
      int i, n=fg_get_num_generators(G);
      fprintf(fp, "%% Generators\n%d\n\n", n); 

      for (i=0; i<n; i++) 
	snappea_print_o31(fg_gen(G,i), fp); 
      
      fclose(fp); 
    }
    break;

  case SnappeaIsometry:
    j = choose_manifold("Find isometry with "); 
    if (!get_TT(j)) {
      cout << "There is no triangulation " << (j+1) << endl; 
      break; 
    }
    check_if_isometric(T->M(), get_TT(j)->M(), 1); 
    break;

  case PrintCoreLength:
    {
      Complex len;
      int n_cusps = get_num_cusps(T->M());
      int s_index;
      int i; 
      Boolean is_complete; 
      for (i=0; i<n_cusps; i++) {
	get_cusp_info(T->M(), i, 0, &is_complete, 0,0,0,0,0,0,0,0);
	if (is_complete) continue; 
	core_geodesic(T->M(), i, &s_index, &len, NULL);
	cout << "Core " << i << ": " << len << endl;
      }
    }
    break;
    
  case PrintHolonomies:
    print_holonomies(T->M()); 
    break;
    
  case PrintVolume:
    {
      int precision; 
      int old_cout_prec = cout.precision();
      double vol = volume(T->M(), &precision);
      cout.precision(precision); 
      cout << "Volume is: " << vol << "\n";
      cout.precision(old_cout_prec); 
      cout << "Estimated precision: " << precision << "\n"; 
    }
    break;
  case PrintChernSimons:
    {
      int precision, known;
      double cs = T->chern_simons(&precision,&known);
      if (!known) {
	printf("Chern-Simons value is not known (CS fudge not available).\n");
	break;
      }
      printf("Chern-Simons (mod 1/2): %.16g\n", cs);
      printf("Estimated precision: %d digits\n", precision); 
    }
    break; 
  case BarycentricSubdivision:
    T->test_bc_subdivision(); 
    break;
    
#if 0
  case PrintData:
    T->print(); 
    break;
#endif
  case StandardizeCurves:
    j = choose_manifold("Isometric manifold number? "); 
    if (!get_TT(j)) {
      cout << "There is no triangulation " << (j+1) << endl; 
      break; 
    }
    if (!T->standardize_curves(get_TT(j)->M())) 
      cout << "Curves not changed.\n"; 
    break;
    
  case ChangeCurves:
    {
      int n = get_num_cusps(T->M()); 
      char pr[100]; 
      sprintf(pr, "%d basis change matrices in format [[a,b;c,d],...]\n? ",n); 
      if (!get_input(s,pr,-1)) break; 
      MatrixInt22* change_matrices = new MatrixInt22 [n]; 
      if (!get_basis_change_matrices(n, s.c_str(), change_matrices)) {
	printf("Invalid input, curves not changed.\n"); 
	delete [] change_matrices;
	break;
      }
      T->change_peripheral_curves(change_matrices);
      delete [] change_matrices;
    }
    break; 

  case SplitEdge:
    {
      if (!get_input(s,"Edge to split: edge face1 face2 ?",3)) break; 
      int e, f[2];
      if (sscanf(s.c_str(), "%d %d %d", &e, &f[0], &f[1])!=3) break; 
      // show_holonomy_contribs = 1; 
      // compute_holonomies(T->M()); 
      // print_holonomies(T->M());
      T->split_edge(e, f); 
      // compute_holonomies(T->M()); 
      // print_holonomies(T->M());
      // show_holonomy_contribs = 0; 
    }
    break;

  case CheckCurves:
    if (check_peripheral_curves(T->M()))
      cout << "all OK\n";
    else 
      cout << "peripheral curves don't look right\n";
    break; 
    
  case PrintCuspShapes:
    {
      int i, n = get_num_cusps(T->M());
      Boolean is_complete;
      Complex shape; 
      for (i=0; i<n; i++) {
	get_cusp_info(T->M(), i, NULL, &is_complete, NULL, NULL, NULL, &shape, 
		      NULL, NULL, NULL, NULL); 
	if (is_complete) cout << "cusp " << i << ": " << shape << endl; 
      }
    }
    break;
  case CuspVolumes:
    {
      if (num_filled_cusps(T->M())) {
	cout << "Manifold may not have any filled cusps\n"; 
	break; 
      }
      int i, c = get_num_cusps(T->M()); 
      n_vector S = get_cusp_sizes(c); 
      if (!S.dim) break; 
      CuspNeighborhoods *cn = initialize_cusp_neighborhoods(T->M());

      // set the initial displacements. 
      for (i=0; i<c; i++) { 
	set_cusp_neighborhood_displacement(cn, i, log(S[i])); 
      }
      double cvol, total_cvol=0.;
      printf("Bumping order: "); 
      bump_cusps(cn, 1); 
      for (i = 0; i<c; i++) {
	cvol = get_cusp_neighborhood_cusp_volume(cn, i);
	total_cvol += cvol; 
	printf("Cusp %d volume: %.10g\n", i, cvol); 
      }
      printf("Cusp density:  %.10g\n", total_cvol/T->volume()); 
      free_cusp_neighborhoods(cn); 
    }
    break;

  case CuspDensity:
    {
      CuspNeighborhoods *cn = initialize_cusp_neighborhoods(T->M());
      int i, n = get_num_cusps(T->M()), ci = 0; 
      Boolean complete; 
      for (i=0; i<n; i++) {
	get_cusp_info(T->M(),i,0,&complete,0,0,0,0,0,0,0,0);
	if (!complete) continue; 
	set_cusp_neighborhood_tie(cn,ci++,TRUE); 
      }
      double cusp_volume = 0.0; 
      if (ci) {
	set_cusp_neighborhood_displacement(cn, 0, 1e5); // Maximize cusp nbhd. 
	for (i = 0; i<ci; i++) {
	  cusp_volume += get_cusp_neighborhood_cusp_volume(cn, i); 
	}
	printf("Cusp density: %.16g\n", cusp_volume/T->volume()); 
      } else { 
	printf("Manifold has no unfilled cusps.\n"); 
      }
      free_cusp_neighborhoods(cn); 
    }
    break; 

  case PrintOrthodistances:
    {

      GroupPresentation* ug = T->group(0);
      int n = get_num_tetrahedra(T->M()); 
      vector<Complex> orth(n, Infinity);
      get_edge_orthodistances(T->M(), ug, orth); 

      cout << "Cosh orthodistances:\n" << PSeq(orth) << endl; 

      int i;
      Complex od; 
      for (i=0; i<n; i++) {
	od = complex_acosh(orth[i]); 
	if (od.real < 0.) od = -od; 
	orth[i] = od; 
      }
      cout << "Orthodistances: \n" << PSeq(orth) << endl; 

      print_orthoangles(T->M(), ug);
      free_group_presentation(ug); 
    }
    break;

  case PrintEdgeK:
    {
      GroupPresentation* ug = T->group(0);
      int n = get_num_tetrahedra(T->M()); 
      vector<Complex> orth(n, Infinity);
      get_edge_orthodistances(T->M(), ug, orth); 
      free_group_presentation(ug); 

      cout << "Edge K's:\n";
      Complex od; 
      double k; 
      int i; 
      cout << "[";
      for (i=0; i<orth.size(); i++) {
	if (i!=0) cout << ',';
	od = complex_acosh(orth[i]); 
	k = sin(od.imag)/sinh(od.real);
	if (k<0) k = -k; 
	fwprint(cout, k, 0);
      }
      cout << "]\n";
    }
    break; 


  case -1:
    cout << "Input not recognized.\n";
    break; 
    
  default:
    cout << "Command not available.\n";
  }
  
}

