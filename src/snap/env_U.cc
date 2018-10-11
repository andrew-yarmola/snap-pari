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

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <ctype.h>
#include <fstream>
#include <algorithm>
#include "snappea/unix_io.h"
#include "snappea/kernel.h"
#include "kernel_extras.hh"
#include "orthoangles.hh"
#include "helpers.hh"
#include "drill_Dirichlet.hh"
#include "warn.hh"
#include "i_triangulation.hh"
#include "env_U.hh"

using std::max;
using std::setw;
using std::swap;
using std::sscanf;

string env_U::gv_picture_opts("-fbp-ftw-Pbo");
string env_U::picture_options("-fbp-b");

enum {

  SetIdealVertexCutoff =500,
  SetEqsEpsilon,
  SetFaceEpsilon, 
  SetTubeEpsilon,
  SetGvOptions,
  SetPictureOptions,

  SetFaceTranspositions,

  ResetHolonomies, 
  ClearTube,

  // Group needed (for Dirichlet domain).

  ComputeTube,
  AddTubeFaces,
  DrillGeodesic,

  DrillDirichlet,

  // Triangulation needed. 

  ComputeCoreTube,
  PrintCoreWords,
  PrintCoreClippings,
  PrintCoreK,
  LocalCoreTube,
  CheckCoreTube,
  AddCoreWord,
  DeleteCoreWord,

  CurveInfo,
  PeripheralCrv,
  AllPeripheral,

  // Tube domain needed.

  PrintTubeInfo, 
  PrintFaces,
  PrintNaturalFaces,
  PrintConnectedFaces,
  PrintFaceMatrices, 
  CheckFacePairings,
  CheckEdgeMatching,
  DrillTube, 
  PrintDrillingData, 
  PrintCorners,
  SavePicture,
  SaveFacePicture,
  SaveNaturalPicture,
  SaveCurvePicture,
  GvSaveTube
};


void env_U::setup_menu()
{
  env_D::setup_menu(); 

  mm.add_item("set ideal_vertex_cutoff", SetIdealVertexCutoff); 
  mm.add_item("set edge epsilon", SetEqsEpsilon); 
  mm.add_item("set face epsilon", SetFaceEpsilon); 
  mm.add_item("set tube epsilon", SetTubeEpsilon); 
  mm.add_item("set picture", SetPictureOptions); 
  mm.add_item("set gv_options", SetGvOptions); 
  mm.add_item("reset_holonomies", ResetHolonomies); 
  mm.add_item("tube info", PrintTubeInfo); 
  mm.add_item("print core ortholines", PrintCoreWords);
  mm.add_item("add core word", AddCoreWord);
  mm.add_item("delete core word", DeleteCoreWord);
  mm.add_item("print core k", PrintCoreK);
  mm.add_item("curve", CurveInfo); 
  mm.add_item("compute tube", ComputeTube); 
  mm.add_item("drill geodesic", DrillGeodesic);
  mm.add_item("local core tube", LocalCoreTube); 
  mm.add_item("compute core tube", ComputeCoreTube); 
  mm.add_item("check core_tube", CheckCoreTube); 
  // mm.add_item("print core clippings", PrintCoreClippings); 
  mm.add_item("add faces", AddTubeFaces); 
  mm.add_item("clear tube", ClearTube); 
  mm.add_item("drill dirichlet", DrillDirichlet); 

  mm.add_item("peripheral_curve", PeripheralCrv); 
  mm.add_item("all_peripheral", AllPeripheral); 

  mm.add_item("print faces", PrintFaces);
  mm.add_item("print natural", PrintNaturalFaces);
  mm.add_item("print connected_faces", PrintConnectedFaces);
  // mm.add_item("print matrices", PrintFaceMatrices);
  mm.add_item("check face_pairings", CheckFacePairings); 
  mm.add_item("check edge_matching", CheckEdgeMatching); 
  mm.add_item("drill tube", DrillTube); 
  mm.add_item("print drilling_data", PrintDrillingData); 
  mm.add_item("print corners", PrintCorners); 
  mm.add_item("save picture", SavePicture);
  mm.add_item("save face_picture", SaveFacePicture);
  mm.add_item("save natural_picture", SaveNaturalPicture);
  mm.add_item("save curve_picture", SaveCurvePicture);
  mm.add_item("gv_save tube", GvSaveTube);
}

void env_U::validate_event(int& what)
{
  env_D::validate_event(what); 

  if (what < SetIdealVertexCutoff || what >= 1000) return; 

  if (what >= PrintTubeInfo && !tube_is_valid()) { 
    cout << "Please compute a tube first.\n"; 
    what = Nothing; 
  }

  if (what >= ComputeCoreTube && what < PrintTubeInfo && !T) { 
    cout << "This function requires a triangulation.\n";
    what = Nothing; 
  }

  if (what >= ComputeTube && !(T||G)) {
    cout << "This function requires a manifold.\n";
    what = Nothing; 
  }
}

void env_U::print_settings() const
{
  env_D::print_settings(); 

  // Want to use scientific notation for all epsilons. 
  ios_fmtflags old_fmt = cout.flags();
  cout.unsetf(ios::fixed);

  cout << "Tube edge epsilon:   " << eqs_interval::epsilon << endl;
  cout << "Tube face epsilon:   " << tube_face::epsilon << endl; 
  cout << "Tube epsilon:        " << tube_face::boundary_epsilon << endl; 
  cout << "Ideal vertex cutoff: " << tube::ideal_cutoff << endl; 

  cout.flags(old_fmt); 

  cout << "2D picture options: " << picture_options << endl; 
  cout << "Geomview picture options: " << gv_picture_opts << endl; 
}

void env_U::T_shape_changed()
{
  if (core_ort.size()) compute_core_ortholines(); 
  clear_tube();
  env_D::T_shape_changed(); 
}

void env_U::T_changed()
{
  if (tube_is_valid()) clear_tube(); 
  core_ort.clear();
  env_D::T_changed(); 
}

void env_U::clear_tube()
{ 
  U.clear(); 

  tube_tile_rad = 0.0; 
  tube_ort.clear(); 
  used_ort.clear(); 
  g_num = -1; 
}

void env_U::print_core_words() const
{
  if (!core_ort.size()) {
    cout << "You must call compute_core_ortholines() before print_core_words().\n";
    return; 
  }
  
  // First print the meridian and longitude words.
  cout << "Meridian  :" << holonomy[0] << ' ' << fg_meridian(G, 0) << endl; 
  cout << "Longitude :" << holonomy[1] << ' ' << fg_longitude(G, 0) << endl; 

  int i, n = core_ort.size(); 
  for (i=0; i<n; i++) {
    cout << core_ort[i] << endl; 
  }
}

void env_U::add_core_word(FGWord const& w)
{
  Ortholine O(One, CGRay(Zero,-2), CGRay(Zero, -2)); 
  O.word = w; 
  core_ort.push_back(O);
}

void env_U::delete_core_word(FGWord const& w)
{
  int i, n = core_ort.size(); 
  for (i=0; i<n; i++) {
    if (w==core_ort[i].word) {
      core_ort.erase(core_ort.begin()+i);
      break; 
    }
  }
}

void env_U::delete_core_word(int n)
{
  if (n < 0 || n >= core_ort.size()) return; 
  core_ort.erase(core_ort.begin()+n);
}

double env_U::core_tube_radius() const
{
  if (!core_ort.size()) {
    cout << "You must call compute_core_ortholines() before core_tube_radius().\n";
    return 0.0; 
  }
  double rad = core_ort[0].distance().real/2.0, r;
  int i; 
  for (i=1; i<core_ort.size(); i++) {
    r = core_ort[i].distance().real/2.0;
    if (r < rad) rad = r; 
  }
  return rad;
}


// We initialize the tube if necessary for geodesic(gn). 
// Then we compute it using all ortholines available from the 
// current tiling. If prev_tile_rad is nonzero we assume that 
// the tube has already had all faces computed which arise from
// tiling out to prev_tile_rad (and was presumably incomplete) 
// so we only use ortholines coming from new layers of the
// tiling. 

int env_U::compute_tube(int gn, bool report)
{
  if (!get_Dirichlet_domain()) return 3; 

  if (gn < 0 && g_num < 0) return 3;

  if (gn >=0 && gn != g_num) {
    clear_tube(); 
    g_num = gn; 

    if (g_num >= num_geodesics())
      if (!compute_geodesics(g_num+1)) return 3; 

    U.set_holonomies(TwoPiI, geodesic(g_num).length);
  }

  list<interval> ivls;
  get_crossing_lifts(domain, GSpec(g_num, &geodesic(g_num)), ivls, true);

  double rad = 3.;
  int res; 
  while (rad < 5.51) {
    new_ortholines(domain, ivls, geodesics, rad, tube_ort);
    res = U.compute_tube(tube_ort, used_ort, rad, report); 
    if (res != 1) break; 
    rad += 0.5; 
  }

  if (res==0) cout << "Tube has " << U.num_faces() << " faces.\n"; 
  return res;
}

static int standardize_drilled(Triangulation* ref, Triangulation* new_m, bool redo=true)
{
  int nc = get_num_cusps(new_m); 
  if (nc==1) return 1; 
  
  // If original manifold had cusps we want cusp bases to agree with 
  // those of the new manifold. 
  
  // Give up if original manifold has filled cusps.. just too 
  // confusing at the moment!
  
  if (get_num_cusps(ref)!=(nc-1)) {
    cout << "Sorry, can't restore peripheral curves";
    cout << " when starting manifold has filled cusps.\n";
    return 0; 
  }
  
  // First do (1,0) dehn filling on the newly drilled cusp. 
  
  set_cusp_info(new_m, 0, FALSE, 1., 0.);
  do_Dehn_filling(new_m);

  // Seek an isometry. 
  
  IsometryList *isometry_list = 0;
  Boolean isometric; 
  compute_isometries(ref,new_m,&isometric,&isometry_list,NULL); 
  if (!isometric) {
    printf("Oops - filled drilled manifold is not isometric!\n");
    return 0; 
  }
  if (!isometry_list) {
    printf("no isometry list found for filled drilled manifold!\n");
    return 0; 
  }

  int i;
  MatrixInt22 *pcb = new MatrixInt22 [nc];
  for (i=0; i<nc; i++) {
    pcb[i][0][0] = 1; 
    pcb[i][0][1] = 0; 
    pcb[i][1][0] = 0; 
    pcb[i][1][1] = 1; 
  }

  // Look for an isometry in which change mats have det 1. 

  int d, img_cusp; 
  for (i=1; i<nc; i++) {
    isometry_list_cusp_action(isometry_list, 0, i-1, &img_cusp, pcb[i]);
    swap(pcb[i][0][1],pcb[i][1][0]);
    d = pcb[i][0][0]*pcb[i][1][1] - pcb[i][1][0]*pcb[i][0][1];
    if (d==-1) break; 
  }

  // If d==-1 we need to reorient new_m first. 
  if (d==-1) {
    delete [] pcb; 
    if (redo) {
      reorient(new_m); 
      int res = standardize_drilled(ref, new_m, false); 
      if (!res) reorient(new_m);
      return res; 
    }
    return 0; 
  }

  // Change curves on new_m to agree with those on m. 

  if (change_peripheral_curves(new_m, pcb)!=func_OK) {
    delete [] pcb;
    return 0; 
  }
  delete [] pcb;
  return 1; 
}


void env_U::drill_tube(int manifold_num)
{
  Triangulation* new_m = U.drill(); 
  if (!new_m) return; 

  // Fill in the new triangulation name. 
  char newname[250]; 
  sprintf(newname, "%s-[%d]", name().c_str(), geodesic_number());
  set_triangulation_name(new_m, newname); 

  if (T) standardize_drilled(T->M(), new_m);

  get_T(manifold_num)->set_T(new_manifold(new_m));
  
  double vol = volume(new_m, NULL);
  cout << "Volume of complement: " << vol << endl; 
}

// This should be called any time tube_face::boundary_epsilon or
// tube_face::epsilon are changed.

bool env_U::recheck_face_pairings()
{
  double urad; 
  bool match = U.check_face_pairings(urad);
  if (match) {
    U.compute_connected_faces();
  }
  return match;
}


static Complex signed_complex_length(MoebiusTransformation const& mt, Complex const& z)
{
  if (!same_point(z, mt*z)) {
    cout << "signed_complex_length called with incorrect fixed point.\n";
    return Zero; 
  }

  if (complex_big(z))
    return -complex_log(mt.matrix[0][0]/mt.matrix[1][1], 0.0);

  // In usual notation w = cz + d. 
  Complex w = mt.matrix[1][0] * z + mt.matrix[1][1];
  return -complex_log(w*w, 0.0);
}

bool env_U::compute_core_ortholines()
{
  if (simplify) {
    cout << "Switching to unsimplified fundamental group\n";
    simplify = false; 
    clear_group(); 
  }

  if (get_num_cusps(T->M())!=1) return false; 

  if (!get_group()) return false; // eg. if structure is degenerate

  FGWord mw = fg_meridian(G, 0);
  FGWord lw = fg_longitude(G, 0);

  MoebiusTransformation m = word_to_Moebius(G, mw); 
  MoebiusTransformation l = word_to_Moebius(G, lw); 
  Complex basepoint = cusp_vertex(T->M(), 0); 

  Boolean complete; 
  get_cusp_info(T->M(), 0, 0, &complete, 0, 0, 0, 0, 0, 0, 0, 0);

  int i, n; 
  if (core_ort.size()) {
    if (complete) {
      holonomy[0] = Zero; 
      holonomy[1] = Zero; 
      n = core_ort.size();
      for (i=0; i<n; i++) 
	core_ort[i].set_to(Infinity, Zero, Zero); 
    } else {
      holonomy[0] = nearby_complex(signed_complex_length(m, basepoint), 
				   holonomy[0]); 
      holonomy[1] = nearby_complex(signed_complex_length(l, basepoint), 
				   holonomy[1]);
      
      vector<Ortholine> old_core_ort(core_ort); 
      words_to_ortholines(G, T->M(), core_ort, true); 
    }
  } else {
    holonomy[0] = signed_complex_length(m, basepoint); 
    holonomy[1] = signed_complex_length(l, basepoint);

    vector<FGWord> core_ort_words; 
    get_end_pairing_words(T->M(), core_ort_words); 
    n = core_ort_words.size();
    core_ort.resize(n); 
    for (i=0; i<n; i++) 
      core_ort[i].word = fg_word_from_original(G, core_ort_words[i]); 

    words_to_ortholines(G, T->M(), core_ort, false); 
  }

  return true;
}

void env_U::print_tube_ortholines() const
{
  int i, n = U.the_ortholines().size();
  FGWord wd;
  Ortholine o; 
  for (i=0; i<n; i++) {
    o = U.the_ortholines()[i];
    cout << "Faces " << (2*i) << ',' << (2*i+1);

    if (!core_ort.size()) { 
      cout << endl; 
    } else if (ortholine_to_word(o,wd)) {
      cout << ": " << wd << endl; 
    } else {
      cout << "\nProblem computing ortholine word.\n";
    }
    cout << o << endl;
  }
}

bool env_U::ortholine_to_word(Ortholine const& o, FGWord& wd) const
{
  if (!G) return false; 

  R_matrix<2> H;
  H[0][0] = holonomy[0].real;
  H[1][0] = holonomy[0].imag;
  H[0][1] = holonomy[1].real;
  H[1][1] = holonomy[1].imag;
  H.invert(); 

  // H now maps the fundamental torus to the unit square. 

  R_vector<2> pos, p2;

  double m, l; 

  int j, k, e, f, mi, li;
  FGWord w, pw;

  for (e=0; e<2; e++) {
    pos[0] = o.position(e).real;
    pos[1] = o.position(e).imag;
    pos = H * pos; 

    for (f=0; f<2; f++) {
      for (j=0; j<core_ort.size(); j++) {
	p2[0] = core_ort[j].position(f).real;
	p2[1] = core_ort[j].position(f).imag;
	p2 = H * p2; 
	
	p2 -= pos; 
	// p2 should now be integral m,l if ortholine positions 
	// were equal on the torus. 
	
	m = floor(p2[0] + .5);
	l = floor(p2[1] + .5); 
	if (fabs(p2[0]-m) + fabs(p2[1]-l) < 1e-5) break;
      }
      if (j<core_ort.size()) break; 
    }
    if (f == 2) return false; 

    // core_ort[j] is now a match for o.

    mi = int(m);
    li = int(l); 

    w = FGWord(); 
    pw = (mi >= 0) ? inverse(fg_meridian(G, 0)) : fg_meridian(G, 0); 
    if (mi < 0) mi = -mi;
    for (k=0; k<mi; k++) w *= pw;
    pw = (li >= 0) ? inverse(fg_longitude(G, 0)) : fg_longitude(G, 0); 
    if (li < 0) li = -li;
    for (k=0; k<li; k++) w *= pw;

    // what's this, this looks like rubbish!
    if (e==0) {
      wd = w * core_ort[j].word;
    } else {
      wd *= inverse(w); 
    }
  }
  return true;
}


void env_U::compute_core_tube(bool report) 
{
  if (geodesic_number()!=-1) clear_tube(); 
  if (core_ort.size()==0 && !compute_core_ortholines()) {
    cout << "Problem computing core ortholines.\n";
    return; 
  }

  U.set_holonomies(holonomy[0], holonomy[1]);
  int res = U.compute_tube(core_ort, report);
  if (res==0) cout << "Tube has " << U.num_faces() << " faces.\n"; 
  g_num = -2; 
}

int env_U::verify_core_tube()
{
  static int call=1; 
  // if (geodesic_number()!=-1) clear_tube(); 
  if (core_ort.size()==0 && !compute_core_ortholines()) {
    cerr << "Problem computing core ortholines.\n";
    return 2; 
  }

  tube* t = new tube(); 

  t->set_holonomies(holonomy[0], holonomy[1]);
  // g_num = -2; 
  // if (call==6) cout << endl; 
  int res = t->compute_tube(core_ort, false);
  delete t; 
  call++; 
  return(res);
}


void env_U::print_clip_list() const
{
  tube S;
  if (!core_ort.size()) {
    cout << "Please call compute_core_ortholines() first.\n";
    return;
  }
  S.set_holonomies(holonomy[0], holonomy[1]);
  S.print_clip_list(core_ort);
}

void env_U::save_picture(string const& file, string const& options) const
{
  Complex mer, lon; 
  U.get_holonomies(mer, lon); 

  double diam = max(complex_modulus(mer+lon), complex_modulus(mer-lon));

  // Open a picture file and set the scale so that the diameter of the 
  // fundamental parallelogram is 4 inches.
  ps_picfile pic = ps_picfile(file.c_str());
  pic.set_scale(72.0 * 4.0/diam); 
  double inch = diam/4.0;

  U.save_picture(options, pic); 

  char buf[255]; 
  sprintf(buf, "manifold: %s, geodesic: %d",
	  name().c_str(), geodesic_number());
  pic.print_text(inch * Complex(-3.25,-4.5), string(buf));
  
  pic.close();
  return; 
}

void env_U::process_event(int what)
{
  string s; 

  switch (what) {

  case SetEqsEpsilon:
    {
      if (!get_input(s,"tube edge epsilon ? ")) break; 
      sscanf(s.c_str(), "%lf", &eqs_interval::epsilon); 
    }
    break;
  case SetFaceEpsilon:
    {
      double new_face_epsilon;
      if (!get_input(s,"tube face epsilon ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &new_face_epsilon) != 1) break; 
      tube_face::epsilon = new_face_epsilon; 
      if (tube_is_valid()) recheck_face_pairings(); 
    }
    break; 
  case SetTubeEpsilon:
    {
      double new_epsilon;
      if (!get_input(s,"tube epsilon ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &new_epsilon) != 1) break; 
      tube_face::boundary_epsilon = new_epsilon; 
      if (tube_is_valid()) recheck_face_pairings(); 
    }
    break; 
  case SetIdealVertexCutoff:
    if (!get_input(s,"Ideal vertex cutoff (default is 5.) ? ")) break; 
    sscanf(s.c_str(), "%lf", &tube::ideal_cutoff);
    break;
  case SetFaceTranspositions:
    {
      transpose_no_faces(); 
      vector<int> to_transpose;
      get_input(s, "faces to transpose (return for none) ? ", -1);
      if (s.length() > 0) to_transpose = read_numbers(s.c_str()); 
      int i; 
      for (i=0; i<to_transpose.size(); i++)
	transpose_face(to_transpose[i]); 
    }
    break;
    
  case SetGvOptions:
    if (get_input(s, "Geomview picture options required? "))
      gv_picture_opts = s;
    break;
    
  case SetPictureOptions:
    if (get_input(s, "2D picture options required? "))
      picture_options = s;
    break;
    
  case PrintCoreK:
    {
      Ortholine S;
      int i, n = core_ort.size(); 
      for (i=0; i<n; i++) {
	S = core_ort[i];
	cout << S.distance() << ' '; 
	fwprint(cout, sin(S.distance().imag)/sinh(S.distance().real), 0);
	cout << ' ' << S.word << endl; 
      }
    }
    break;
    
  case DrillDirichlet:
    {
      int face_num;
      if (!get_input(s,"faces to pair ? ")) {
	print_dirichlet_info(domain); 
	break; 
      }
      if (sscanf(s.c_str(), "%d", &face_num) != 1) break; 
      
      // Find where user wants to save resulting manifold. 
      int j = choose_save_manifold(); 
      if (j<0) break; 
      
      triangulation_builder TB;
      get_dirichlet_info(domain, TB); 
      TB.glue_face(face_num); 
      TB.remove_redundant_edges(0);
      TB.remove_redundant_vertices(); 
      TB.compute_cell_complex(); 
      
      TriangulationData TD;
      if (!TB.get_triangulation_data(TD)) {
	cout << "Problem drilling this curve from Dirichlet domain.\n";
	break; 
      }
      
      // Fill in the new triangulation name. 
      char newbit[10]; 
      sprintf(newbit, "-D%d", face_num);
      string nn = name() + newbit;
      TD.name = new char [nn.length()+1];
      strcpy(TD.name, nn.c_str()); 
      
      // Store the new manifold. 
      Triangulation *new_m; 
      data_to_triangulation(&TD, &new_m);
      if (new_m) {
	get_T(j)->set_T(new_manifold(new_m));
      }
      
    }
    break; 
    
  case ResetHolonomies:
    { 
      // Analytic continuation of holonomies may not
      // always be successful. In this case the user needs
      // to do a large Dehn surgery and then set imaginary 
      // parts of the resulting holonomies close to zero. 
      int i; 
      for (i=0; i<2; i++) 
	holonomy[i] = nearby_complex(holonomy[i], Zero); 
    }
    break; 

  case CurveInfo:
    {
      if (!get_input(s, "Which curve? ", 2)) break;
      
      double cm, cl;
      if (sscanf(s.c_str(), "%lf %lf", &cm, &cl) != 2) break;
      
      int cusp = 0; 
      if (get_num_cusps(T->M()) > 1) {
	if (!get_input(s, "Which cusp? ")) break;
	if (sscanf(s.c_str(), "%d", &cusp) != 1) break;
	if (cusp < 0 || cusp >= get_num_cusps(T->M())) {
	  cout << "Cusp must be in the range 0 to " << get_num_cusps(T->M()) << "\n"; 
	  break;
	}
      }
      
      Complex mh, lh;
      get_holonomy(T->M(), cusp, &mh, &lh, 0, 0);
      Complex len = cm * mh + cl * lh;
      
      cout << "Complex length: " << len << "\n"; 
      
      double rad = core_tube_radius(); 
      cout << "Tube radius: " << rad << "\n";
      
      double s = sinh(rad);
      double c = cosh(rad);
      double L = sqrt(s * s * len.imag * len.imag + c * c * len.real * len.real); 
      double K = sqrt(c * c * len.imag * len.imag + s * s * len.real * len.real); 
      
      cout << "Length: " << L << "\n";
      cout << "Curvature: " << K << "\n"; 
    }
    break;
    
  case AllPeripheral:
    {
      double limit;
      if (!get_input(s, "Length limit? ")) break; 
      if (sscanf(s.c_str(), "%lf", &limit)!=1) break; 

      // get cusp number if more than one cusp
      int cusp = 0, n = get_num_cusps(T->M()); 
      if (n > 1) {
	if (!get_input(s, "Which cusp? ")) break;
	if (sscanf(s.c_str(), "%d", &cusp) != 1) break;
	if (cusp < 0 || cusp >= get_num_cusps(T->M())) {
	  cout << "Cusp must be in the range 0 to " << 
	    get_num_cusps(T->M()) << "\n"; 
	  break;
	}
      }

      // get shape
      Complex shape; 
      Boolean complete; 
      get_cusp_info(T->M(),cusp,0,&complete,0,0,0,&shape,0,0,0,0);

      cout << "surgery npc-len   pc-len           core-len  ac-core   volume       dV     L/~L   dV/~dV  dV/c~dV   tube-r\n";

      int me,lo;
      double maxab = sqrt(shape.imag) * limit; 
      double mmid, memax;
      double ln; 
      for (lo=0; lo<maxab; lo++) {
	mmid = -lo*shape.real; 
	me = (int)floor(mmid - maxab -.1);
	if (lo==0) me = 1; // avoids duplications. 
	memax = mmid + maxab + .1; 
	for (; me < memax; me++) {
	  if (me==0 && lo==0) continue; 
	  ln = sqrt(complex_modulus_squared(lo*shape + Complex(me,0.))/
		   shape.imag); 
	  if (ln > limit) continue; 
	  print_peripheral_curve_info(cusp, me, lo); 
	}
      }

    }
    break;

  case PeripheralCrv:
    {
      if (!get_input(s, "Which curve? ", 2)) break;
      
      double cm, cl;
      if (sscanf(s.c_str(), "%lf %lf", &cm, &cl) != 2) break;
      
      // get cusp number if more than one cusp
      int cusp = 0, n = get_num_cusps(T->M()); 
      if (n > 1) {
	if (!get_input(s, "Which cusp? ")) break;
	if (sscanf(s.c_str(), "%d", &cusp) != 1) break;
	if (cusp < 0 || cusp >= get_num_cusps(T->M())) {
	  cout << "Cusp must be in the range 0 to " << 
	    get_num_cusps(T->M()) << "\n"; 
	  break;
	}
      }

      print_peripheral_curve_info(cusp, cm, cl); 

    }
    break;

  case AddTubeFaces:
    compute_tube(-1,true);
    break;
    
  case ComputeTube:
    {
      int gnum; 
      if (!(get_input(s,"which geodesic ? ") &&
	    (sscanf(s.c_str(), "%d", &gnum) == 1))) break;
      compute_tube(gnum,false);
    }
    break;
    
  case ComputeCoreTube:
    compute_core_tube(false); 
    break;
  case LocalCoreTube:
    {
      U.clear();
      if (!get_group()) break; 
      int res = U.compute_tube(T->M(), G, true);
      if (res == 2) {
	cout << "tube computation failed\n";
	break; 
      } else if (res == 1) {
	cout << "non-tubal triangulation\n";
      } else if (res==0) { 
	cout << "Tube has " << U.num_faces() << " faces.\n"; 
      }
      g_num = -2; 
    }
    break; 
    
  case CheckCoreTube:
    {
      if (!get_group()) break; 
      int res = triangulation_is_tubal(T->M(), G);
      if (res==0) {
	cout << "triangulation is tubal\n";
      } else if (res==1) {
	cout << "non-tubal triangulation\n";
      } // error message already given when res==2
    }
    break; 

  case PrintCoreWords:
    if (!get_group()) break; 
    if (!core_ort.size()) compute_core_ortholines(); 
    print_core_words();
    break; 
    
  case AddCoreWord:
    if (!get_input(s,"word to add ? ")) break; 
    add_core_word(FGWord(s.c_str())); 
    break; 
    
  case DeleteCoreWord:
    if (!get_input(s,"word to delete ? ")) break; 
    delete_core_word(FGWord(s.c_str())); 
    break; 
    
  case PrintCoreClippings:
    if (!core_ort.size()) compute_core_ortholines(); 
    print_clip_list(); 
    break;
    
  case ClearTube:
    clear_tube(); 
    break;
    
    // EVERYTHING FROM HERE ON REQURES m TO HAVE A VALID TUBE
  case PrintFaces:
    the_tube().print_faces(false); 
    break;
    
  case PrintNaturalFaces:
    the_tube().print_faces(true); 
    break;
    
  case PrintConnectedFaces:
    the_tube().print_connected_faces(); 
    break;
    
  case PrintFaceMatrices:
    {
      int j; 
      if (!get_input(s,"Which face ? ")) break; 
      if (sscanf(s.c_str(), "%d", &j)!=1) break; 
      the_tube().print_face_matrices(j); 
    }
    break;
    
  case PrintTubeInfo:
    {
      Complex me, lo; 
      if (!tube_is_valid()) {
	cout << "No tube currently computed in this manifold.\n";
	break;
      }
      cout << "Geodesic: " << geodesic_number() << endl; 
      the_tube().get_holonomies(me, lo);
      cout << "Holonomies: " << me << ", " << lo << endl;
      if (geodesic_number()==-2) {
	print_tube_ortholines();
      } else {
	the_tube().print_ortholines(); 
      }
    }
    break;
    
  case CheckFacePairings:
    {
      double urad; 
      the_tube().check_face_pairings(urad, 1); 
    }
    break;

  case CheckEdgeMatching:
    {
      U.tidy_edges(true); 
    }
    break; 
    
  case PrintDrillingData:
    U.print_triangulation_data(); 
    break; 
    
  case PrintCorners:
    the_tube().print_face_corners(); 
    break; 
    
  case DrillGeodesic:
    {
      int gnum; 
      if (!(get_input(s,"which geodesic ? ") &&
	    (sscanf(s.c_str(), "%d", &gnum) == 1))) break;
      int mnum = choose_save_manifold(); 
      if (mnum<0) break; 
      
      if (compute_tube(gnum)!=0) { 
	cout << "tube computation failed!\n";
	break; 
      }
      drill_tube(mnum); 
    }
    break;
  case DrillTube:
    {
      int mnum = choose_save_manifold(); 
      if (mnum<0) break; 
      drill_tube(mnum);
    }
    break;
    
  case SavePicture:
    save_picture("tube.ps", picture_options); 
    cout << "wrote file tube.ps\n";
    break;
    
  case SaveFacePicture:
    {
      ps_picfile pic = ps_picfile("tube.ps"); 
      pic.set_scale(72.0);
      the_tube().save_facepairing_picture(pic);
      pic.close(); 
      cout << "wrote file tube.ps\n";
    }
    break;
    
  case SaveNaturalPicture:
    {
      ps_picfile pic = ps_picfile("tube.ps"); 
      pic.set_scale(72.0);
      the_tube().save_natural_face_picture(pic);
      pic.close(); 
      cout << "wrote file tube.ps\n";
    }
    break;
    
  case SaveCurvePicture:
    {
      if (!get_input(s,"meridian or longitude ? ")) break; 
      if (s[0]!='m' && s[0]!='M' && s[0]!='l' && s[0]!='L') break;
      
      Complex curve;
      if (s[0]=='m' || s[0]=='M')
	curve = TwoPiI;
      else 
	curve = geodesic(geodesic_number()).length;
      
      ps_picfile pic = ps_picfile("tube.ps"); 
      pic.set_scale(2.0 * 72.0/complex_modulus(curve)); // Make curve 2 inches. 
      the_tube().covered_curve_picture(pic, curve);
      pic.close(); 
      cout << "wrote file tube.ps\n";
    }
    break;
    
  case GvSaveTube:
    the_tube().save_gv_picture(gv_picture_opts); 
    break;
    
  default:
    env_D::process_event(what); 
  }
}

vector<Ortholine> env_U::get_core_ortholines() const 
{
  if (get_num_cusps(T->M())!=1) 
    return vector<Ortholine>(); 

  FGWord mw = fg_meridian(G, 0);
  FGWord lw = fg_longitude(G, 0);

  MoebiusTransformation m = word_to_Moebius(G, mw); 
  MoebiusTransformation l = word_to_Moebius(G, lw); 

  vector<FGWord> core_ort_words; 
  get_end_pairing_words(T->M(), core_ort_words); 
  
  int i, n = core_ort_words.size(); 
  vector<Ortholine> core_ort(n); 

  for (i=0; i<n; i++) 
    core_ort[i].word = fg_word_from_original(G, core_ort_words[i]); 

  words_to_ortholines(G, T->M(), core_ort, false); 
  return core_ort;
}

void env_U::print_peripheral_curve_info(int cusp, double cm, double cl)
{
  if (!get_group()) return;

  // get shape
  Complex shape; 
  Boolean complete; 
  get_cusp_info(T->M(),cusp,0,&complete,0,0,0,&shape,0,0,0,0);
  
  // quit if cusp is filled
  if (!complete) {
    // cout << "Cusp must be complete\n"; 
    return; 
  }
      
  // get normalized length
  // cout << "Cusp shape: " << shape << endl; 

  double ln = sqrt(complex_modulus_squared(cl * shape + Complex(cm,0.))/
		   shape.imag); 

  // get cusp area
  CuspNeighborhoods *cn = initialize_cusp_neighborhoods(T->M());
  int i, ci = 0, n=get_num_cusps(T->M()); 
  for (i=0; i<n; i++) {
    get_cusp_info(T->M(),i,0,&complete,0,0,0,0,0,0,0,0);
    if (!complete) continue; 
    set_cusp_neighborhood_tie(cn,ci++,TRUE); 
  }
  set_cusp_neighborhood_displacement(cn, 0, 1e5); // Maximize cusp nbhd. 
  double cusp_area = 2.0*get_cusp_neighborhood_cusp_volume(cn, cusp); 
  free_cusp_neighborhoods(cn); 
  
  // get absolute length
  double la = ln * sqrt(cusp_area); 
  
  // print out basic filling curve info.
  // cout << "Length: " << la << endl; 
  // cout << "Normalized length: " << ln << endl;

  cout << setw(3) << int(cm) << ' ' << setw(2) << int(cl) << ' ';
  fwprint(cout, ln, 8); 
  cout << ' '; 
  fwprint(cout, la, 8); 
  cout << ' '; 

  // get the complete volume
  double vol=volume(T->M(),0); 

  // get existing surgery coeffs
  double sm, sl; 
  vector<double> coeffs(2*n); 
  for (i=0; i<n; i++) {
    get_cusp_info(T->M(),i,0,0,&sm,&sl,0,0,0,0,0,0);
    coeffs[2*i] = sm; 
    coeffs[2*i+1] = sl; 
  }

  // do the surgery 
  int cone=int(6./ln); // do it gradually if ln < 6. 
  if (cone < 1) cone=1; 
  SolutionType sol; 
  for (; cone > 0; cone--) {
    coeffs[2*cusp] = cm*cone;
    coeffs[2*cusp+1] = cl*cone; 
    sol = T->do_Dehn_surgery(coeffs); 
    if (sol!=geometric_solution && sol!=nongeometric_solution) break; 
  }

  if (sol==geometric_solution || sol==nongeometric_solution) {

    // get the new volume
    double svol=T->volume(); 
    
    // get core geodesic length
    int sing_index;
    Complex core_length;
    core_geodesic(T->M(), cusp, &sing_index, &core_length, 0);
    double ac_core_len = core_length.real/sing_index; 

    // print core length
    // cout << "Core length: " << core_length << endl; 
    // cout << "A-C core length: " << ac_core_len << endl; 
    fwprint(cout, core_length, 8); cout << ' '; 
    fwprint(cout, ac_core_len, 8); cout << ' '; 


    // print out remaining info
    double vdiff = vol - svol; 
    double asy_vdiff = PI*PI/(ln*ln); 
    double asy_clen = TWO_PI/(ln*ln); 
    double core_vdiff_est = PI_OVER_2*ac_core_len; 

    // cout << "Volume: " << svol << endl; 
    // cout << "Volume change: " << vdiff << endl;
    // cout << "Len/Asymptotic-Len: " << (ac_core_len/asy_clen) << endl; 
    // cout << "VolCh/Asymptotic-VolCh: " << vdiff/asy_vdiff << endl; 
    // cout << "VolCh/Core-estimate: " << vdiff/core_vdiff_est << endl; 
    fwprint(cout, svol, 8); cout << ' '; 
    fwprint(cout, vdiff, 8); cout << ' '; 
    fwprint(cout, (ac_core_len/asy_clen), 8); cout << ' '; 
    fwprint(cout, (vdiff/asy_vdiff), 8); cout << ' '; 
    fwprint(cout, (vdiff/core_vdiff_est), 8); cout << ' '; 


    // Get tube radius assuming the triangulation contains a 
    // shortest ortholine. 
    double tube_rad, ctr; 
    vector<Ortholine> OL = get_core_ortholines(); 
    tube_rad = OL[0].distance().real/2.0;
    for (i=1; i<OL.size(); i++) {
      ctr = OL[i].distance().real/2.0;
      if (ctr < tube_rad) tube_rad = ctr; 
    }


    fwprint(cout, tube_rad, 8); 
    int vct = verify_core_tube(); 
    if (vct==1) cout << " ?";
    else if (vct==2) cout << " ??";
    cout << '\n';

  } else {
    cout << endl;
  }

  // restore original dehn fillings. 
  T->find_complete_hyperbolic_structure(); 
  if (n>1) {
    coeffs[2*cusp] = 0.0;
    coeffs[2*cusp+1] = 0.0; 
    T->do_Dehn_surgery(coeffs); 
  }
}

