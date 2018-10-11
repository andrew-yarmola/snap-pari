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
#include <fstream>
#include <cstdio>

using std::sscanf;

#include "env_A.hh"

#include "snappea/kernel_prototypes.h"
#include "snappea/unix_io.h"
#include "snappea/symmetry_group.h"
#include "kernel_extras.hh"
#include "closed_geodesics.hh"
#include "get_equations.hh"
#include "algebraic.hh"
#include "alg_group.hh"
#include "helpers.hh"
#include "record.hh"
#include "printable.hh"
#include "polish_rep.hh"
#include "pari_code.hh"
#include "dirichlet.hh"
#include "p_VS.hh"

#undef normalize

double env_A::acc_vertex_eps = 1e-20;
int env_A::cc_poly = 2; 
int env_A::print_fields_long = 0; 
int env_A::show_field_comp = 0; 
int env_A::hilb_print = 0; 
bool env_A::cs_diagnostic = false; 
text_menu env_A::fm; 
int env_A::sally_verbosity = 1; 
bool env_A::save_sally_results = false; 

using std::setw; 

const pari_printform pp_short = pari_printform(1,1,'g',16,0); 

static record fudge_record = string("M");


enum {
  // Print and change settings. 
  SetStackSize = 1000,
  SetPrintstyle, 
  SetHilbertPrintstyle,
  SetCSDiagnostic,

  SetAccurateVertexEpsilon,
  SetPrecision,
  SetDegree,
  SetLLLFudge,
  SetSallyVerbosity,
  SetSaveSally,
  CCPoly, 

  FieldPrintStyle,
  ShowFieldComp,

  // Misc. non manifold requiring stuff. 
  Status,
  TestMats,
  FactorMatrix, 
  EtaCuspAdj, 
  DedekindSum,

  // Group is optional, triangulation is optional.

  PrintManifold, 
  ComputeFields, 

  PrintFields,
  WriteField,

  VerifyEquations,

  // Group required.

  WordToAlgMatrix,
  EvaluateTrace,
  NormalizeGroup,
  NormalizeParabolics,
  SallysFunction, 
  PolishRep,

  AccurateDirichletDomain,

  ComputeInvariantTraceField, 
  ComputeTraceField, 
  ComputeGroupField, 

  PrintArithmetic,
  TestArithmeticTraces,

  ComputeEigenvalueFields, // Requires invt. trace field.
  ComputeUnits,

  // Anything after here requires a triangulation.

  Holonomies,
  PrintCSCorrection,

  PariEquations,
  PariFillingEquations,

  // PrintIsometries,
  // SetIsometry, 
  // ClearIsometries, 

  PrintNonPerZ2Homology,
  PrintCuspHomologyKernel, 

  PrintSignature, 
  CheckSignature, 

  Polish, 
  // Inaccurate, 

  // The following require accurate shape info. but no Dirichlet domain. 

  PrintFlattening,

  PrintAccurateShapes,
  AccurateLogShapes, 
  AccurateVolume,
  AccurateChernSimons, 
  ComplexVolume, 
  ComplexVolume2,

  BootstrapCS,

  PrintEta,
  PrintEtaInfo,
  BootEta,
  SetEta,
  SetEtaFudge,
  SaveEtaFudge, 
  CheckEta, 
  PrintEtaFudge,
  PrintEtaDifference, 
  PrintOldEtaFudge,
  PrintOldEtaDifference, 
  ExtendEtaFudge, 

  ComputeShapeField, 
  ComputeCuspFields, 

  // Stuff after here makes use of exact numbers
  PrintExactShapes,

  PrintInvariants,

  PrintBorel,
  PrintCuspInfo,
  PrintCuspCommensurability,

};

void env_A::setup_menu()
{
  env_U::setup_menu(); 

  mm.add_item("set stack_size", SetStackSize);
  mm.add_item("set precision", SetPrecision); 
  mm.add_item("set degree", SetDegree); 
  mm.add_item("set canonical_poly_mode", CCPoly); 
  mm.add_item("set lll_significant_digits", SetLLLFudge); 
  mm.add_item("set digits_printed", SetPrintstyle); 
  mm.add_item("set hilbert_printstyle", SetHilbertPrintstyle); 
  mm.add_item("set chern_simons_diagnostic", SetCSDiagnostic); 

  mm.add_item("set accurate_dirichlet epsilon", SetAccurateVertexEpsilon); 
  mm.add_item("set sally verbosity", SetSallyVerbosity); 
  mm.add_item("set sally save", SetSaveSally);
  mm.add_item("toggle field_print_style", FieldPrintStyle); 
  mm.add_item("toggle show_field_computation", ShowFieldComp); 

  mm.add_item("status", Status); 
  // mm.add_item("test mats", TestMats); 
  mm.add_item("factor", FactorMatrix); 
  mm.add_item("adjustment", EtaCuspAdj); 
  mm.add_item("dedekind", DedekindSum); 

  mm.add_item("print manifold", PrintManifold); 

  mm.add_item("print -pari gluing_equations", PariEquations);
  mm.add_item("print -pari filling_equations", PariFillingEquations);
  mm.add_item("print cusp homology_kernel", PrintCuspHomologyKernel);
  mm.add_item("print nonperipheral", PrintNonPerZ2Homology);

  mm.add_item("print fields", PrintFields);
  mm.add_item("gp_save_field", WriteField); 

  mm.add_item("print arithmetic", PrintArithmetic); 
  mm.add_item("print invariants", PrintInvariants); 
  mm.add_item("print borel_regulator", PrintBorel); 

  mm.add_item("print cusp invariant", PrintCuspInfo); 
  mm.add_item("print cusp commensurability", PrintCuspCommensurability); 
  mm.add_item("print eigenvalue_fields", ComputeEigenvalueFields); 
  mm.add_item("print units", ComputeUnits); 

  mm.add_item("print accurate shapes", PrintAccurateShapes); 
  mm.add_item("print accurate log_shapes", AccurateLogShapes); 
  mm.add_item("print exact shapes", PrintExactShapes); 

  mm.add_item("print cscorrection", PrintCSCorrection); 

  mm.add_item("test arithmetic", TestArithmeticTraces); 

  mm.add_item("compute invariant_trace_field", ComputeInvariantTraceField); 
  mm.add_item("compute trace_field", ComputeTraceField); 
  mm.add_item("compute shape_field", ComputeShapeField); 
  mm.add_item("compute group_field", ComputeGroupField); 
  mm.add_item("compute cusp_fields", ComputeCuspFields); 
  mm.add_item("compute fields", ComputeFields); 

  mm.add_item("print flattening", PrintFlattening); 

  mm.add_item("evaluate exact", WordToAlgMatrix); 
  mm.add_item("evaluate trace", EvaluateTrace); 

  mm.add_item("print volume", AccurateVolume); 
  mm.add_item("print complex_volume", ComplexVolume); 
  mm.add_item("print -2 complex_volume", ComplexVolume2); 
  mm.add_item("print chern_simons", AccurateChernSimons); 
  mm.add_item("accurate dirichlet", AccurateDirichletDomain); 

  mm.add_item("normalize group", NormalizeGroup); 
  mm.add_item("normalize parabolics", NormalizeParabolics); 
  mm.add_item("sally", SallysFunction); 


  mm.add_item("print signature", PrintSignature); 
  mm.add_item("check signature", CheckSignature); 
//  mm.add_item("print cusp maps", PrintCuspMaps); 

  mm.add_item("set eta invariant", SetEta); 
  mm.add_item("set eta fudge", SetEtaFudge); 
  mm.add_item("save eta_fudge", SaveEtaFudge); 
  mm.add_item("print eta invariant", PrintEta); 
  mm.add_item("print eta information", PrintEtaInfo); 
  mm.add_item("check eta", CheckEta); 
  mm.add_item("print eta fudge", PrintEtaFudge); 
  mm.add_item("print eta difference", PrintEtaDifference); 
  //  mm.add_item("print old eta fudge", PrintOldEtaFudge); 
  //  mm.add_item("print old eta difference", PrintOldEtaDifference); 
  //  mm.add_item("print eta rat_diff", PrintEtaRatDiff); 

  // mm.add_item("print isometries", PrintIsometries); 
  // mm.add_item("clear isometries", ClearIsometries); 

  mm.add_item("extend_fudge", ExtendEtaFudge); 

  mm.add_item("verify exact", VerifyEquations);

  mm.add_item("polish triangulation", Polish); 
  mm.add_item("polish group", PolishRep);
//  mm.add_item("inaccurate", Inaccurate); 

  mm.add_item("bootstrap chern_simons", BootstrapCS); 
  mm.add_item("bootstrap eta", BootEta); 

  // set up the menu of field types. 
  fm.add_item("invariant", 1); 
  fm.add_item("trace", 2); 
  fm.add_item("shape", 3); 
  fm.add_item("group_coeff", 4); 

  // set up the fudge_record while we're at it. 
  fudge_record.add_tag("EF");
  fudge_record.add_tag("F");
  fudge_record.add_tag("STD"); 
}

void env_A::validate_event(int& what)
{
  // if (!AG) AG = new alg_group; 
    
  alg_snap *m = get_AS();

  env_U::validate_event(what); 

  if (what < 1000) return; 

  if (what >= Holonomies && !m) {
    cout << "This function requires a triangulation.\n";
    what = Nothing; 
    return; 
  }

  if (what >= PrintManifold && !(m||G)) {
    cout << "This function requires a manifold or group.\n";
    what = Nothing; 
    return; 
  }

  if (what < PrintManifold) return; 

  // compute accurate info for all functions. 

  if (m) {
    if (!m->accurate_info_set()) {

      m->compute_accurate_info(0,0);

      // Force group to be recomputed. 
      clear_group(); 
    }

    if (what >= WordToAlgMatrix && what < Holonomies) {
      if (!get_group()) {
	what = Nothing; 
	return; 
      }
    }
  }
   
  // G is accurate when read in. 
  sync_AG();
}

void env_A::print_settings() const
{
  env_U::print_settings(); 

  printf("Simplify fundamental group: %s\n", simplify ? "yes":"no");
  printf("Stack space available in MB: %.2f of %.2f\n", 
	 stack_free()/1.E6, stack_size()/1.E6);
  printf("Significant digits printed: %ld\n", def_print.dec);
  printf("Show hilbert symbol factorization/embeddings: %s/%s\n", 
	 (hilb_print & 0x1) ? "yes" : "no", (hilb_print & 0x2) ? "yes" : "no"); 
  printf("Epsilon for accurate Dirichlet domain: %g\n", acc_vertex_eps);
#ifdef LONG_IS_64BIT
  printf("Precision in 64-bit words: %ld (= %ld decimal digits)\n", prec-2, digits); 
#else
  printf("Precision in 32-bit words: %ld (= %ld decimal digits)\n", prec-2, digits); 
#endif
  printf("Maximum degree of polynomials: %d\n", field::maxdeg); 
  printf("Fraction of significant digits to use in LLL: %.2f\n", field::fudge); 
  printf("Canonical poly mode: %d\n", cc_poly); 
  printf("Show Chern-Simons diagnostic: %s\n", cs_diagnostic ? "yes" : "no"); 
  printf("Field print style is: %s\n", print_fields_long ? "long" : "short"); 
  printf("Show field computation is: %s\n", show_field_comp ? "on" : "off"); 
}

pari get_surgery_coeffs_from_user(int n)
{
  if (n==0) return rvector(0); 

  char pr[100];
  string s;
  sprintf(pr, "Please enter %d pair(s) of surgery coefficients\n? ", n);
  get_input(s, pr, 0); 

  int cusp_num;
  pari pm, pl; 
  pari sc = rvector(2*n); 
  for(cusp_num=0; cusp_num < n; ++cusp_num) {
    if (!get_input(s, "? (surgery coefficients) ", 1)) return pZERO;
    pm = p_lisexpr(s.c_str()); 
    if (!get_input(s, "? (longitude coefficient) ", 1)) return pZERO;
    pl = p_lisexpr(s.c_str()); 
    if (pm.type()!=t_INT || pl.type()!=t_INT) {
      cout << "Can't handle non-integral surgeries yet\n";
      return pZERO;
    }
    sc[cusp_num] = pm; 
    sc[cusp_num + n] = pl;
  }
  return sc; 
}

pari get_pari_expr_from_user(const string& pr)
{
  string s; 
  get_input(s, pr, -1); /* get a whole line of input */ 
  if (s.length() > 0) {
    return p_flisexpr((char*)s.c_str()); 
  }
  return pZERO;
}

void convert_to_int_mat(pari const& mx, MatrixInt22 int_mat)
{
  int_mat[0][0] = mx[0][0].int_value();
  int_mat[0][1] = mx[1][0].int_value();
  int_mat[1][0] = mx[0][1].int_value();
  int_mat[1][1] = mx[1][1].int_value();
}

void gpwrite(const field& fld)
{
  string name;
  if (!get_input(name, "Filename to write ? ", 1)) return; 
  switchout((char*)name.c_str());
  fld.print(2); 
  switchout(0);
}

static void set_printform(pari_printform& pf)
{
  string s; 
  long d; 

  if (get_input(s, "significant digits required (-1=all) ? ")) {
    if (sscanf(s.c_str(), "%ld", &d)) {
      pf.dec = d; 
    } else {
      printf("number of significant digits not changed\n");
    }
  }
  if (get_input(s, "format (e=scientific, f=fixed, g=pari's choice) ? ")) {
    if (s[0] >= 'e' && s[0] <= 'g') {
      pf.format = s[0]; 
    } else {
      printf("format not changed\n");
    }
  }
#if 0
  if (get_input(s, "pretty or raw (0=pretty, 1=raw) ? ")) {
    if (sscanf(s.c_str(), "%d", &i) && i >= 0 && i <= 1) {
      pf.raw = i; 
    } else {
      printf("pretty or raw not changed\n");
    }
  }
  if (get_input(s, "integer field width (0=to fit) ? ")) {
    if (sscanf(s.c_str(), "%ld", &d)) {
      pf.field = d; 
    } else {
      printf("integer field width not changed\n");
    }
  }
#endif
}  

// This sets the flattening for a manifold even if none is found in the 
// eta_fudges file. 

static int load_eta_fudge(snap* m, const string& search_path, int report = 1)
{
  string_map info;
  char* name = get_triangulation_name(m->M()); 
  if (!load_record(search_path, "eta_fudges", fudge_record, name, info)) {
    if (report) printf("Warning: eta fudge not known for this manifold.\n"); 
    m->find_a_flattening(); 
    return 0; 
  }
  pari flattening = p_flisexpr(info["F"].c_str());

  /* check the equations are satisfied */ 
  pari equations = gtrans(get_complete_equations(m->M())); 
  pari sol = concat(flattening, pONE); 
  if (vecmax(abs(equations * sol))!=pZERO) {
    // Deal with invalid flattening. 
    printf("Invalid flattening found in eta fudge file.\n");
    m->find_a_flattening(); 
    return 0; 
  }

  m->set_flattening(flattening); 
  pari fudge = p_flisexpr(info["EF"].c_str())/integer(18); 
  m->set_eta_fudge(fudge);
  return 1; 
}

static string field_db_file(const alg_snap* m)
{
  /* decide which filename to use */ 
  string name; 
  if (all_cusps_are_filled(m->M()))
    name = "closed.fields"; 
  else if (all_cusps_are_complete(m->M()))
    name = "cusped.fields"; 
  else 
    name = "other.fields"; 
  return name; 
}

void env_A::update_precision()
{
  alg_snap* m = get_AS(); 
  if (m) {
    m->update_precision(); 
    clear_group(); 
    if (get_group()) sync_AG();
  } else {
    if (AG) AG->update_precision(); 
  }
}

void env_A::sync_AG()
{
  if (!G) {
    if (AG) {
      delete AG;
      AG = 0;
    }
    return;
  }

  if (!AG) AG = new alg_group;

  if (AG->group() != G) 
    AG->set_group(G);
}

void env_A::set_G(GroupPresentation* g)
{
  env_U::set_G(g);
  sync_AG();
}

bool env_A::read_group()
{
  if (!env_T::read_group()) return false; 
  polish(G, false, false);
  MoebiusTransformation c; 
  normalize(G, c, false, true); // try for a nice nearby normalization. 
  return true; 
}

void env_A::normalize_group(vector<FGWord> const& wl, bool report)
{
  MoebiusTransformation C; 
  if (!normalize(G, C, report, false, wl)) return; 
  conjugate(C); 
  if (AG) AG->clear_gf(); 
}

i_triangulation* env_A::new_manifold(Triangulation* tri) const
{
  alg_snap* m = new alg_snap(tri);
  SolutionType sol = get_filled_solution_type(m->M());
  if (sol < 1 || sol > 3) {
    cout << "Warning: solution type is " << solution_types[sol] << endl; 
  }
  m->make_orientable(); 
  choose_generators(m->M(), TRUE, FALSE); 
  load_eta_fudge(m, path, 0);
  return m; 
}

void env_A::eigenvalue_fields(bool canonical, int report, bool find_unit, double max_len) const
{
  list<named_field> inv_ev_fields;

  if (!num_geodesics()) {
    warn("eigenvalue_fields called but no geodesics computed\n"); 
    return; 
  }
  if (!AG || !AG->itfield().is_set()) {
    warn("eigenvalue_fields requires the invariant trace field\n"); 
    return;
  }

  const field& itfield(AG->itfield()); 

  int i = 0, j, n = num_geodesics();
  bool ok; 
  char buf[30]; 
  MoebiusTransformation mt; 
  pari tr, extr, ev, exev, tcp, ecp, denom, emp; 
  named_field new_ief; 
  list<named_field>::iterator ief; 
  vector<pari> bignf; 
  pari x_plus_xinv = p_lisexpr("x+x^-1");
  string fname; 
  pari rv0 = rvector(0); 

  Complex len;
  double eps = 1e-6; 
  int k = -1, ii=0; 

  while (ii < n) {

    // First print out complex length and list all equal (up to conjugacy)
    // length geodesics. 
    i = ii; 
    len = geodesic(i).length; 
    if (len.real > max_len) break; 
    k++; 
    if (report) cout << '(' << k << ')' << ' ' << len << ' ' << ii; 
    ii++; 
    while (ii < n) {
      if (fabs(geodesic(ii).length.real - len.real) > eps) break;
      if (fabs(geodesic(ii).length.imag - len.imag) < eps) {
	if (report) cout << ',' << ii; 
      } else if (fabs(geodesic(ii).length.imag + len.imag) < eps) {
	if (report) cout << ',' << ii << '*'; 
      } else break; 
      ii++; 
    }
    cout << endl; 

    mt = word_to_Moebius(G, geodesic(i).word);
    if (mt.acc_matrix.type()==1) {
      warn("problem computing accurate traces in eigenvalue_fields\n"); 
      break;
    }

    // First we simply compute and print out the trace 
    // of the square of this element, as an element of the invariant
    // trace field. 

    tr = (gsqr(gtrace(mt.acc_matrix)) - pTWO); 

    if (!itfield.contains(tr, extr)) {
      // Should never happen!
      if (report) {
	cout << "tsq(" << geodesic(i).word << ") = "; 
	tr.print();
      }
      continue;
    }

    if (report) {
      cout << "T(g^2) = "; lift(extr).print(); 
    }

#if 0
    // Keep going if this was the same as the last one. 
    if (extr == last_extr) {
      // Insert a ,i into the name of the field.
      fname = (*last_ief).name();
      sprintf(buf, ",%d", i); 
      fname.insert(fname.length()-1,buf); 
      (*last_ief).set_name(fname); 

      continue; 
    }
#endif

    // Now compute the eigenvalue and look for the eigenvalue field.
    ev = (tr + sqrt(tr * tr - pFOUR))/pTWO; 

    // First see if ev is in any of the fields we have found so far.
    for (j = 1, ief = inv_ev_fields.begin(); ief != inv_ev_fields.end(); j++, ief++) {
      if ((*ief).contains(ev, exev)) {

#if 0
	// Insert a ,i into the name of the field.
	fname = (*ief).name();
	sprintf(buf, ",%d", i); 
	fname.insert(fname.length()-1,buf);
	(*ief).set_name(fname); 
#endif

	if (report) {
	  cout << "IEF[" << k << "] = " << (*ief).name() << endl; 
	}

	// last_extr = extr; 
	// last_ief = ief; 
	break; 
      } else if ((*ief).contains(gconj(ev), exev)) {
	if (report) {
	  cout << "IEF[" << k << "] = " << (*ief).name() << '*' << endl; 
	}
	break; 
      }
    }

    if (ief == inv_ev_fields.end()) {

      // Did not find it in any field already computed so we have to make a
      // new field. 
      // The eigenvalue could actually be in the invariant trace field. 
      // Check that first: 

      if (itfield.contains(ev, exev)) {
	new_ief = itfield; 
	ok = true; 
	
      } else {

	// The eigenvalue lives in a proper extension of the inv. trace field. 
	// Try finding it algebraically: a minimum polynomial for the 
	// eigenvalue can be obtained from one for the trace. 

	tcp = caradj0(extr, 1);
	ecp = gsubst(tcp, 1, x_plus_xinv)*pow(polx[0], integer(lgeff(tcp)-1));
	denom = ggcd(ecp, deriv(ecp,0));
	emp = (ecp * denom[lgeff(denom)-1])/denom; 

	if (lgeff(emp)-1 == 2 * itfield.degree()) {
	  if (canonical) {
	    ok = new_ief.set_field(emp, ev, exev, 1); 
	  } else {
	    ok = new_ief.set_field(emp, ev, 0) && new_ief.contains(ev, exev); 
	  }
	} else {
	  // Could not find extension algebraically so do it numerically.
	  new_ief = itfield; 
	  ok = new_ief.extend(ev, exev); 
	}
      }

      if (!ok) {
	warn("problem finding an invariant eigenvalue field.\n");
	continue;
      }
    
      // Now we have the new field, name it, save it, and print it. 
      
      sprintf(buf, "IEF[%d]", k); 
      new_ief.set_name(string(buf)); 
    
      inv_ev_fields.push_back(new_ief); 
      ief = inv_ev_fields.end(); --ief; 

      if (report) {
	(*ief).print(0); cout << endl; 
      }

      if (find_unit) {
	bignf.push_back(bnfinit0((*ief).min_poly(), 0, rv0)); 
	j = bignf.size(); 
      }

      // last_extr = extr; 
      // last_ief = ief; 
    }

    cout << "E = "; lift(exev).print(); 
    if (find_unit) {
      if (j<1 || j>bignf.size()) {
	cout << "j out of range!!\n"; 
	continue;
      }
      pari u = isunit(bignf[j-1], exev); 
      cout << "Unit = "; u.print();
    }
    cout << endl; 
  }
  return;
}

static bool get_rational(pari const& re, pari& rat)
{
  if (numerical_zero(re)) { rat = pZERO; return true; }
  pari lin = algdep0(re, 1, (long)(1+digits*field::fudge));
  if (length(lin)!=2 || lin[1]==pZERO) return false; 
  rat = -lin[0]/lin[1];
  return numerical_zero(re-rat); 
}

static pari find_vanishing_factor(pari const& pol, pari const& root, bool report)
{
  pari factors = factor(pol); 

  if (report) cout << "  Finding vanishing factor:\n";

  pari eps, mineps=real(1); 
  pari f; 
  int i, nf = length(factors[0]), best_i=0; 
  for (i=0; i<nf; ++i) {
    f = factors[0][i]; 
    eps = gabs(gsubst(f,0,root)/gsubst(deriv(f,0),0,root)); 
    if (eps < mineps) {
      mineps = eps; 
      best_i = i; 
    }

    if (report) {
      cout << "  "; factors[0][i].print(); 
      cout << "  "; eps.print(); 
    }

  }

  if (!numerical_zero(mineps)) {
    cout << "  Insufficient accuracy in find_vanishing_factor\n";
    return pZERO;
  }

  return factors[0][best_i]; 
}

static pari get_min_poly(pari const& z, bool report)
{
  pari pol = algdep0(z, field::maxdeg, (long)(1+digits*field::fudge)); 
  pari mp = find_vanishing_factor(pol, z, report); 
  return mp;
}


static pari min_poly_z_on_zbar(pari z, pari const& mpz, bool report)
{
  // Below we assume z != conj(z), therefore.. 
  if (numerical_zero(gimag(z))) return p_lisexpr("x-1"); 

  pari x(polx[0]);
  pari y(polx[1]);

  pari res = polresultant0(gsubst(mpz, 0, x*y), gsubst(mpz, 0, y), 1, 0); 

  // Make z more accurate (double the precision). 
  pari rl = roots(mpz, 2 * prec);
  int i, n = length(rl); 
  for (i=0; i<n; ++i) {
    if (numerical_zero(rl[i] - z)) {
      z = rl[i]; 
      break; 
    }
  }

  if (i==n) {
    cout << "Problem trying to increase the precision of the given root\n"; 
  }

#if 0
  // Find a polynomial for z/conj(z).
  pari pol = pONE;
  pari lin = polynomial(1); 
  lin[1] = pONE; 
  for (i=0; i<n; ++i) {
    for (j=0; j<n; ++j) { 
      if (i==j) continue; // Assumes z!=conj(z), skips factors of x-1. 
      lin[0] = -rl[i]/rl[j]; 
      pol *= lin; 
    }
  }

  // Make coeffs rational. 
  pari rat; 
  n = length(pol); 
  for (i=0; i<n; ++i) {
    if (!get_rational(pol[i], rat)) break; 
    pol[i] = rat; 
  }

  pari mp0; // mp0 is initially 0. 

  if (i==n) { // If found rational approximations for all the coefficients. 
    mp0 = find_vanishing_factor(pol, z/gconj(z), report); 
  }


  if (mp0.type()==t_INT) {
    cout << "Original min_poly_z_on_zbar fails\n";
  } else {
    mp0 /= content(mp0); 
  }

  // cout << "Resultant:\n";
  // res.print(); 

  if (mp0.type()!=t_INT) {
    if (mp0 != mp) {
      cout << "DIFFERENT RESULTS FROM THE TWO METHODS\n";
      cout << "Original method gives:\n"; 
      mp0.print(); 
    }
  }
#endif

  pari mp = find_vanishing_factor(res, z/gconj(z), report); 
  if (mp.type()==t_INT) return pZERO; 

  mp /= content(mp); 

  return mp; 
}

static void get_wordlist(string const& prompt, vector<FGWord>& wl)
{
  string s; 
  if (get_input(s, prompt, -1)) {
    vector<string> words;
    int i, n = split(s, words); 
    wl.resize(n); 
    for (i=0; i<n; ++i) 
      wl[i] = FGWord(words[i].c_str()); 
  } 
}

// test if v == u^n for some integer n. 
// 1=yes, v == u^n
// 0=no,  v != u^n
// <0 =unknown. 
//  -1: Q(u) field info not found
//  -2: LLL failed to find Q(u,v) (probably not a power)
//  -3: two non-units, test not implemented

static int test_is_power(pari const& u, pari const& mpu,
			 pari const& v, pari const& mpv, bool report)
{
  // If so Q(v) \subseteq Q(u) so deg(mpv) <= deg(mpu);

  int du=lgeff(mpu)-1, dv=lgeff(mpv)-1;
  if (dv > du) {
    if (report) cout << "Degree of Q(v) is bigger\n";
    return 0; 
  }
  if (du % dv != 0) {
    if (report) cout << "Degree of Q(v) does not divide degree of Q(u)\n";
    return 0; 
  }
  if (nfisincl(mpv, mpu).type()==t_INT) {
    if (report) cout << "Q(u) does not contain Q(v)\n";
    return 0; 
  }

  bool uint = (mpu[du]==pONE);
  bool vint = (mpv[dv]==pONE);
  if (uint != vint) {
    if (report) cout << "Have one integer and one non-integer\n";
    return 0;
  }

  bool uunit = (uint && gabs(mpu[0])==pONE);
  bool vunit = (vint && gabs(mpv[0])==pONE);
  if (uunit != vunit) {
    if (report) cout << "Have one unit and one non-unit\n";
    return 0; 
  }

  // Ensure that we can find a field of the required degree.
  int mdeg = field::maxdeg; 
  if (du > field::maxdeg)
    field::maxdeg = du; 

  field Qu;
  pari xu, xv;
  if (Qu.extend(u,xu,2)==0 &&
      Qu.extend(u,xu,0)==0) {
    if (report) cout << "Problem finding field info for Q(u)\n";
    field::maxdeg = mdeg; 
    return -1; 
  }

  int res = Qu.extend(v,xv,0);

  field::maxdeg = mdeg; // Restore original field degree limit. 

  if (res==2) { // field containing v is bigger
    if (report) cout << "LLL says v is not in Q(u) because Q(u,v) is bigger\n";
    return 0; 
  } 

  if (res==0) {
    if (report) {
      cout << "LLL failed to find Q(u,v)\n";
      cout << "This probably means v is not a power of u\n";
    }
    return -2; 
  }

  // res==1, i.e. v belongs to Q(u).

  if (report) {
    cout << "LLL says v belongs to Q(u)\n"; 
    cout << "Q(u): "; Qu.print(0); cout << endl; 
    cout << "u = "; lift(xu).print(); 
    cout << "v = "; lift(xv).print(); 
  }

  if (!uunit) {
    if (report) cout << "Two non-units in Q(u), power test not implemented\n";
    return -3;
  }
    
  pari rv0 = rvector(0);
  pari bnf = bnfinit0(Qu.min_poly(),0,rv0);

  pari uub = isunit(bnf,xu);
  pari vub = isunit(bnf,xv);
  if (report) {
    cout << "unit(u) = "; uub.print();
    cout << "unit(v) = "; vub.print(); 
  }
  int i, r = length(uub)-1;
  for (i=0; i<r; ++i)
    if (uub[i]!=pZERO) break;
  bool ucyclotomic = (i==r);
  if (!ucyclotomic) {
    
    // uub, vub are something like [-2, 4, Mod(0, 2)]
    // and we want to see if v is an integer multiple of u. 
    
    if (vub[i] % gabs(uub[i]) != pZERO) 
      return 0; 
    pari n = vub[i]/uub[i];
    pari check = vub - n*uub;
    for (i=0; i<r; ++i)
      if (check[i]!=pZERO) return 0; 
    
    return (lift(check[r])==pZERO) ? 1 : 0; 
    
  } 
  
  // u is cyclotomic.
  
  for (i=0; i<r; ++i)
    if (vub[i]!=pZERO) return 0; // u is cyclotomic, v isn't. 
  
  // uub, vub are of the form [0, 0,...,Mod(a,m)]
  if (report) cout << "Both are roots of unity\n";
  
  int m = uub[r][0].int_value();
  
  // Check if vub[r] is multiple of uub[r] (mod m)..
  //  surely there is a neater way?
  for (i=0; i<m; ++i)
    if (lift(integer(i)*uub[r]-vub[r])==pZERO) return 1;
  
  return 0; 
}

// An ortholine is significant if each end is 
// either on the spiraling geodesic in the interval
// e, or on a closed geodesic. 

static bool significant_ortholine(Ortholine const& o, const double e[2])
{
  const double eps=1e-3; 
  int i;
  double x; 
  for (i=0; i<2; ++i) {
    if (o.geodesic_num(i) != -3) continue; // not on spiraling geodesic. 
    x = o.position(i).real;
    if (x < e[0]-eps || x > e[1]+eps) return false; // not in interval. 
  }
  return true; // in the interval 
}


static void find_interval_outside_tube(const vector<line>& l, const O31_line& L0, double t_rad, double e[2])
{
  // Print out the points where the spiraling geodesic L0
  // leaves the tube around l[0] for the first time and enters 
  // the tube around l[1] for the last time. 

  // For the second case where L0 = [infinity,1]
  // these would be the points where L0 leaves the
  // tube around l[0] and later enters the tube around l[1]. 
  
  // Given two parallel geodesics l0, l and points p,q
  // on l0 distance t apart, (with d(q,l) < d(p,l) -- i.e.
  // q is closer to the common endpoint than p), we have 
  // the formula e^t = sinh(d(p,l))/sinh(d(q,l)). 
  
  // The line L0 has its own (arbitrary) coordinate system in which 
  // the real parameter of L0.point() (on the line) is zero. 
  // We use the above formula to find the point p on L0, in this 
  // coordinate system where d(p, l[1]) = t_rad, in terms of 
  // d(L0.point(), l[1]), and t_rad. 

  // e records the parameter values [leaving l[0] tube, entering l[1] tube]. 
  int i; 
  double dist;
  for (i=0; i<2; ++i) {
    dist = distance_to_org(O31_line(l[i]), L0.point());
    e[i] = (i?1.:-1.) * log(sinh(dist)/sinh(t_rad));
  }
}

static int sally_test_arithmetic(const pari& z, const pari& u, const pari& mpu, int report)
{
  pari mp = get_min_poly(z, report > 2); 
  
  if (report > 1) cout << "MP(Z): ";
  if (mp.type()==t_INT) { 
    if (report > 1) cout << "not found\n"; 
    return -5; 
  }
  if (report > 1) mp.print(); 

  pari v = z/gconj(z); 
  pari mpv = min_poly_z_on_zbar(z, mp, report > 2);

  if (report > 1) cout << "MP(v = Z/conj(Z)): "; 

  if (mpv.type()!=t_POL) {
    if (report > 1) cout << "not found\n";
    return -5; 
  }
  if (report > 1) mpv.print(); 

  // Now we try to determine if v=z/conj(z) 
  // might be a power of u=lam^2/conj(lam^2).
  // If all's well so far we have minimum 
  // polynomials for u and v in mpu and mpv. 

  return test_is_power(u, mpu, v, mpv, report > 1);
}

/* 
 * The purpose of this function is to test conditions 
 * discovered by Sally Kuhlmann for the manifold to contain
 * infinitely many simple geodesics. Given a geodesic in
 * the manifold we have to check two things: an arithmetic
 * condition satisfied by the complex length of the geodesic
 * and its smallest orthodistance (half of whose real part
 * is the injectivity radius of the geodesic); and a geometric
 * condition on a "spiraling" geodesic, one whose lift 
 */

void env_A::sally(int gn, int report)
{
  if (gn < 0) {
    cout << "Invalid geodesic number\n"; 
    return; 
  }

  vector<FGWord> wl(2);
  if (!compute_geodesics(gn+1)) {
    cout << "Problem computing geodesic " << gn << endl; 
    return; 
  }

  vector<Ortholine> oll; 
  vector<int> vgn(1); 
  double rad = 1.; 
  vgn[0] = gn; 
  while (!oll.size() && rad < 4.1) {
    oll = ortholines(vgn, rad);
    rad += .5; 
  }
  if (!oll.size()) {
    cout << "Problem computing smallest orthodistance\n"; 
    return; 
  }

  // oll[0] is the shortest ortholine. 

  // wl[0] is word for (a lift) geodesic gn. 
  // wl[1] is word for a nearest lift of geodesic gn.  
  wl[0] = geodesics[gn].word;
  wl[1] = oll[0].word * wl[0] * inverse(oll[0].word); 

  // Set t_rad to the radius of the maximal tube. 
  double t_rad = oll[0].distance().real/2.; 

  // Print length of geodesic and tube radius. 
  if (report) {
    cout << "Geodesic: " << gn << ", length:" << geodesics[gn].length;
    cout << ", tube radius: " << t_rad << endl; 
  }

  // Print words for wl[0] and wl[1].
  if (report > 1) cout << "Words: " << PSeq(wl) << endl; 

  // Conjugate the fundamental group such that 
  // ends of geodesic wl[0] are at [0,infinity] an an end of
  // wl[1] is at 1. 
  normalize_group(wl, false); 

  // Print the Moebius transformations of wl[0] and wl[1]
  // in the normalized fundamental group, and also the fixed 
  // points. l[0] will store the fixed points of wl[0] in
  // the order (repelling, attracting). gens[0] stores the
  // Moebius transformation (SL(2,C) matrix) for wl[0]. 
  int i, n_gens = wl.size(); 
  vector<MoebiusTransformation> gens(n_gens); 
  vector<line> l(n_gens); 
  for (i=0; i<n_gens; ++i) { 
    gens[i] = word_to_Moebius(G, wl[i]); 
    fixed_points(gens[i].matrix,l[i]); 

    if (report < 2) continue; 
    cout << gens[i] << endl; 
    cout << l[i] << endl; 
  }

  // THE ARITHMETIC CONDITION

  // lambda_sq = square of eigenvalue for wl[0]. 
  pari u, mp, mpu, lambda_sq = gsqr(gens[0].acc_matrix[0][0]); 

  // get the invariant trace field (subfield gen by gens[0]) degree. 
  pari itf_mp = get_min_poly(lambda_sq + ginv(lambda_sq), false);
  int deg = 0; 
  int olddeg = field::maxdeg;

  // increase the degree temporarily if necessary. 
  if (itf_mp.type()!=t_INT) {
    deg = lgeff(itf_mp)-1;
    cout << "Inv trace field degree = " << deg << endl; 

    if (2*deg > olddeg) 
      field::maxdeg = 2*deg; 
  }

  // find a minimum polynomial satisfied by lambda^2
  mp = get_min_poly(lambda_sq, report > 2);
  if (report > 1) {
    cout << "MP(Lam^2): "; 
    if (mp.type()!=t_INT) mp.print();
    else cout << "not found\n";
  }

  if (mp.type()!=t_INT) { 

    // u is lambda^2 divided by its complex conjugate. 
    u = lambda_sq/gconj(lambda_sq);
    mpu = min_poly_z_on_zbar(lambda_sq, mp, report > 2); 

    if (report > 1) {
      cout << "MP(u = Lam^2/conj(Lam^2)): "; 
      mpu.print(); 
    }

  } 
  
  // oll[0].word was the word conjugating wl[0] to wl[1]. 
  // we compute its Moebius transformation. [a,b;c,d]
  pari mx = word_to_Moebius(G, oll[0].word).acc_matrix;

  pari z;
  int res[2]; 
  int spir; 

  int k, tc_ok[2]; 
  O31_line L0, L1;
  line l0; 
  double e[2]; 

  for (spir=0; spir<2; ++spir) {

    if (spir==0) {
      // Create the line from Zero (repelling for l[0]) to attracting of l[1]. 
      l0 = line(Zero,l[1].end[1]);
    } else {
      // Create the geodesic [infinity, 1], attracting of l[0] to repelling of l[1]. 
      l0 = line(Infinity,One); 
    }
    if (report) cout << "\nSpiraling geodesic: " << l0 << endl; 


    // THE ARITHMETIC CONDITION

    if (mpu.type() != t_POL) {
      res[spir] = -4; // Failed to find MP(u).

    } else {
      if (spir==0) {
	z = -mx[0][1]*mx[1][0]/(mx[0][0]*mx[0][0]);
      } else {
	z = -mx[1][1]*mx[1][1]/(mx[1][0]*mx[0][1]); 
      }

      // Return value is 
      // 0 if definitely not a power
      // 1 if definitely a power
      // <0 for all inconclusive results (see test_is_power). 

      res[spir] = sally_test_arithmetic(z, u, mpu, report);
    }

    if (report) {
      switch (res[spir]) {
      case 0:
	cout << "Arithmetic test satisfied, spiral ends are disjoint\n";
	break; 
      case 1:
	cout << "Arithmetic test fails, spiral ends may intersect\n";
	break; 
      case -2:
      case -5:
	cout << "Arithmetic test probably satisfied (" << res[spir] << "), ends probably disjoint\n";
	break; 
      default:
	cout << "Arithmetic test inconclusive (" << res[spir] << ")\n";
	break; 
      }
    }

    // Now look at the portion of the spiraling geodesic l0
    // lying outside the maximal tube around our geodesic. 

    L0 = O31_line(l0); 

    find_interval_outside_tube(l, L0, t_rad, e); 

    if (report > 1) cout << "Spiraling interval: [" << e[0] << ',' << e[1] << "]\n";

  
    // We want to know, in the manifold, whether the above segment 
    // meets the original geodesic, and whether it meets itself, 
    // including either of its spiraling "arms". 

    // To do this we will essentially be computing some ortholines
    // and checking whether any of them have length zero. 
    // So far we have been working in a certain normalized 
    // conjugate of the fundamental group (image of holonomy rep.)
    // The ortholine code needs to use a Dirichlet domain, which
    // was computed with respect to a different conjugate of the 
    // fundamental group. Therefore we now have to start working
    // in Dirichlet domain coordinates. 

    L1 = inverse(domain->conjugacy) * L0; 

    // Create the interval. (It doesn't have any fundamental group
    // information, hence the empty FGWord().) We label it -3. 
    // (Closed geodesics are labelled, 0,1,...)

    interval Spir(L1,e[0],e[1],FGWord(),-3);
    list<interval> ivls; 

    // Chop the interval up into pieces inside one copy 
    // of the Dirichlet domain. Also get the chopped up pieces
    // of the closed geodesic. 
    all_crossing_lifts(domain, Spir, ivls, false); 
    get_crossing_lifts(domain, GSpec(gn, &geodesics[gn]), ivls, false); 

    // Compute the closest approaches between the interval and 
    // itself, the interval and the closed geodesic (and the closed
    // geodesic to itself, which we're not interested in.) 
    Ortholine_set os; 
    compute_ortholines(domain, tiling, ivls, t_rad+.1, os, false); 

    // cout << PSeq(os,"\n") << endl; 

    // First look for any places where the interval on the spiraling
    // geodesic enters the tube. 
    Ortholine_set::const_iterator ol;
    double small_t_rad = t_rad; 
    for (ol=os.begin(); ol!=os.end(); ++ol) {
      if (ol->geodesic_num(0)==-3 && ol->geodesic_num(1)==gn && 
	  ol->distance().real < small_t_rad)
	small_t_rad = ol->distance().real;
    }

    if (small_t_rad < t_rad) { 
      if (report > 1) cout << "Interval enters tube, closest approach: " << small_t_rad << endl; 

      // Find a bigger interval on the spiraling geodesic. 
      // This will ensure that any time the original interval
      // met either of the spiraling ends this new bigger 
      // interval will find that point as a point of self
      // intersection. 

      find_interval_outside_tube(l, L0, small_t_rad, e); 

      if (report > 1) cout << "Expanded interval: [" << e[0] << ',' << e[1] << "]\n";

      Spir = interval(L1,e[0],e[1],FGWord(),-3);

      // Get the chopped up interval & geodesic again. 
      ivls.clear(); 
      all_crossing_lifts(domain, Spir, ivls, false); 
      get_crossing_lifts(domain, GSpec(gn, &geodesics[gn]), ivls, false); 
    }

    // Look for any of the kinds of intersection that
    // we are interested in: intersections between the 
    // interval and the closed geodesic or between the interval
    // and itself (or possibly between the closed geodesic
    // and itself if it is not the shortest geodesic say.) 

    tc_ok[spir] = 1; 
    double eps=1e-3, dist, min_dist = 2*t_rad;
    Ortholine min_ort; 
    bool found_significant = false;
    for (ol=os.begin(); ol!=os.end(); ++ol) {

      // Each party to the crossing has to be either the 
      // closed geodesic, or the spiraling geodesic with 
      // the crossing point being inside the interval. 

      if (significant_ortholine(*ol,e)) {
	dist = ol->distance().real; 
	if (dist < eps) tc_ok[k] = 0; 

	// Let's save the smallest orthodistance we find. 
	if (dist < min_dist) {
	  found_significant = true; 
	  min_dist = dist; 
	  min_ort = *ol; 
	}
      }
    }

    // Print the shortest significant orthodistance we found. 
    if (found_significant && report > 1) {
      cout << "Smallest orthodistance:\n";
      cout << min_ort << endl;
    } 
    if (report) {
      if (tc_ok[spir]) cout << "No intersections on interval (distance bound " << min_dist << ")\n";
      else cout << "Intersection found on interval\n";
    }
  }

  field::maxdeg = olddeg; 

  if (!save_sally_results) 
    return; 


  FILE* slog = fopen("sally_results", "a");
  if (slog) { 
    fprintf(slog, "%5d %3d %2d %2d %2d %2d\n", closed_num, gn, res[0], tc_ok[0], res[1], tc_ok[1]);  
    fclose(slog); 
  } else {
    cout << "Problem writing record to log file!!\n"; 
  }

}

#if 0
// deleted from sally() above. 

  vector<string> names; 
  pari gf_gens, ex_gf_gens; 
  get_groupfield_gens(gens, n_gens, names, gf_gens); 
  
  // try to find the group coefficient field. 
  field gf;
  int res = gf.generated_by(gf_gens, ex_gf_gens, 2); 
  if (res) {
    cout << "Field\n";
    gf.print(0);
    cout << endl; 
  } else {
    cout << "Group coefficient field not found\n"; 
  }
  if (wl.size() < 2) {
    delete[] gens; 
    return; 
  }
#endif

void env_A::process_event(int what)
{
  alg_snap *m = get_AS();

  string s; 
  int j; 

  switch (what) {
  case PrintAccurateShapes:
    m->acc_shapes().print(); 
    break;

  case AccurateLogShapes:
    m->accurate_log_shapes().print(); 
    break;

  case PrintExactShapes:
    if (m->alg_shapes().type()==t_INT) {
      cout << "Please load or compute the shape field first.\n"; 
      break; 
    }
    lift(m->alg_shapes()).print(); 
    m->tetfield().print(0); cout << "\n"; 
    break;

  case Holonomies:
    holonomies(m->M()).print(); 
    break;


#if 0
  case PrintIsometries:
    m->print_isometries();
    break;
  case ClearIsometries:
    m->clear_isometries(); 
    break;
#endif
  case PrintFlattening:
    {
      pari fl = concat(m->the_flattening(), pONE);
      pretty_print(to_int_matrix(fl)); 
      pari equations = gtrans(get_complete_equations(m->M())); 
      cout << "Complete equations:\n"; 
      show_matprod(equations, fl); 
    }
    break;

  case ComplexVolume2:
    {
      vector<int> e_or; 
      bool has_acy = seek_acyclic_edge_orientations(T->M(), e_or);
      if (has_acy) {
	char pr[100]; 
	sprintf(pr, "Edge orientations to use\n? "); 
	pari new_eor = get_pari_expr_from_user(pr);
	if (new_eor.type()!=t_INT) {
	  int j;
	  for (j=0; j<e_or.size(); ++j) 
	    e_or[j] = new_eor[j].int_value(); 
	}
      }

      complex_volume3(T->M(), e_or);

    }
    break;

  case PrintManifold:
    if (get_group(false)) sync_AG();
    if (AG) {
      AG->set_tracefield_gens(); 
      AG->print(); 
    }
    if (m) m->print(); 
    break;

  case PariEquations:
    gtrans(get_complete_equations(m->M())).print(); 
    break;
  case PariFillingEquations:
    gtrans(get_filled_equations(m->M())).print(); 
    break; 

  case WordToAlgMatrix:
    {
      FGWord wd;
      if (!get_word(wd, geodesics)) break; 
      if (AG) lift(AG->alg_group_element(wd)).print(); 
    }
    break;

  case EvaluateTrace:
    {
      FGWord wd;
      if (!get_word(wd, geodesics)) break; 
      if (!AG) break;
      pari xtr = AG->exact_trace(wd);
      if (xtr.type()==t_POLMOD) lift(xtr).print(); 
    }
    break; 

  case Polish:
    m->compute_accurate_info(1,1); 
    break; 

  case PolishRep:
    {
      if (!get_group()) break; 
      bool filled = false; 
      bool find_field=true;

      if (T) {
	int nf = num_filled_cusps(T->M());
	filled = (nf==get_num_cusps(T->M()));
	if (nf!=0 && !filled) {
	  cout << "Warning: parabolicity not enforced at unfilled cusps!\n";
	  cout << "Use the \"fill\" function first, for accurate results.\n";
	  find_field=false; 
	}
      }

      if (!polish(G, filled, true)) break; 
      if (!find_field) break; 

      // AG->compute_fields(cc_poly, show_field_comp | 2, 0x1); 
    }
    break; 

  case Status:
    {
      cout << "Triangulation: " << T << endl; 
      if (T) cout << "  num tet: " << get_num_tetrahedra(T->M()) << endl; 
      cout << "Group: " << G << endl;
      if (G) cout << "  num gens: " << fg_get_num_generators(G) << endl; 
      cout << "Domain: " << domain << endl;
      if (domain) {
	cout << "  volume: " << domain->approximate_volume << endl; 
	cout << "  conj ht.: " << domain->conjugacy(0,0) << endl; 
	cout << "  tiling rad: " << tiling.radius << endl; 
	cout << "  n geod: " << geodesics.size() << endl; 
	cout << "  n symms: " << symmetries.size() << endl; 
      }
      int nf = the_tube().num_faces(); 
      cout << "Tube: " << nf << (nf ? " faces":"") << endl; 
      if (nf)
	cout << "  geod num.: " << geodesic_number() << endl; 

      if (num_core_ortholines())
	cout << "  n core ort.: " << num_core_ortholines() << endl; 

      cout << "Alg. Group: " << AG << endl; 
      if (AG) {
	cout << "  group: " << AG->group() << endl; 
	pari eg = AG->get_exact_gens(); 
	cout << "  exact gens: ";
	if (eg.type()==t_INT) cout << 0 << endl; 
	else cout << length(eg) << endl; 
      }
    }
    break; 

  case TestMats:
    testmats(); 
    break;

  case NormalizeGroup:
    {
      vector<FGWord> wl;
      get_wordlist("Elements to normalize ? ", wl);
      normalize_group(wl, true); 
    }
    break; 

  case NormalizeParabolics:
    {
      // get the list of parabolic words from the group
      vector<FGWord> wl;
      int i, n = fg_get_num_cusps(G); 
      for (i=0; i<n; i++) {
	wl.push_back(fg_meridian(G,i));
      }
      normalize_group(wl, true);
    }
  break;
  case SetSallyVerbosity:
    if (!get_input(s, "verbosity level ? ")) break; 
    if (!sscanf(s.c_str(), "%d", &sally_verbosity)) break; 
    if (sally_verbosity < 0) sally_verbosity = 0; 
    break; 

  case SetSaveSally:
    save_sally_results = ask("append results to sally_results?", 1);
    break; 

  case SallysFunction:
    {
      if (!get_Dirichlet_domain()) break; 

      if (!get_input(s,"cutoff radius (or 'g') ? ")) break; 

      // do a single geodesic.
      if (s=="g") {
	if (!get_input(s, "geodesic number ? ", 1)) break; 
	int gn;
	if (sscanf(s.c_str(), "%d", &gn)!=1) break; 
	sally(gn, sally_verbosity); 
	break; 
      }

      if (s=="c") {
	if (!compute_geodesics(2.2)) {
	  cout << "Problem computing geodesics\n";
	  break;
	}
	Complex clen;
	int gn, ori; 
	FGWord conj; 
	if (!find_core_geodesic(0, clen, gn, ori, conj)) {
	  cout << "Problem finding the core geodesic\n"; 
	  break; 
	}
	sally(gn, sally_verbosity); 
	break; 
      }

      double len;
      if (sscanf(s.c_str(), "%lf", &len) != 1) break; 
      
      // do symmetry orbits up to given length
      vector<geo_orbit> orbits; 
      if (!symmetry_orbits(orbits, 0, 0, len, 2)) break;
      int k, no=orbits.size();
      Ortholine o; 
      for (k=0; k<no; ++k) {
	if (k) cout << endl; 
	sally(orbits[k].gnum[0], sally_verbosity); 
      }
      env_D::clear(); // normalizations make this stuff inaccurate

    }
    break; 

  case AccurateDirichletDomain:
    {
      if (!polish(G, num_cusps()==0, true)) break; 

      if (AG) AG->set_tracefield_gens(); 
      cout << "Accuracy: " << accuracy(G) << endl; 
      compute_domain(G, &domain);
    }
    break; 

#if 0
  case Inaccurate:
    m->compute_accurate_info(2,1); 
    break;
#endif
    // Things relating to NUMBER FIELDS. 
  case ComputeInvariantTraceField:
    if (AG) AG->compute_fields(cc_poly, show_field_comp | 2, 0x1); 
    break;
  case ComputeTraceField:
    if (AG) AG->compute_fields(cc_poly, show_field_comp | 2, 0x2); 
    break;
  case ComputeGroupField:
    if (AG) AG->compute_fields(cc_poly, show_field_comp | 2, 0x8); 
    break;
  case ComputeShapeField:
    m->compute_fields(cc_poly, show_field_comp | 2, 0x4); 
    break;
  case ComputeCuspFields:
    m->compute_fields(cc_poly, show_field_comp | 2, 0x10); 
    break;

  case ComputeFields:
    if (get_group(false)) sync_AG();
    if (AG) {
      if (ask("Invariant trace field?",1))
	AG->compute_fields(cc_poly, show_field_comp | 2, 0x1); 
      if (ask("Trace field?",1)) 
	AG->compute_fields(cc_poly, show_field_comp | 2, 0x2); 
      if (ask("Group coeff field?",1))
	AG->compute_fields(cc_poly, show_field_comp | 2, 0x8); 
    }
    if (!m) break; 
    if (ask("Shape field?",1)) 
      m->compute_fields(cc_poly, show_field_comp | 2, 0x4);

    if (!m->is_closed() && ask("Cusp fields?",1)) 
      m->compute_fields(cc_poly, show_field_comp | 2, 0x10); 

    break;

#if 0
  case SaveFields:
    {
      string name = field_db_file(m); 
      string file = choose_save_file(".", name); 
      if (!file.length()) break; 

      m->save_fields(file.c_str());
    }
    break; 
  case LoadFields:
    m->load_fields(path.c_str(), field_db_file(m).c_str());
    break; 
#endif
      
  case WriteField:
    switch (fm.read_input("?")) {
    case 1:
      if (AG) gpwrite(AG->itfield());
      break; 
    case 2:
      if (AG) gpwrite(AG->tracefield());
      break;
    case 3:
      if (AG) gpwrite(AG->groupfield()); 
      break;
    case 4:
      if (m) gpwrite(m->tetfield());
      break;
    default:
      printf("invalid field name\n"); 
    }
    break; 

  case PrintFields:
    if (AG) AG->print_fields(print_fields_long);
    if (m) m->print_fields(print_fields_long);
    break;

  case ComputeEigenvalueFields:
    {
      if (!get_Dirichlet_domain()) break; 
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      if (!compute_geodesics(radius)) break;

      eigenvalue_fields(1,1,0,radius); 
      // (1,1,0) = (get canonical polys, report on it, don't find units).
    }
    break; 

  case ComputeUnits:
    {
      if (!get_Dirichlet_domain()) break; 
      double radius;
      if (!get_input(s,"cutoff radius ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &radius) != 1) break; 
      if (!compute_geodesics(radius)) break;

      eigenvalue_fields(1,1,1,radius); 
      // (1,1,1) = (get canonical polys, report on it, find units).
    }
    break; 

  case VerifyEquations:
    if (m) m->verify_and_report(); 
    if (AG) AG->verify(); 
    break;

  case TestArithmeticTraces:
    {
      // bool cusped = get_num_cusps(m->M())!=num_filled_cusps(m->M());
      if (!AG) break;
      AG->set_tracefield_gens(); 
      if (AG->imquad_integer_itgens()) {
	cout << "Manifold appears to be arithmetic.\n";
      } else {
	cout << "Traces squared are not quadratic integers\n"; 
      }
    }
    break; 

  case PrintInvariants:
    {
      printf("Volume: "); m->accurate_volume().print(); 

      int corr;
      pari cs = m->accurate_chern_simons(corr); 
      if (corr == 0) {
	printf("Chern-Simons (mod 1/24): "); 
      } else if (corr == 1) {
	printf("Chern-Simons (mod 1/2): "); 
      } else {
	printf("Chern-Simons (mod 1): "); 
      }
      cs.print(); 

      if (!m->set_eta_info()) {
	printf("Eta invariant: unknown\n"); 
      } else {
	printf("Eta invariant: "); 
	m->eta_invariant().print(); 
      }
    }

    if (!AG || !AG->itfield().is_set()) break; 
    /* else FALL THROUGH to PrintArithmetic */

  case PrintArithmetic:
    if (!AG || !AG->itfield().is_set()) {
      printf("Arithmetic invariants not available: invariant trace field required.\n"); 
      break; 
    }

    printf("Invariant trace field: "); 
    AG->itfield().print(0); // Terse format. 
    printf("\n"); 
    printf("Integer traces: "); 
    switch (AG->integer_traces()) {
    case 0:
      printf("NO\n"); 
      break;
    case 1:
      printf("YES\n"); 
      break;
    default: /* -1 */ 
      printf("unknown (this should never happen)\n"); 
    }

    // Force computation of hilbert symbol. 
    printf("Invariant quaternion algebra: "); 
    if (!AG->compute_hilbert_symbol()) {
      printf("unknown\nproblem computing the Hilbert symbol.\n");
    } else {
      AG->print_hilbert_symbol(hilb_print); 

      if (m && get_num_cusps(m->M()) != num_filled_cusps(m->M())) {
	printf("Ramification: none (manifold has cusps)\n");
      } else {
	printf("Real ramification: "); 
	AG->real_ramification().print();

	printf("Finite ramification: "); 
	pari pl = ramification(AG->hilbert_symbol(), AG->itfield().p_nf(), 
			       AG->integer_traces(), 0 /* no report */, 30);

	if (length(pl)==0) printf("none"); 
	printf("\n"); 
	for (int i=0;i<length(pl);i++)
	  pl[i].print(); // Print each place on a separate line. 
      }

      printf("Arithmetic: "); 
      switch (AG->arithmetic()) {
      case 0:
	printf("NO\n"); 
	break;
      case 1:
	printf("YES\n"); 
	break;
      default: /* -1 */ 
	printf("unknown (this should never happen)\n"); 
      }
    }

    printf("Borel Regulator: "); 
    if (m && m->tetfield().is_set() && AG->itfield().is_set()) {
      m->borel_regulator(AG->itfield()).print(); 
    } else {
      printf("not computed (shape field required)\n"); 
    }
    break; 

  case PrintBorel:
    {
      if (!m->tetfield().is_set() || !AG || !AG->itfield().is_set()) {
	printf("Shape field and invt. trace field are required for this function.\n"); 
	break; 
      }

      char pr[100]; 
      sprintf(pr, "Field to use in format [min-poly, canonical-root-num]\n? "); 
      pari fspec = get_pari_expr_from_user(pr);

      if (fspec.type() == t_INT) {
	// No field supplied, just do default (invt. trace field). 
	printf("Borel regulator: ");
	m->borel_regulator(AG->itfield()).print();
	break;
      }

      field F(fspec);
      pari rt = AG->itfield().root(), xr; 
      if (!F.contains(rt)) {
	printf("Field doesn't contain invariant trace field, computing join...\n"); 
	if (!F.extend(rt + F.root(), xr, cc_poly, show_field_comp)) {
	  printf("Sorry, unable to compute join.\n");
	  break; 
	}
	F.print(0); printf("\n"); 
      }
      printf("Borel regulator: ");
      m->borel_regulator(AG->itfield(), &F).print();
    }     
    break;

  case PrintCuspInfo:
    if (!m->tetfield().is_set()) {
      printf("Shape field required for this function.\n"); 
      break; 
    }
    printf("For each cusp: [deg[shape-field/cusp-field], cross-ratio invt, shape]\n"); 
    m->cusp_info().print(); 
    break;

  case PrintCuspCommensurability:
    {
      if (!m->tetfield().is_set() && get_num_cusps(m->M()) > 1) 
	cout << "numerical commensurability check\n";
      vector<int> ccc = m->cusp_commensurability_classes(); 
      print_vector(cout, ccc); 
      cout << endl; 
    }
    break;

#if 0 
  case AccurateDirichletDomain:
    verbose_Dirichlet = TRUE; 
    {
      AccWEPolyhedron* P = AccDirichlet(m->M(),acc_vertex_eps,FALSE,FALSE);
      if (!P) cout << "Dirichlet domain computation failed\n"; 
      else {
	free_Dirichlet_domain(P); 
      }
    }
    verbose_Dirichlet = FALSE;
    break; 
#endif

    // Transcendental 3-manifold invariants. 

  case AccurateVolume:
    printf("Volume is: "); 
    m->accurate_volume().print();
    break;

  case ComplexVolume:
    {
      pari vol = m->complex_volume2(); 
      printf("Complex volume: "); vol.print(); 
      pari cs = gimag(vol)/(pTWO*pPI*pPI);
      cs -= pHALF * ground(pTWO * cs);
      printf("Chern-Simons (mod 1/2): ");
      cs.print(); 
    }
    break;

  case AccurateChernSimons:
    {
      // The 0 says don`t bother making a fuss if you can`t find it. 

      int corr;
      pari cs = m->accurate_chern_simons(corr); 
      if (corr == 0) {
	printf("Chern-Simons (mod 1/24): "); 
      } else if (corr == 1) {
	printf("Chern-Simons (mod 1/2): "); 
      } else {
	printf("Chern-Simons (mod 1): "); 
      }
      cs.print(); 

      if (corr && cs_diagnostic) {
	int precision, known; 
	pari diff = integer(24) * (m->raw_chern_simons() - 
				   real(m->chern_simons(&precision,&known)));
	printf("24 * (new-CS-no-correction - old-CS) = %d + ", gfloor(diff + pHALF).int_value()); 
	(diff - gfloor(diff + pHALF)).print(pp_short); 
	printf("Non-peripheral Z/2 Homology: "); m->print_npz2_homology(); 
      }
    }
    break; 

  case BootstrapCS:
    load_standard_set(path,1);
    if (!standard_set_present) {
      printf("Problem loading standard set.\n"); break; }

    double value; 
    if (m->bootstrap_CS_value(value, 1)) // 1 says report on it. 
      printf("Chern-Simons: %.16g\n", value);
    else printf("Sorry, call to bootstrap_CS_value failed for this manifold.\n"); 
    break;

    // Diagnostic. 
  case PrintCSCorrection:
    {
      int csc = face_mismatch_cs_term(m->M(), 1); // 1 says give a report on it. 
      cout << "Total correction " << csc << endl;
    }
    break;

  case PrintNonPerZ2Homology:
    printf("Non-peripheral Z/2 Homology: "); 
    m->print_npz2_homology(); 
    break; 

#if 0 
    // Obscure diagnostic stuff. 
  case PrintHomology:
    m->print_homology_matrix(); 
    break; 
#endif
  case PrintCuspHomologyKernel:
#if 0
    printf("cusp homology relations:\n"); 
    cusp_homology_relations(m->M()).print(); 
    printf("edge homology relations:\n"); 
    edge_homology_relations(m->M()).print();
    printf("cusp homology kernel:\n"); 
#endif
    m->cusp_homology_kernel().print();
    break;

  case PrintEtaDifference:
    printf("Eta inv - eta fudge: "); 
    m->eta_difference().print(); 
    break; 

  case PrintEta:
    {
      if (!m->set_eta_info()) {
	printf("Eta invariant not available for this manifold.\n"); 
	break; 
      }

      pari eta = m->eta_invariant(); 
      printf("Eta invariant: "); 
      eta.print(); 

      int correct; /* mod 1/2 */ 
      pari cs = m->accurate_chern_simons(correct); 
      if (!correct) break;
      printf("3*eta - 2*CS: "); (integer(3) * eta - pTWO * cs).print(); 
    }
    break;

  case PrintEtaFudge:
    if (!m->set_eta_info()) {
      printf("Eta invariant not available for this manifold.\n"); 
      break; 
    }
    printf("Eta fudge: "); 
    m->eta_fudge().print(); 
    break;

#if 0
  case PrintOldEtaDifference:
    printf("Eta inv - eta fudge: "); 
    m->eta_difference(1).print(); 
    break; 
  case PrintOldEta:
    {
      if (!load_eta_fudge(m, path)) break; 

      pari eta = m->eta_invariant(); 
      printf("Eta invariant: "); 
      eta.print(); 

      int correct; /* mod 1/2 */ 
      pari cs = m->accurate_chern_simons(correct); 
      if (!correct) break;
      printf("3*eta - 2*CS: "); (integer(3) * eta - pTWO * cs).print(); 
    }
    break;
  case PrintOldEtaFudge:
    if (!load_eta_fudge(m, path)) break;
    printf("Eta fudge: "); 
    m->eta_fudge().print(); 
    break;
#endif

  case SetEta:
    m->set_eta_invariant(get_pari_expr_from_user("value for eta invariant ? "));
    break;
  case SetEtaFudge:
    m->set_eta_fudge(get_pari_expr_from_user("value for eta fudge ? "));
    break;
  case SaveEtaFudge:
    { 
      string file = choose_save_file(".", "eta_fudges"); 
      if (!file.length()) break; 

      if (!m->eta_available()) break; 
      switchout((char*)file.c_str());
      fprintf(outfile, "M:%s\n", get_triangulation_name(m->M())); 
      fprintf(outfile, "EF:");
      (ground(integer(1800) * m->eta_fudge()) / integer(100)).print(); 
      fprintf(outfile, "F:"); m->the_flattening().print();
      switchout(0);
    }
    break; 

  case PrintEtaInfo:
    m->print_eta_info(); 
    break; 

  case BootEta:
    {
      load_standard_set(path,1); /* ensure that these are present */ 
      if (!standard_set_present) {
	printf("Problem loading standard set.\n"); break; }

      Triangulation* result; 
      int which_std = get_filled_standard(m->M(), &result, 1); 

      if (which_std == -1) {
	printf("Unable to find an ancestor manifold in the standard set.\n");
	break; 
      }
	  
      alg_snap *new_m = new alg_snap(result);

      // Get the eta fudge for new_m with its standard peripheral curves. 
      if (!load_eta_fudge(new_m, path)) break; 

      new_m->compute_accurate_info(0);

      // OK, now it is safe to alter the peripheral curves. 
      new_m->standardize_curves(m->M());
 
      m->set_eta_invariant(new_m->eta_invariant()); 

      delete new_m; 

      // Finally print it out. 
      pari f = m->eta_fudge();
      printf("Eta fudge: "); 
      (ground(integer(1800)*f)/integer(1800)).print(); 
    }
    break; 
#if 0
  case ExtendEtaFudge:
    {
      if (m->isometries().size()==0) {
	printf("Manifold must be isometric to some other.\n"); 
	break; 
      }
      if (!m->is_closed()) {
	printf("Manifold must be closed.\n"); 
	break; 
      }

      it_isometry const& iso = m->isometries().front(); 
      snap* t = (snap*)iso.target(); 
      t->compute_accurate_info(0,0); // Target must be accurate. 

      pari diff; 
      if (iso.orientation_preserving()) {
	printf("%s fudge - %s fudge = ", 
	       get_triangulation_name(m->M()), 
	       get_triangulation_name(t->M()));
	diff = t->eta_difference() - m->eta_difference(); 
      } else {
	printf("%s fudge + %s fudge = ", 
	       get_triangulation_name(m->M()), 
	       get_triangulation_name(t->M()));
	diff = - (t->eta_difference() + m->eta_difference()); 
      }
      (ground(integer(1800)*diff)/integer(100)).print(); 

    }
    break; 

  case PrintEtaRatDiff:
    {
      int me, l, j;
      k = get_num_cusps(m->M()); 

      if (cusp_hom_ker.type()==1) 
	cusp_hom_ker = gtrans(m->cusp_homology_kernel()); 
      if (cusp_maps.type()==1) 
	cusp_maps = m->get_isometry()->pari_cusp_mappings(); 
      if (surgeries.type()==1)
	surgeries = rvector(2 * k); 

      /* get surgery coeffs to use .. large */ 
      if (!get_input(s,"(large) value for surgery coeffs ? ", 2)) break; 
      if (sscanf(s.c_str(), "%d %d", &me, &l)!=2) break; 
      for (j = 0; j < k; j++) {
	surgeries[j] = integer(me); 
	surgeries[j+k] = integer(l); 
      }

      cusp_hom_ker.print(); 
      cusp_maps.print();
      surgeries.print(); 

      eta_rat_diff(cusp_hom_ker, surgeries, cusp_maps).print(); 
    }
    break; 
#endif

  case PrintSignature:
    {
      pari chk = m->cusp_homology_kernel();
      printf("Cusp homology kernel: "); chk.print();
      pari sc = surgery_coeffs(m->M());
      printf("sign(Y(p,q)): "); signature_Y(chk, sc).print();
    }
    break;
	
  case CheckSignature:
    {
      pari chr = cusp_homology_relations(m->M()); 
      pari ehr = edge_homology_relations(m->M()); 
      pari sc1 = surgery_coeffs2(m->M());
      pari sc2 = get_surgery_coeffs_from_user(get_num_cusps(m->M()) - num_filled_cusps(m->M())); 

      pari sc3 = gcopy(sc1);
      int l = 0; 
      int k = get_num_cusps(m->M()); 
      for (j=0; j < k; j++) {
	if (cusp_is_complete(m->M(),j)) {
	  sc3[j] = sc2[l]; 
	  sc3[j+k] = sc2[l+length(sc2)/2];
	  l++; 
	}
      }
      printf("(p,q) = "); sc1.print(); 
      printf("(r,s) = "); sc2.print(); 
      printf("sign(Y(p,q,r,s))  = ");
      signature_Y(cusp_homology_kernel(chr,ehr),sc3).print(); 
      printf("sign(Y(p,q)(r,s)) = ");
      signature_Y(cusp_homology_kernel(chr,ehr,sc1),sc2).print(); 
    }
    break; 

    // Diagnostics. 
  case EtaCuspAdj:
    { 
      pari p; 
      p = get_pari_expr_from_user("cusp change matrix ? "); 
      cusp_term_adjustment(p[0],p[1]).print();
    }
    break; 
  case CheckEta:
    m->check_eta(); 
    break; 

    // Diagnostics.. simple calculations (not manifold related).  
  case FactorMatrix:
    {
      pari p; 
      p = get_pari_expr_from_user("matrix to factor ? "); 
      factor_sl2Z_matrix(p).print(); 
    }
    break; 
  case DedekindSum:
    {
      int p,q; 
      if (!get_input(s, "please enter two integers ? ", 2)) break; 
      if (sscanf(s.c_str(), "%d %d", &p, &q)!=2) break; 
      printf("s(%d,%d) = %d\n",p,q,dedekind(p,q)); 
    }
    break;

  case SetStackSize:
    long mb;
    if (!get_input(s, "New stack size in MB ? ")) break; 
    if (!sscanf(s.c_str(), "%ld", &mb)) break;
    pari::resize_stack(1000000*mb);
    break;

  case ShowFieldComp:
    show_field_comp = (show_field_comp ? 0 : 1); 
    break; 

  case FieldPrintStyle:
    print_fields_long = (!print_fields_long); 
    break; 

  case CCPoly:
    if (!get_input(s, "Find canonical polys 0=never, 1=trace/inv, 2=always ? ")) 
      break; 
    sscanf(s.c_str(), "%d", &cc_poly);
    if (cc_poly < 0) cc_poly = 0; 
    if (cc_poly > 2) cc_poly = 2; 
    break; 

  case SetHilbertPrintstyle:
    hilb_print = (ask("show factorization ?", hilb_print & 0x1)!=0) + 
      2 * (ask("show embeddings ?", hilb_print & 0x2)!=0); 
    break; 
  case SetCSDiagnostic:
    cs_diagnostic = ask("show diagnostic ?", cs_diagnostic); 
    break; 
  case SetAccurateVertexEpsilon:
    {
      double new_vertex_epsilon;
      if (!get_input(s,"accurate vertex epsilon ? ")) break; 
      if (sscanf(s.c_str(), "%lf", &new_vertex_epsilon) != 1) break; 

      acc_vertex_eps = new_vertex_epsilon; 
    }
    break; 
  case SetPrecision:
    long d; 
    if (!get_input(s, "precision in 32-bit words ? ")) break; 
    if (sscanf(s.c_str(), "%ld", &d) != 1) {
      printf("Precision not changed\n"); 
      break;
    } 
    if (d > 30) {
      cout << "Warning: this function has changed. Precision is now specified\n";
      cout << "as a number of 32-bit words, so 1 = 9.63... decimal digits.\n"; 
      if (!get_input(s, "precision in 32-bit words ? ")) break; 
      if (sscanf(s.c_str(), "%ld", &d) != 1) {
	printf("Precision not changed\n"); 
	break;
      } 
    }
    set_precision(d);
    printf("Precision in decimal digits: %ld\n", digits); 
    init_real_globals(); /* Recompute constants like pi.. */ 
    for (j=0; j<10; j++) { 
      if (get_A(j)) get_A(j)->update_precision(); 
    }
    break;
  case SetDegree:
    if (!get_input(s, "maximum degree of polynomials to consider ? ")) 
      break; 
    int deg; 
    if (sscanf(s.c_str(), "%d", &deg) != 1) {
      printf("maximum degree not changed\n");
      break; 
    }
    if (deg < 0) deg = 1; 
    field::maxdeg = deg; 
    break;
  case SetLLLFudge:
    if (!get_input(s, "fraction of significant digits to use in lll ? ")) 
      break; 
    double f; 
    if (sscanf(s.c_str(), "%lf", &f) != 1) {
      printf("significant digits not changed\n");
      break; 
    }
    if (f < 0.1) f = 0.1; 
    if (f > 1.0) f = 1.0; 
    field::fudge = f; 
    break;
  case SetPrintstyle:
    set_printform(def_print); 
    break; 

  case -1:
    cout << "input not recognized\n";
    break; 

  default:
    env_U::process_event(what); 
  }
}


