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
#include "algebraic.hh"
#include "get_equations.hh"
#include "kernel_extras.hh"
#include "record.hh"
#include "pari_code.hh"

#ifdef PARI_2_2_OR_LATER
#include "pari_oldnames.hh"
#endif

using std::cout;
using std::cerr;
using std::endl; 

// static say xx("algebraic.cc"); 

static int manifold_record_set = 0; 
static record manifold_record = string("M"); 

const pari_printform v_short = pari_printform(1,1,'g',6,0); 

void set_manifold_record()
{
  if (manifold_record_set) return; 

  manifold_record.add_tag("IF"); // Invariant trace field. 
  manifold_record.add_tag("TF"); // Trace field. 

  manifold_record.add_tag("SF"); // Shape field. 
  manifold_record.add_tag("GF"); // Group coeff field. 

  manifold_record.add_tag("TS"); // Tetrahedron shapes. 
  manifold_record.add_tag("FG"); // Fundamental group gens. 

//  manifold_record.add_tag("TG"); 
//  manifold_record.add_tag("IG"); 
//  manifold_record.add_tag("TIG"); 
//  manifold_record.add_tag("IIT"); 
//  manifold_record.add_tag("IIS"); 
//  manifold_record.add_tag("CS"); 

//  record* fld; 

//  fld = manifold_record.find_tag("GF"); 
//  fld->add_tag("MP"); 
//  fld->add_tag("CR"); 

//  fld = manifold_record.find_tag("SF"); 
//  fld->add_tag("MP"); 
//  fld->add_tag("CR"); 

//  fld = manifold_record.find_tag("TF"); 
//  fld->add_tag("MP"); 
//  fld->add_tag("CR"); 

//  fld = manifold_record.find_tag("IF"); 
//  fld->add_tag("MP"); 
//  fld->add_tag("CR"); 
}

void alg_snap::clear_algebraic_info() 
{
  shapef.clear(); 
  _cusp_fields.resize(0); 
  _alg_shapes = pZERO; 
}

void alg_snap::clear_shape_dependencies()
{
  clear_algebraic_info(); 
  snap::clear_shape_dependencies(); 
}

void alg_snap::triangulation_changing()
{
  clear_algebraic_info(); 
  snap::triangulation_changing(); 
}

void alg_snap::clear()
{
  clear_algebraic_info(); 
  snap::clear(); 
}

void alg_snap::print() const
{
  int i, n_shapes = get_num_tetrahedra(manifold); 
  pari values; 
  char** names; 

  printf("\nTetrahedron shapes\n"); 
  pari alg = alg_shapes(); 
  pari the_shapes = (alg.type()!=1 && alg[0].type()!=1) ? alg : acc_shapes(); 

  for (i=0; i<n_shapes; i++) {
    printf("shape(%d) = ", i); the_shapes[i].print(); 
  }

  if (is_closed()) return; 

  printf("\nCusp shapes\n"); 

  int n = get_num_cusps(manifold); 

  set_accurate_shapes(manifold, the_shapes); 
  values = get_accurate_cusp_shapes(manifold); 
  clear_accurate_shapes(manifold); 

  for (i=0; i<n; ++i) {
    if (values[i]==pZERO) continue; 
    printf("cusp(%d) = ", i);
    values[i].print(); 
  }
}

void alg_snap::print_fields(int print_fields_long) const
{
  if (shapef.is_set()) { 
    printf("Shape field\n"); 
    shapef.print(print_fields_long); 
    printf("\n"); 
  }
  if (!is_closed() && _cusp_fields.size()) { 
    int j, nc = _cusp_fields.size();
    printf("Cusp fields\n");
    for (j=0; j < nc; j++) {
      _cusp_fields[j].print(print_fields_long); 
      if (!print_fields_long) printf("\n"); 
    }
  }
}

alg_snap* alg_snap::orientable_double_cover() const
{
  if (get_orientability(manifold) != nonorientable_manifold) return 0; 

  printf("passing to orientable double cover\n"); 
    Triangulation* dc = double_cover(manifold); 

  /* update the name */ 
  char name[100]; 
  strcpy(name, get_triangulation_name(manifold)); 
  strcat(name, ".or"); 
  set_triangulation_name(dc, name);

  return new alg_snap(dc); 
}
/* FINDING FIELDS */

void alg_snap::update_precision()
{
  snap::update_precision(); 
  shapef.update_precision(); 
  int i, n=_cusp_fields.size();
  for (i=0; i<n; i++) _cusp_fields[i].update_precision(); 
}

static pari pari_cross_ratio(const pari& p, const pari& q, const pari& r, const pari& s)
{
  return ((p-r)/(p-s))/((q-r)/(q-s)); 
}

static pari polymod_minpoly(const pari& elt)
{
  pari cp = caradj0(elt, 0); 
  pari denom = ggcd(cp, deriv(cp,0));  

  return (cp * denom[lgeff(denom)-1])/denom; 
}

pari alg_snap::cusp_info() const
{
  set_accurate_shapes(manifold, alg_shapes()); 
  pari c_shapes = get_accurate_cusp_shapes(manifold); 
  clear_accurate_shapes(manifold); 

  int c, nc = length(c_shapes), k=0;
  int nf = num_filled_cusps(manifold); 
  pari info = rvector(nc - nf);
  if (nf == nc) return info; // Nothing more to do. 

  int i, j; 
  pari pqrs = rvector(4);

  // get signature of shape field 
  int r1, r2, deg; 
  r1 = shapef.p_nf()[1][0].int_value(); 
  r2 = shapef.p_nf()[1][1].int_value(); 
  deg = r1 + 2*r2; 

  Boolean is_complete; 

  for (c=0; c<nc; c++) {

    get_cusp_info(manifold, c, 0, &is_complete, 0, 0, 0, 0, 0, 0, 0, 0); 
    if (!is_complete) continue; 

    info[k] = rvector(3);

    info[k][2] = lift(c_shapes[c]); 

    info[k][0] = deg/(lgeff(polymod_minpoly(c_shapes[c]))-1);

    if (deg < 4)
      info[k][1] = pZERO; 
    else {
      // take first four embeddings of r1+1,...,r1+r2,1,...,r1,-(r1+1),...,-(r1+r2). 
      // ie. take complex embeddings first, then real, then finally complex conjugates. 
      for (i = r1+1, j=0; j<4 && i<=r1+r2; i++, j++)
	pqrs[j] = shapef.numeric_value(c_shapes[c], i); 
      for (i = 1; j<4 && i<=r1; i++, j++)
	pqrs[j] = shapef.numeric_value(c_shapes[c], i); 
      for (i = r1+1; j<4 && i<=r1+r2; i++, j++)
	pqrs[j] = shapef.numeric_value(c_shapes[c], -i); 

      info[k][1] = pari_cross_ratio(pqrs[0],pqrs[1],pqrs[2],pqrs[3]);

      // If invariant is real, get rid of imaginary part resulting 
      // from roundoff error.
      if (numerical_zero(gimag(info[k][1])))
	info[k][1] = greal(info[k][1]); 
      else if (shapef.contains(gconj(shapef.root()))) {
	// Shape field is conjugation invariant so invt. is defined 
	// only up to conjugation. Therefore conjugate if necessary so 
	// imaginary part is positive.
	if (gimag(info[k][1]) < pZERO)
	  info[k][1] = gconj(info[k][1]); 
      }
    }
    k++; 
  }
  return info;
}

bool commensurable_shapes_n(pari const& a, pari const& b)
{
  pari v=rvector(4); 
  v[0]=complex(real(1),real(0)); 
  v[1]=a;
  v[2]=b;
  v[3]=a*b;
  //  cout << "shape vector = "; v.print(); 
  pari dep=lindep(v); 
  pari p1k=integer(10000); 
  if (!numerical_zero(dep * gtrans(v))) return false; 
  // cout << "dependence = "; dep.print(); 
  int i; 
  for (i=0; i<4; i++) {
    if (dep[i] > p1k) break; 
  }
  if (i<4) {
    cout << "numerical cusp commensurability test may be unreliable\n"; 
    cout << "dependency vector: "; dep.print(); 
  }
  return true; 
}

// a and b must be polymods with same min poly. 
bool commensurable_shapes(pari const& a, pari const& b)
{
  if (a.type()!=t_POLMOD) return commensurable_shapes_n(a,b); 

  int l = lgeff(a[0])-1; // = degree(a[1]). 
  if (l < 4) return true; 

  pari M, R;
  int pl, col; 

  // Check if 1,a,b,a*b
  // are linearly dependent over Q.
  M = matrix(l,4);
  M[0][0] = pONE; 

  R = gtovec(lift(a));
  pl = length(R); if (pl>l) pl=l; // just in case. 
  for (col=0; col<pl; col++) M[col][1] = R[pl-col-1];

  R = gtovec(lift(b));
  pl = length(R); if (pl>l) pl=l; // just in case. 
  for (col=0; col<pl; col++) M[col][2] = R[pl-col-1];

  R = gtovec(lift(a*b));
  pl = length(R); if (pl>l) pl=l; // just in case. 
  for (col=0; col<pl; col++) M[col][3] = R[pl-col-1];

  return rank(M.P()) < 4;
}

// Returns a vector of integers of length equal to 
// the number of cusps of the manifold. The shape field
// must be available and the manifold must not have any filled
// cusps. Entries are -1 if the function is for any reason
// unable to decide on a particular cusp. Otherwise each 
// cusp is assigned an integer >= 0 denoting its commensurability
// class. 

vector<int> alg_snap::cusp_commensurability_classes()
{
  int nc = get_num_cusps(manifold); 
  vector<int> classes(nc,-1);

  // Check no filled cusps. 
  if (num_filled_cusps(manifold)!=0) return classes; 

  // Check shape field available. 
  if (!accurate_info_set()) return classes; 

  // Get the cusp shapes as elements of the shape field. 
  pari alg = alg_shapes(); 
  bool exact = (alg.type()!=1 && alg[0].type()!=1);
  pari the_shapes = exact ? alg : acc_shapes(); 
  pari rconj; 
  bool chiral = false; 
  if (exact) chiral = !shapef.contains(gconj(shapef.root()),rconj); 

  set_accurate_shapes(manifold, the_shapes); 
  pari c_shapes = get_accurate_cusp_shapes(manifold); 
  clear_accurate_shapes(manifold); 

  // Find commensurability classes. 
  int i; 
  int c=0, j;
  for (i=0; i<nc; i++) {
    if (classes[i]!=-1) 
      continue; 
    classes[i] = c; 
    for (j=i+1; j<nc; j++) {
      if (commensurable_shapes(c_shapes[i], c_shapes[j]))
	{ classes[j] = c; continue; }
      if (chiral) continue; 
      if (exact) {
	if (commensurable_shapes(c_shapes[i], gsubst(lift(c_shapes[j]),0,rconj)))
	  classes[j] = c; 
      } else {
	if (commensurable_shapes(c_shapes[i], gconj(c_shapes[j])))
	  classes[j] = c; 
      }
      if (classes[j]==c) cout << "orientation reversing cusp commensurability!\n";
    }
    c++;
  }

  return classes;
}

/*
  which is bitwise AND of
   0x4 for shape field
   0x10 for cusp fields

  report is bitwise AND of 
   0x1 for trace of the field computation
   0x2 for printout of the answers 
*/

void alg_snap::compute_fields(int cc_poly, int report, unsigned int which)
{
  int talk = report & 1; 

  /* compute accurate info if necessary */ 
  compute_accurate_info(0);

  if (which & 0x4) {

    if (talk) printf("Computing shape field\n"); 
    pari acs = gtrans(acc_shapes()); 

    shapef.generated_by(acs, _alg_shapes, cc_poly>1, talk); 

    if (report & 2) {
      if (shapef.is_set()) {
	printf("Shape field: "); shapef.print(0); printf("\n"); 
      } else {
	printf("Shape field not found.\n"); 
      }
    }
  }

  char buf[30];

  // Find cusp fields. 
  if (!is_closed() && (which & 0x10)) {

    if (talk) printf("Computing cusp fields\n"); 

    pari values = cusp_shapes(manifold); 

    // Use accurate or exact tetrahedron shapes to get cusp shapes. 
    pari tet_shapes = (alg_shapes().type()!=1) ? alg_shapes() : acc_shapes(); 
    set_accurate_shapes(manifold, tet_shapes); 
    values = get_accurate_cusp_shapes(manifold); 
    clear_accurate_shapes(manifold); 

    int i, nc = length(values); 
    _cusp_fields.resize(nc); 

    named_field f; /* An uninitialized field. */ 
    pari c_shape = rvector(1); 

    pari exact_shape; 
    for (i=0; i<nc; i++) {
      _cusp_fields[i].clear(); 
      if (values[i].type()==1) continue; 
      
      sprintf(buf, "cusp(%d)", i); 
      _cusp_fields[i].set_name(string(buf)); 
      c_shape[0] = values[i]; 
      if (shapef.is_set())
	_cusp_fields[i].generated_by(c_shape, exact_shape, cc_poly, &shapef, 0, talk);
      else
	_cusp_fields[i].generated_by(c_shape, exact_shape, cc_poly, talk);

      if ((report & 2) && _cusp_fields[i].is_set())
	printf("Cusp %d field: ", i); _cusp_fields[i].field::print(0); printf("\n"); 
    }
  }
}

#if 0
void alg_snap::save_fields(const char* filename) const
{
  pari iit, iis, tig; 

  switchout((char*)filename); 
  string name = get_filling_name(manifold);
  fprintf(outfile, "M:%s\n", name.c_str()); 

  /* invariant trace field */ 
  if (it.is_set()) {
    fprintf(outfile, "IF:"); 
    it.print(3); 
  }

  /* trace field */ 
  if (trace.is_set()) {
    fprintf(outfile, "TF:"); 
    trace.print(3); 
  }

  /* shape field */ 
  if (shapef.is_set()) {  
    fprintf(outfile, "SF:"); 
    shapef.print(3);
    fprintf(outfile, "TS:"); 
    lift(alg_shapes()).print(); 
  }

  /* group field */ 
  if (groupf.is_set()) { 
    fprintf(outfile, "GF:"); 
    groupf.print(3); 
    fprintf(outfile, "FG:"); 
    lift(exact_gens).print(); 
  }

  switchout(0); 
}
#endif

static pari read_expr(const string& s)
{
  return p_flisexpr(s.c_str()); 
}

/* which says which fields to load. 
   set bit 0 for the invariant trace field, 
   set bit 1 for the trace field,
   set bit 2 for the shape field, 
   set bit 3 for the group coefficient field. 
*/ 
#if 0
int alg_snap::load_fields(const char *search_path, const char* name, int which)
{

  /* Look for a db_file record for this manifold. */ 
  string_map m; 
  set_manifold_record(); 
  string man_name = get_filling_name(manifold); 
  int res = load_record(search_path, name, manifold_record, man_name, m);
  if (!res) {
    printf("Unable to find an entry for %s in %s.\n", man_name.c_str(), name);
    return 0; 
  }

  compute_accurate_info(0); 
  set_tracefield_gens(); 

  /* Load the fields. If two fields coincide we simply make a copy
     rather than creating the same field again from scratch. */

  pari p_it, p_tf, p_sf, p_gf, exact;
  int i; 

  if (m.find("IF")!=m.end() && !it.is_set()) {
    p_it = read_expr(m["IF"]);
    it = field(p_it);
    printf("Invariant trace field\n"); 
    it.print(0); printf("\n");

    itgens[1] = rvector(n_itgens); 
    for (i=0; i<n_itgens; i++) {
      if (!it.contains(itgens[0][i],exact)) {
	warn("invariant trace field appears not to contain all of its generators.\n");
	itgens[1] = pZERO;
	it.clear(); 
	break; 
      }
      itgens[1][i] = exact; 
    }
  }
  if (m.find("TF")!=m.end() && !trace.is_set()) {
    p_tf = read_expr(m["TF"]);
    trace = (p_tf==p_it) ? it : field(p_tf);
    printf("Trace field\n"); 
    trace.print(0); printf("\n"); 

    tracegens[1] = rvector(n_itgens); 
    for (i=0; i<n_itgens; i++) {
      if (!trace.contains(tracegens[0][i],exact)) {
	warn("trace field appears not to contain all of its generators.\n");
	tracegens[1] = pZERO;
	trace.clear(); 
	break; 
      }
      tracegens[1][i] = exact; 
    }
  }

  if (m.find("SF")!=m.end() && !shapef.is_set()) {
    p_sf = read_expr(m["SF"]);
    shapef = (p_sf==p_it) ? it : field(p_sf);
    printf("Shape field\n"); 
    shapef.print(0); printf("\n"); 

    if (m.find("TS")!=m.end()) 
      _alg_shapes = gmodulcp(read_expr(m["TS"]),shapef.min_poly()); 
  }
  if (m.find("GF")!=m.end() && !groupf.is_set()) {
    p_gf = read_expr(m["GF"]);
    groupf = (p_gf==p_tf) ? trace : field(p_gf);
    printf("Group coefficient field\n"); 
    groupf.print(0); printf("\n"); 

    if (m.find("FG")!=m.end())
      exact_gens = gmodulcp(read_expr(m["FG"]), groupf.min_poly()); 
  }

  return 1; 
}
#endif

pari bloch_wigner_dilog(pari x)
{
  return gimag(dilog(x)) + glog(gabs(x)) * garg(pONE-x);
}

pari alg_snap::borel_regulator(const field& it, const field* e) const
{
  int s1 = (shapef.p_nf())[1][0].int_value(); 
  int s2 = (shapef.p_nf())[1][1].int_value();
  int r1 = (it.p_nf())[1][0].int_value(); 
  int r2 = (it.p_nf())[1][1].int_value();

  int e1, e2; 

  if (e) {
    e1 = (e->p_nf())[1][0].int_value();
    e2 = (e->p_nf())[1][1].int_value();
  } else {
    e1 = r1; 
    e2 = r2; 
  }

  pari reg = rvector(e2);

  pari itgen_in_sf; 
  if (!shapef.contains(it.root(), itgen_in_sf)) {
    warn("invt. trace field not contained in shape field in borel_regulator().\n"); 
    return reg; // All zeros. 
  }
  itgen_in_sf = lift(itgen_in_sf); 

  pari itgen_in_e; 
  if (e && !e->contains(it.root(), itgen_in_e)) {
    warn("invt. trace field not contained in ext field in borel_regulator().\n"); 
    return reg; // All zeros. 
  }
  if (e) itgen_in_e = lift(itgen_in_e); 

  // For each place of the invariant trace field, find a corresponding 
  // place of the shape field. 
  vector<int> sfplace(r1+r2+1);

  int i, ind, j, n;

  for (i=0; i<r1+r2+1; i++) sfplace[i] = 0; 

  pari rt; 
  for (i = s1+1; i <= s1+s2; i++) {
    rt = gsubst(itgen_in_sf, 0, shapef.root(i)); 
    ind = it.root_number(rt); 

    if (ind < 0) {
      j = -i; ind = -ind; 
    } else j = i; 

    if (ind <= r1 || ind > r1+r2) continue; 
    if (!sfplace[ind]) sfplace[ind] = j;
  }

  pari sh; 
  int n_shapes = get_num_tetrahedra(manifold); 

  for (i = e1+1; i < e1+e2+1; i++) {

    if (e) {
      rt = gsubst(itgen_in_e, 0, e->root(i)); 
      ind = it.root_number(rt); 

      if (ind < 0) {
	n = -sfplace[-ind]; 
      } else {
	n = sfplace[ind];
      }

    } else {
      n = sfplace[i]; 
    }

    if (n==0) {
      // This will happen if this particular place of e happens to 
      // induce a real place of the inv. trace field. The corresponding 
      // entry in the borel regulator is then simply 0. 
      continue; 
    }

    for (j = 0; j < n_shapes; j++) {
      sh = shapef.numeric_value(alg_shapes()[j], n);
      if (gimag(sh) >= pZERO)
	reg[i-e1-1] += bloch_wigner_dilog(sh); 
      else 
	reg[i-e1-1] -= bloch_wigner_dilog(gconj(sh));
    }
  }

  return reg; 
}

static void print_row(const pari& row, int n_tet, const pari& shapes, const pari& log_res, int& OK)
{
  pari alg_res = monomial(shapes,row)[1];

  int j; 
  for (j=0; j<n_tet; j++) {
    printf("%d", row[j].int_value());
    if (j<n_tet-1) printf(", "); 
    else printf("; ");
  }
  for (j=0; j<n_tet; j++) {
    printf("%d", row[j+n_tet].int_value());
    if (j<n_tet-1) printf(", "); 
    else printf("; ");
  }
  printf("%d -> ", row[2*n_tet].int_value());

  alg_res.print(no_newline); 
  printf(" : ");
  abs(log_res).print(v_short); 

  if (alg_res != pONE || !small(log_res)) 
    OK = 0; 
}

static void check_row(const pari& row, const pari& shapes, const pari& log_res, int& OK)
{
  pari alg_res = monomial(shapes,row)[1];
  if (alg_res != pONE || !small(log_res)) OK = 0; 
}

// -2 - shape field unknown
// -1 - alg shapes do not agree with numeric shapes
//  0 - verification failed
//  1 - verification OK
//  2 - gluing OK but some flat simplices
//  3 - gluing OK but non-geometric. 

int alg_snap::verify() const
{
  int OK = 1; 

  /* verify the exact hyperbolic structure and print a full report */

  if (!shapef.is_set()) {
    return -2; 
  }

  int i, n = length(alg_shapes());

  pari error; 
  for (i=0; i<n; i++) {
    error = abs(shapef.numeric_value(alg_shapes()[i]) - acc_shapes()[i]);
    if (!numerical_zero(error)) return -1; 
  }

  int c = get_num_cusps(manifold);
  pari eqns,alg_res,log_res;
  if (num_filled_cusps(manifold)==0) {
    eqns = get_complete_equations(manifold); 
    log_res = gtrans(eqns) * concat(accurate_log_shapes(), pPI*pI);
    for (i=0; i<c; i++)
      check_row(eqns[i], alg_shapes(), log_res[i],OK);
    for (i=0; i<c; i++)
      check_row(eqns[i+c], alg_shapes(), log_res[i+c],OK);
    for (i=0; i<n; i++)
      check_row(eqns[i+2*c], alg_shapes(), log_res[i+2*c],OK);
  } else {
    eqns = get_filled_equations(manifold); 
    log_res = gtrans(eqns) * concat(accurate_log_shapes(), pPI*pI);
    for (i=0; i<c; i++)
      check_row(eqns[i], alg_shapes(), log_res[i],OK);
    for (i=0; i<n; i++)
      check_row(eqns[i+c], alg_shapes(), log_res[i+c],OK);
  }

  if (!OK) return 0;

  pari acc_log_shapes = accurate_log_shapes();
  pari angle; 
  bool flat=false, nongeometric=false; 
  for (i=0; i<n; i++) {
    angle = gimag(acc_log_shapes[i]);
    if (numerical_zero(angle) || numerical_zero(angle-pPI)) {
      flat=true;
    } else {
      if (angle < pZERO || angle > pPI) nongeometric=true; 
    }
  }

  if (nongeometric) return 3; // verified but non-geometric. 
  if (flat) return 2; // verified and (at least partially) flat
  return 1; // verified and geometric
}

int alg_snap::verify_and_report() const
{
  int OK = 1; 

  /* verify the exact hyperbolic structure and print a full report */

  if (!shapef.is_set()) {
    printf("Shape field must be computed before structure can be verified\n");
    return 0; 
  }

  printf("Shapes (Numeric)\n"); 
  int i, n = length(alg_shapes());
  for (i=0; i<n; i++) {
    printf("shape(%d) = ", i);
    acc_shapes()[i].print();
  }

  printf("\nShape Field\n"); 
  printf("min poly: ");
  shapef.min_poly().print(); 
  printf("root number: %d\n", shapef.root_number());  
  printf("root: "); shapef.root().print();

  printf("\nShapes (Exact)\n");
  pari error; 
  for (i=0; i<n; i++) {
    printf("shape(%d) = ", i);
    alg_shapes()[i][1].print(no_newline);
    printf(" "); 
    error = abs(shapef.numeric_value(alg_shapes()[i]) - acc_shapes()[i]);
    error.print(v_short); 
  }

  /* where are alg_shapes set?? */

  printf("\nGluing Equations\n"); 

  int c = get_num_cusps(manifold);
  pari eqns,alg_res,log_res;
  if (num_filled_cusps(manifold)==0) {
    eqns = get_complete_equations(manifold); 
    log_res = gtrans(eqns) * concat(accurate_log_shapes(), pPI*pI);
    printf("Meridians:\n");
    for (i=0; i<c; i++)
      print_row(eqns[i], n, alg_shapes(), log_res[i],OK);
    printf("Longitudes:\n");
    for (i=0; i<c; i++)
      print_row(eqns[i+c], n, alg_shapes(), log_res[i+c],OK);
    printf("Edges:\n");
    for (i=0; i<n; i++)
      print_row(eqns[i+2*c], n, alg_shapes(), log_res[i+2*c],OK);
  } else {
    eqns = get_filled_equations(manifold); 
    log_res = gtrans(eqns) * concat(accurate_log_shapes(), pPI*pI);
    printf("Filling curves:\n");
    for (i=0; i<c; i++)
      print_row(eqns[i], n, alg_shapes(), log_res[i],OK);
    printf("Edges:\n");
    for (i=0; i<n; i++)
      print_row(eqns[i+c], n, alg_shapes(), log_res[i+c],OK);
  }

  if (OK) 
    printf("Verification OK\n"); 
  else 
    printf("Verification failed!\n"); 

  return OK; 
}

