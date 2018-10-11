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
#include "snap.hh"
#include "get_equations.hh"
#include "kernel_extras.hh"
#include "snappea/unix_io.h"
#include "snappea/kernel.h"
#include "pariwrap.hh"
#include "pari_code.hh"
#include "printable.hh"

#ifdef PARI_2_2_OR_LATER
#include "pari_oldnames.hh"
#endif

extern "C" {
#include "pari_ext.h"
}

using std::cout;
using std::endl;
using std::setw;

// static prt_init xx("snap.cc\n"); 
// static say xx("snap.cc"); 

#define parent i_triangulation

/* BASIC STUFF */

snap::snap()
: _eta_invariant_known(false), _eta_fudge_known(false)
{}

snap::snap(Triangulation *m)
: _eta_invariant_known(false), _eta_fudge_known(false)
{
  set_manifold(m); 
}

snap::~snap()
{}

string snap::filling_name() const
{
  return get_filling_name(manifold); 
}

int snap::is_closed() const
{ 
  return get_num_cusps(manifold) == num_filled_cusps(manifold); 
}

int snap::is_complete() const
{ 
  return num_filled_cusps(manifold)==0; 
}

void snap::clear_shape_dependencies()
{
  // doing surgery so save the fudge. 
  if (eta_available()) { 
    compute_accurate_info(0); 
    if (flattening.type()==t_INT) find_a_flattening(); 
    set_eta_fudge(eta_fudge()); 
  }

  _acc_shapes = pZERO; 
  _eta_invariant_known = false; 
  parent::clear_shape_dependencies(); 
}

static bool cone_surgery(Triangulation *manifold)
{
  int i, n_cusps = get_num_cusps(manifold);
  pari sc = surgery_coeffs(manifold);
  for (i=0; i<n_cusps; i++) {
    if (cusp_is_complete(manifold,i)) continue; 
    if (sc[i+2*n_cusps]!=pONE) return true; 
  }
  return false; 
}

void snap::triangulation_changing()
{
  if (cone_surgery(manifold)) {
    _eta_invariant_known = false; 

  } else if (eta_available()) { 

    // keeping manifold so save the invariant. 
    compute_accurate_info(0); 
    set_eta_invariant(eta_invariant()); 
  }

  _acc_shapes = pZERO; 
  _eta_fudge_known = false; 
  flattening = pZERO; 
  parent::triangulation_changing(); 
}

void snap::clear()
{
  _acc_shapes = pZERO; 
  _eta_invariant_known = false; 
  _eta_fudge_known = false; 
  flattening = pZERO; 
  parent::clear(); 
}

/* STUFF TO DO WITH READING AND MODIFYING NUMERIC MANIFOLD STRUCTURES */

void snap::change_peripheral_curves(MatrixInt22 change_matrices[])
{
  int i, n = get_num_cusps(manifold); 

  if (!eta_available()) {
    parent::change_peripheral_curves(change_matrices); 
    return; 
  }

  compute_accurate_info(0);
  if (flattening.type()==t_INT) find_a_flattening(); 

  if (is_closed()) {
    /* eta invariant will be unchanged if manifold is closed, therefore
       force computation of the eta invariant and discard the fudge */ 
    set_eta_invariant(eta_invariant());

  } else {
    
    /* if the manifold has cusps we change the fudge factor instead */ 
    pari fudge = eta_fudge(); 
      
    /* get surgery coeffs replacing cusps with "arbitrary" fillings */
    pari sc0 = ::surgery_coeffs(manifold); 
    pari sc = cvector(2*n); 
    for (i=0; i<n; i++) {
      if (sc0[i+2*n]==pZERO) {
	sc[i] = integer(2); 
	sc[i+n] = integer(5); 
      } else {
	if (sc0[i+2*n]!=pONE) {
	  warn("Don't know how to deal with eta invariant at this point\n");
	  _eta_invariant_known = _eta_fudge_known = false; 
	  parent::change_peripheral_curves(change_matrices); 
	  return;
	}
	sc[i] = sc0[i];
	sc[i+n] = sc0[i+n];
      }
    }
    
    /* get change matrices into a slightly different form */ 
    pari chmx = matrix(2*n,2*n); 
    for (i=0; i<n; i++) {
      chmx[i][i] = change_matrices[i][1][1]; 
      chmx[i+n][i] = -change_matrices[i][1][0];
      chmx[i][i+n] = -change_matrices[i][0][1];
      chmx[i+n][i+n] = change_matrices[i][0][0];
    }
    
    /* get change to signature terms resulting from change_matrices */ 
    pari k0 = cusp_homology_kernel(); 
    fudge += signature_Y(chmx * k0, chmx * sc) - signature_Y(k0, sc);  

    /* get change to dedekind sum terms resulting from change_matrices */ 
    pari pq = cvector(2); 
    pari chmx2 = matrix(2,2); 
    for (i=0; i<n; i++) {
      pq[0] = sc[i]; 
      pq[1] = sc[i+n]; 
      chmx2[0][0] = change_matrices[i][1][1]; 
      chmx2[1][0] = -change_matrices[i][1][0];
      chmx2[0][1] = -change_matrices[i][0][1];
      chmx2[1][1] = change_matrices[i][0][0];
      fudge += cusp_term_adjustment(chmx2, pq); 
    }
    set_eta_fudge(fudge); 
  }

  parent::change_peripheral_curves(change_matrices); 
}

#if 0
SolutionType snap::do_Dehn_surgery(const vector<double>& coeffs)
{
  return parent::do_Dehn_surgery(coeffs); 
}

SolutionType snap::find_complete_hyperbolic_structure()
{
  return parent::find_complete_hyperbolic_structure(); 
}
#endif

void snap::reorient()
{
  parent::reorient(); 
  if (eta_available()) _eta_invariant = -_eta_invariant;
  find_a_flattening(); 
}

bool snap::canonize(double area_multipliers[]) 
{
  int res = parent::canonize(area_multipliers); 
  find_a_flattening(); 
  return res; 
}

void snap::simplify()
{
  parent::simplify(); 
  find_a_flattening(); 
}

void snap::randomize()
{
  parent::randomize(); 
  find_a_flattening();
}

void snap::fill()
{
  parent::fill();
  find_a_flattening(); 
}

void snap::test_bc_subdivision()
{
  compute_accurate_info(0); 
  set_accurate_shapes(manifold, _acc_shapes); 
  Triangulation* bc = barycentric_subdivision(manifold); 
  clear_accurate_shapes(manifold); 
  free_triangulation(bc); 
}

/* PRINTING, LOADING AND SAVING MANIFOLD INFO */

void snap::print() const
{
  int i, n = get_num_tetrahedra(manifold); 
  pari values; 
  char** names; 

  printf("\nTetrahedron shapes\n"); 
  for (i=0; i<n; i++) {
    printf("shape(%d) = ", i); 
    _acc_shapes[i].print(); 
  }

  printf("\nCusp shapes\n"); 
  values = cusp_shapes(manifold); 
  n = length(values); 
  for (i=0; i<n; i++) {
    if (values[i]==pZERO) continue; 
    printf("cusp(%d) = ", i); 
    values[i].print(); 
  }
}


/* STUFF TO DO WITH 3-MANIFOLD INVARIANTS */

void snap::print_homology_matrix() const
{
  int n_cusps = get_num_cusps(manifold);
  int n_shapes = get_num_tetrahedra(manifold); 
  pari mx = cusp_homology_relations(manifold); 
  pari mxe = edge_homology_relations(manifold); 
  pari sc = ::surgery_coeffs(manifold); 

  /* now print it all out */ 
  printf("number of generators: %d\n", manifold->num_generators); 
  printf("meridians\n");
  int i; 
  for (i = 0; i < n_cusps; i++) 
    gtrans(mx[i]).print(); 
  printf("longitudes\n");
  for (; i < 2 * n_cusps; i++) 
    gtrans(mx[i]).print(); 
  printf("filled homology relations\n");
  for (i=0; i < n_cusps; i++) {
    if (sc[i+2*n_cusps]!=pZERO)
      gtrans((sc[i] * mx[i] + sc[i+n_cusps] * mx[i+n_cusps])/sc[i+2*n_cusps]).print(); 
  }
  printf("edge homology relations\n");
  for (i=0; i < n_shapes; i++) 
    gtrans(mxe[i]).print();
}

void snap::print_npz2_homology() const
{
  AbelianGroup *ab = non_per_z2_homology(manifold);
  if (ab->num_torsion_coefficients == 0) 
    printf("0\n"); 
  else {
    printf("Z/%ld", ab->torsion_coefficients[0]);
    int i=0;
    while (++i < ab->num_torsion_coefficients)
      printf(" + Z/%ld", ab->torsion_coefficients[i]);
    printf("\n"); 
  }
  free_abelian_group(ab); 
}

/* COMPUTING ACCURATE MANIFOLD INFO */

/* each column of m gives a list of the powers of 
   first z[i] then (1-z[i]). there should be as many
   columns as there are variables, z[1],...,z[n] in
   order for the derivative matrix to be invertible 
   in newton's method. */ 
   
static int newton_step(pari& z, pari m, pari const& eps, int report)
{
  int i, j, n = length(z), ddig; 
  pari ersz; 

  /* pari pz = cvector(n); */
  pari er = cvector(n); 
  pari d = matrix(n,n); 

  /* first double the precision of z */ 
  ddig = (int)((double)(length(abs(z[1]))-2) * pariK * 2.0);

  if (ddig < 16) ddig = 16; 
  if (ddig > digits + 9) ddig = digits + 9; 
  pari Z = gprec(z,ddig); 

#if 0
  /* check to avoid division by zero */ 
  pari zeps = real(1.e-10); 
  for (j=0; j<n; j++) {
    if (gcmp(gabs(z[j]),zeps)<0 || 
	gcmp(gabs(z[j]-pONE),zeps)<0) 
      return true; // stop iteration. 
  }
#endif

  /* compute the error and derivative */ 
  for (i=0; i<n; i++) {
    er[i] = monomial(Z,m[i]) - pONE; 
    for (j=0; j<n; j++) {
      d[j][i] = (m[i][j]/z[j] - m[i][j+n]/(pONE - z[j]));
    }
  }

  /* check how close we are */ 
  pari adj = gauss(d, er);
  ersz = vecmax(abs(adj)); 
  if (report) ersz.print();

  z = Z - adj; 

  return (ersz < eps);
}


/* polish = 0 .. get accurate shapes if none previously computed.
   polish = 1 .. try to improve the accuracy of any existing shapes. 
   polish = 2 .. get approximate shapes direct from the manifold. 
*/ 

void snap::compute_accurate_info(int polish, int report)
{
  /* do nothing if we already have it and 
     we're not being asked to polish it */ 
  if (accurate_info_set() && polish!=1) return; 

  /* accurate_shapes.type()!=1 means that accurate shapes
     has already been set and we're just polishing it, eg. 
     after an increase in precision */ 
  if (!accurate_info_set()) {
    _acc_shapes = tet_shapes(manifold); 
  }

  /* don't try to polish bad solutions */
  int sol = get_filled_solution_type(manifold);
  if (sol==0 || sol > 3) return; 

  pari eps = pow(real(10),-digits); 

  /* compute shapes as accurately as we can */ 
  if (polish!=2) {
    if (report) 
      printf("Polishing hyperbolic structure to %ld decimal places...\n", digits);
    int i = 0, done = 0; 
    pari eqns = get_filled_equations(manifold); 
    pari lisubmatrix = indexrank(eqns); 
    eqns = extract(eqns, lisubmatrix[1]); 
    
    while (!done && ++i < 20) {
      done = newton_step(_acc_shapes, eqns, eps, report); 
    }
    if (done) // do it one more time for luck! 
      newton_step(_acc_shapes, eqns, eps, report);
  }

  /* Force (re)computation of the fundamental group. */
  clear_groups(); 
}

void snap::update_precision()
{
  if (accurate_info_set()) compute_accurate_info(1);
  // discard any inaccurate eta_invariant.
  if (eta_available()) set_eta_fudge(eta_fudge()); 
}



GroupPresentation* snap::group(int s) const
{
  GroupPresentation* result;
  if (accurate_info_set()) {
    set_accurate_shapes(manifold, _acc_shapes); 
    result = i_triangulation::group(s);
    clear_accurate_shapes(manifold); 
  } else {
    result = i_triangulation::group(s); 
  }
  return result; 
}

pari monomial(pari z, pari p)
{
  int i, n = length(z); 

  pari prod = pow(-pONE, p[2*n]); 

  for (i=0; i<n; i++) {
    prod *= (pow(z[i], p[i]) * pow(pONE - z[i], p[i+n])); 
  }
  return prod; 
}

/* Returns log(z1),...,log(zn),log(1-z1),...,log(1-zn) */ 
pari snap::accurate_log_shapes() const
{
  int n_shapes = get_num_tetrahedra(manifold); 
  pari approx_log_shapes = ::log_shapes(manifold); 
  pari acc_log_shapes = cvector(2 * n_shapes); 
  pari two_pi = pTWO * pPI;
  int i; 
  long e; 
  for (i=0; i<n_shapes; i++) {
    acc_log_shapes[i] = glog(_acc_shapes[i]); 
    acc_log_shapes[i] += two_pi * pI * 
      grndtoi(gimag(approx_log_shapes[i] - acc_log_shapes[i]) / two_pi, &e);
    if (e > -10) printf("error bits %ld\n", e); 
    acc_log_shapes[i+n_shapes] = glog(pONE - _acc_shapes[i]); 
    acc_log_shapes[i+n_shapes] += two_pi * pI * 
      grndtoi(gimag(approx_log_shapes[i+n_shapes] - acc_log_shapes[i+n_shapes]) / two_pi, &e);
    if (e > -10) printf("error bits %ld\n", e); 
  }
  return acc_log_shapes; 
}

pari snap::acc_holonomies() const
{
  pari log_shapes = accurate_log_shapes(); 
  log_shapes = concat(log_shapes, pPI * pI); 
  return cusp_equations(manifold) * log_shapes; 
}

pari snap::core_geodesics() const
{
  int i, n_cusps = get_num_cusps(manifold); 
  pari cg = rvector(n_cusps); 
  pari sc = ::surgery_coeffs(manifold); 
  pari h = acc_holonomies(); 
  pari p, q, qq, n, d; 
  
  /* Complex length of the core geodesic of a (p,q)-surgered 
     cusp with logarithmic holonomies u, v is given by:

     v/p + 2 pi i q'/p  when p != 0 and with q' satisfying

     qq' = -1 mod(p) and 0 <= q'/p < 1. This is the formula 
     used when neither p nor q are zero. If q = 0 I think we 
     can still use this formula with q' = 0. If p = 0 then
     we use the alternate formula: 

     2 pi i p'/q - u/q 

     Here is where I'm not sure: I think that in this case
     generally we have pp' = 1 mod(q) and 0 < p'/q <= 1.  
     My guess is that when p = 0 we want p'/q = 1. 

     In what follows u = cusp->holonomy[ultimate][M] and 
     v = cusp->holonomy[ultimate][L]. 
   */

  for (i=0; i<n_cusps; i++) {
    d = sc[i+2*n_cusps];
    if (d==pZERO) continue;
    p = sc[i]; 
    q = sc[i+n_cusps]; 
    n = ggcd(p,q);
    p /= n;
    q /= n;

    if (p != pZERO) {
      if (q == pZERO) qq = pZERO;
      else qq = gmod(-pONE/q, p);
      cg[i] = d*(qq * pTWO_PI * pI + h[i+n_cusps])/(p*n); 
    } else { /* p == 0 */ 
      cg[i] = d*(pTWO_PI * pI - h[i])/(q*n); 
    }
  }
  return cg; 
}


/*
 * dilog(z) goes down by 2pi i log(z) on branch cut (1,infinity).
 * dilog(z) is analytic elsewhere.
 * 
 * log(z) goes down by 2pi i on branch cut (-infinity,0)
 * log(z) is analytic elsewhere.
 * 
 * log(1-z) goes up by 2pi i on branch cut (1, infinity)
 * log(1-z) is analytic elsewhere.
 *
 * rogers(z) := dilog(z) + (1/2) log(z) log(1-z)
 *
 * Principal branch..
 * rogers(z) goes down by pi i log(z)   on branch cut (1,infinity)
 * rogers(z) goes down by pi i log(1-z) on branch cut (-infinity,0)
 *
 * Another branch might be..
 * rogers(z;p,q) := rogers(z) + p pi i log(z) + q pi i log(1-z)
 *
 * rogers(z;p,q) goes down by pi i log(z)   + 2 pi^2 q on cut (1,infinity)
 * rogers(z;p,q) goes down by pi i log(1-z) - 2 pi^2 p on cut (-infinity,0)
 * 
 * The branch in full is..
 * rogers(z;p,q,r) := rogers(z) + p pi i log(z) + q pi i log(1-z) + 2 r pi^2
 * 
 * To get analytic continuation.. 
 * crossing (1,infinity):  p goes up by 1, r goes up   by q
 * crossing (-infinity,0): q goes up by 1, r goes down by p
 */

/*
 * A little structure to keep track of which branch of the
 * rogers dilogarithm we are at. 
 */

struct sh_branch {
  int p, q, r, upper;

  sh_branch() : p(0),q(0),r(0),upper(1) {}

  void update(int wide_angle); 
};

void sh_branch::update(int wide_angle)
{
  switch (wide_angle) {
  case 0: /* crossing [-inf,0] */ 
    q += upper; 
    r -= upper * p;
    break;
  case 1: /* crossing [1,+inf] */ 
    p += upper;
    r += upper * q; 
    break; 
  default:
    break;
  }
  upper = -upper; 
}

/*
 * ShapeInversion* history is a stack, the most recent
 * shape inversion being at the top. To find the branch 
 * we have to read it from the bottom. This is most
 * easily done by recursion. 
 */

static sh_branch get_branch(ShapeInversion *history)
{
  if (!history) return sh_branch();

  sh_branch branch = get_branch(history->next);
  branch.update(history->wide_angle);
  return branch;
}

static pari rogers_dilog_with_history(pari z, ShapeInversion *history)
{
  sh_branch branch = get_branch(history);

  /* 
   * Don't end up on the wrong branch when imaginary part is small.
   * This can only happen when the accurate shape z and SnapPea's
   * shape have opposite signs for their imaginary parts. 
   */ 
  if (numerical_zero(gimag(z)) && 
      branch.upper != gsigne(gimag(z))) z = gconj(z); 

  return dilog(z) + pHALF * log(z) * log(pONE-z) + 
    integer(branch.p) * pPI * pI * log(z) + 
    integer(branch.q) * pPI * pI * log(pONE-z) + 
    integer(branch.r) * pTWO_PI * pPI; 
}

void snap::find_a_flattening()
{
  choose_generators(manifold, TRUE, FALSE);
  // Get flattening with parity condition, for the complete structure. 
  flattening = compute_flattening(manifold, true, false); 
}

pari snap::dilog_sum() const 
{
  /* make sure we have a flattening */ 
  if (flattening.type()==t_INT) {
    ((snap*)(this))->find_a_flattening(); 
  }

  /* get the shapes */ 
  pari Z0 = accurate_log_shapes(); 
  int i, n_shapes = get_num_tetrahedra(manifold); 
  pari d_sum; 
  ShapeInversion **sh = shape_histories(manifold); 

  pari pi_i_over_2 = pPI_OVER_2 * pI;

  for (i = 0; i < n_shapes; i++) {
    d_sum += rogers_dilog_with_history(_acc_shapes[i],sh[i]) - 
      pi_i_over_2 * (flattening[i]          * Z0[i+n_shapes] - 
		     flattening[i+n_shapes] * Z0[i]);
  }

  delete [] sh; 
  return d_sum; 
}

pari snap::complex_volume() const
{
  if (!accurate_info_set()) {
    warn("complex_volume() requires accurate shape info\n");
    return pZERO; 
  }

  /* get the sum of the complex lengths of all core geodesics */
  int n_cusps = get_num_cusps(manifold); 
  pari core_length = core_geodesics(); 
  pari core_len_sum; 
  for (int i = 0; i < n_cusps; i++)
    core_len_sum += core_length[i];

  return (-pPI_OVER_2 * core_len_sum - pI * dilog_sum()); 
}

pari snap::accurate_volume() const
{
  return greal(complex_volume()); 
}

/* To calculate the Chern-Simons invariant mod 1, we use 
   the Atiyah-Patodi-Singer result:
   
   2 CS(M) = 3 \eta(M) + (number of 2-primary summands in H_1(M;Z))  mod 2

   for closed, oriented Riemannian 3-manifolds. 

   The function sets correct to 0 if CS is only known mod 1/24, 
   to 1 if CS is known mod 1/2 and to 2 if CS is known mod 1. 
*/

pari snap::accurate_chern_simons(int& correct) const
{
  pari cs; 
  AbelianGroup* hg;

  if (eta_available() && (hg = homology(manifold))) {

    // Count the number of 2-primary summands in H_1(M;Z). 
    int i, n = hg->num_torsion_coefficients, n2p_summs = 0;
    for (i=0; i<n; i++) {
      if (hg->torsion_coefficients[i] && !(hg->torsion_coefficients[i] & 0x1))
	n2p_summs++; 
    }
    free_abelian_group(hg); 

    cs = pHALF * (integer(3) * eta_invariant() + integer(n2p_summs)); 
    bool cusped = get_num_cusps(manifold)!=num_filled_cusps(manifold);

    // Normalize appropriately: CS is only defined mod 1/2 for
    // cusped manifolds. Of course this is slightly inefficient
    // since it makes the homology computation irrelevant.
    if (cusped) {
      cs -= pHALF * ground(pTWO * cs);
      correct = 1;
    } else {
      cs -= ground(cs);
      correct = 2;
    }

    return cs; 
  }

  cs = gimag(complex_volume())/(pTWO * pPI * pPI); 

  // The line above (conjecturally??) gives CS correct mod 1/24, and 
  // to very high precision. Jeff Weeks` code called below, computes 
  // the correct value to ordinary C-double precision. It is then easy
  // to figure out what multiple of 1/24 should be added to the above
  // CS value to get one which is correct AND very precise. 

  int precision; 
  Boolean is_known, requires_init;
  double value; 
  get_CS_value(manifold, &is_known, &value, &precision, &requires_init); 
  
  if (!is_known) value = 0.;  
  cs += ground((real(value) - cs) * integer(24))/integer(24); 

  correct = is_known; 
  return cs;
}



static pari rogers_dilog(pari const& z, pari const& p, pari const& q)
{
  return pHALF * log(z) * log(pONE-z) + dilog(z) - 
    pPI*pI/pTWO * (p * log(pONE - z) - q * log(z)) -
    pPI*pPI/real(6); 
}

static pari rogers_dilog_sum(Triangulation* manifold, pari const& flattening, bool with_signs)
{
  int i, n = get_num_tetrahedra(manifold); 
  pari sum, z;

  Tetrahedron* tet; 
  for (tet = manifold->tet_list_begin.next; 
       tet!= &manifold->tet_list_end; 
       tet = tet->next) {
    i = tet->index;
    z = tet->acc_shape; 
    if ((!with_signs) || tet->has_correct_orientation) 
      sum += rogers_dilog(z, flattening[i], flattening[i+n]); 
    else 
      sum -= rogers_dilog(z, flattening[i], flattening[i+n]); 
  }
  return sum; 
}

pari snap::complex_volume2() const
{
  set_accurate_shapes(manifold, _acc_shapes); 

  Triangulation* bc = barycentric_subdivision(manifold); 

  pari flat = compute_flattening(bc, true, true, true); 

  pari vol = rogers_dilog_sum(bc, flat, true)/pI; // with signs. 

  clear_accurate_shapes(manifold); 
  free_triangulation(bc); 
  return vol; 
}

// There should be a better place for this, right? 

int_matrix to_int_matrix(pari const& m)
{
  int_matrix M; 
  int r, c, i, j; 

  switch (m.type()) {
  case t_VEC:
  case t_COL:
    r = 1; 
    c = m.length()-1; 
    M.set_dims(1,c);
    for (i=0; i<c; ++i) 
      M[0][i] = m[i].int_value();
    break;
  case t_MAT:
    r = m[0].length()-1; 
    c = m.length()-1; 
    M.set_dims(r,c);
    for (i=0; i<r; i++)
      for (j=0; j<c; j++)
	M[i][j]=m[j][i].int_value();
    break;
  default:
    break;
  }

  return M;
}

void show_matprod(pari const& mat, pari const& v)
{
  int_matrix eq = to_int_matrix(mat); 
  pari r = mat * v; 
  int i, n = eq.rows, w = field_width(eq);
  bool isint = (r[0].type()==t_INT); 
  for (i=0; i<n; i++) {
    pretty_print(eq, w, i);
    if (isint) 
      cout << " -> " << r[i].int_value() << endl; 
    else 
      cout << " -> " << to_complex(r[i]) << endl; 
  }
}

static void print_perm(Permutation perm)
{
  int i; 
  for (i=0; i<4; i++) {
    cout << EVALUATE(perm, i); 
  }
}

static int PE3(Permutation perm, int e)
{
  return edge3_between_vertices[EVALUATE(perm, 0)][EVALUATE(perm, e+1)];
}

static void get_zpq(vector<Complex> const& w, Complex& z, int& p, int& q)
{
  int i, j; 
  double s[] = {-1,1};
  Complex pmz = complex_exp(w[0]);
  Complex pmz1= complex_exp(-w[1]); 

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      if (complex_small(s[i]*pmz - (One - s[j]*pmz1), 1e-5)) break; 
    }
    if (j<2) break; 
  }
  if (i==2) { cout << "problem in get_zpq!\n"; return; }

  z = s[i]*pmz;
  Complex lz = complex_log(z,0.);
  Complex lz1 = complex_log(One - z,0.);

  p = int(floor((w[0].imag - lz.imag)/PI + .5));
  q = int(floor((w[1].imag + lz1.imag)/PI + .5));

  if (!(complex_small(w[0]-lz-p*PI*I) &&
	complex_small(w[1]+lz1-q*PI*I))) {
    cout << "trouble finding p&q in get_zpq\n"; return; 
  }
}

static Complex rogers_dilog(Complex const& z, int p, int q)
{
  return to_complex(rogers_dilog(to_pari(z), pari(p), pari(q)));
}

void complex_volume3(Triangulation* T, vector<int> const& e_or)
{
  // Get flattening with parity for filled manifold 

  pari fl = concat(compute_flattening(T, true, true, false),pONE); 
  cout << "Flattening with parity\n";
  pretty_print(to_int_matrix(fl)); 

  pari eqns = gtrans(get_filled_equations(T)); 
  cout << "Filling equations:\n";
  show_matprod(eqns, fl);

  pari pmat = get_parity_matrix(T);
  cout << "Parity matrix:\n"; 
  show_matprod(pmat, fl); 
  

  cout << "Flattened log parameters\n";
  
  int i, n=get_num_tetrahedra(T); 
  vector<Complex> logsh = get_log_shapes(T); 
  vector<Complex> fsh(3*n); 
  vector<Complex> flat(3); 
  double p, q; 
  Complex lz, lz1; 
  Complex PiI = Complex(0.,PI);
  for (i=0; i<n; i++) {
    lz  = logsh[i  ];
    lz1 = logsh[i+n];
    
    p = -fl[i  ].double_value(); 
    q =  fl[i+n].double_value(); 

    fsh[3*i  ] =        lz +   p   * PiI; 
    fsh[3*i+1] = -lz1      +   q   * PiI;
    fsh[3*i+2] =  lz1 - lz - (p+q) * PiI; 

    flat[0] = fsh[3*i  ];
    flat[1] = fsh[3*i+1]; 
    flat[2] = fsh[3*i+2];
    cout << PSeq(flat) << endl;
  }

  pari fshp = cvector(3*n);
  for (i=0; i<3*n; ++i) fshp[i] = to_pari(fsh[i]);
  
  eqns = gtrans(get_filled_equations(T, false, true)); 
  cout << "Flat equation results\n"; 
  show_matprod(eqns, fshp); 
  
  cout << "Rogers dilog sum contributions\n";
  
  // Calculate complex volume with permuted vertices.
  
  Permutation perm = IDENTITY_PERMUTATION; 
  double psig = 1; 
  Complex z, rd, cv=Zero; 
  int ip, iq; 
  bool has_acy = (e_or.size() > 0);
  for (i=0; i<n; i++) {
    if (has_acy) {
      perm = tet_perm(T, i, e_or);
      psig = parity[perm]?-1:1;
    } 
    
    flat[PE3(perm,0)] = psig*fsh[3*i];
    flat[PE3(perm,1)] = psig*fsh[3*i+1]; 
    flat[PE3(perm,2)] = psig*fsh[3*i+2];
    
    get_zpq(flat, z, ip, iq); 
    rd = psig*rogers_dilog(z, -ip, iq);
    cv += rd/I; 
    
    if (has_acy) { print_perm(perm); }
    cout << " R(" << z << ',' << setw(2) << ip << ',' << 
      setw(2) << iq << ") " << rd << endl;
  }
  
  cout << "\nComplex vol: " << cv << endl; 
  double cs = cv.imag/(2*PI*PI);
  cs -= floor(2. * cs + 0.5)/2.;
  cout << "Chern-Simons (mod 1/2): " << cs << endl; 
}


/* conjecturally correct mod 1/12 (which is equivalent to Pi/6 in Walter's defn.) */ 

pari snap::raw_chern_simons() const
{
  pari corr = integer(face_mismatch_cs_term(manifold,0))/integer(24);
  pari n_shapes = integer(get_num_tetrahedra(manifold)); 

  return gimag(complex_volume())/(pTWO*pPI*pPI) + n_shapes/integer(12) - corr;
}


pari cusp_homology_kernel(pari mxc, pari mxe)
{

  /* 
     mxe = edge_homology_relations(manifold); 
     mxc = cusp_homology_relations(manifold); 

     writing r for the number of rows of mxc (or mxe which 
     is the same), H_1(complete) is the quotient of R^r by the subspace
     spanned by the columns of mxe. 

     [if the manifold were filled in some way then we would also have
     to quotient by the linear combinations of the columns of mxc
     corresponding to the filling curves.] 

     what we actually want is: to find which (linear combinations of the)
     columns of mxc can be expressed in terms of the columns of mxe. 
     these combinations are precisely the elements of the peripheral
     (cusp) homology which map to zero in H_1. we put the two sets of 
     columns side by side and use ker to find all the linear relations 
     between them. */ 

  pari k = ker(concat(mxc,mxe)); 

  /* we only want to know what combinations of the cols of mxc 
     can be expressed in terms of the mxe, not what the relation 
     with the mxe's is, so discard the rows corresponding to the mxe's. */ 
  pari mask = pow(pTWO, integer(length(mxc))) - pONE; 
  k = gtrans(extract(gtrans(k),mask)); 

  /* reduce to a linearly independent set */ 
  return image(k); 
}

/* we do the same calulation as above except that we want 
   to fill some of the cusps first and find the kernel in the peripheral
   homology of the remaining unfilled cusps. */ 

pari cusp_homology_kernel(pari mxc, pari mxe, pari sc)
{
  /* first count the number of filled cusps */ 
  int i, n = length(sc)/2; 
  int n_filled = 0; 
  for (i=0; i < n; i++)
    if (sc[i]!=pZERO || sc[i+n]!=pZERO) n_filled++; /* cusp is filled */ 
  
  /* make a mask to select only the columns of the 
     cusp homology relations which correspond to complete
     cusps. at the same time make a matrix which will 
     map the columns representing the meridian and longitude
     of each filled cusp to a single column representing the 
     filling curve. */ 
  int j=0; 
  pari mask, fc_mx = matrix(n_filled,2*n);  
  for (i=0; i < n; i++) {
    if (sc[i]==pZERO && sc[i+n]==pZERO) /* cusp is complete */ 
      mask += (pow(pTWO, integer(i)) + pow(pTWO, integer(i+n))); 
    else {
      fc_mx[j][i] = sc[i]; 
      fc_mx[j][i+n] = sc[i+n]; 
      j++; 
    }
  }

  /* put the columns of mxc corresponding to meridians
     and longitudes of complete cusps first, followed by 
     the columns corresponding to the filling curves on 
     the filled cusps, followed by the edge relations. */ 
  pari mx = concat(concat(extract(mxc,mask), mxc * fc_mx), mxe);

  pari k = ker(mx); 

  mask = pow(pTWO, integer(2 * (n - n_filled))) - pONE; 
  k = gtrans(extract(gtrans(k),mask)); 

  return image(k); 
}

pari snap::cusp_homology_kernel() const
{
  return ::cusp_homology_kernel(cusp_homology_relations(manifold),
				edge_homology_relations(manifold)); 
}


/* the manifold Y whose signature we compute here is obtained 
   from a cusped manifold (with fixed meridians and longitudes) 
   and a set of fillings for those cusps. 
   
   in fact in order to compute the signature we only need:
   1. a basis for the kernel of the map of the homology of 
   the cusps of the manifold into the homology of the manifold 
   itself (k0), and 
   2. the filling coefficients. 

   the signature of Y is then the signature of a certain form 
   on the space spanned by k0. 
*/ 

pari signature_Y(pari k0, pari sc)
{
  int i, n = length(sc)/2; 

  if (length(k0) != n) {
    warn("dimension of cusp homology kernel != number of cusps in signature_Y\n"); 
    return pZERO; 
  }
  pari k1 = matrix(n, 2 * n); 
  pari k2 = matrix(n, 2 * n); 
  for (i=0; i<n; i++) {
    k1[i][i+n] = pONE; 
    k2[i][i] = sc[i]; /* filling curve meridian */ 
    k2[i][i+n] = sc[i+n]; /* filling curve longitude */ 
  }

  /* now k0, k1 and k2 should all be matrices having 2*n rows and n columns. 
     they represent bases for various subspaces of the homology of the n cusps. 
     k0 is the kernel for the map of cusp homology into the homology of manifold. 
     k1 is the space spanned by the chosen longitudes. 
     k2 is the space spanned by the filling curves. 
   */ 

  pari k3 = intersect(image(concat(k1,k2)),k0); 
  if (length(k3)==0) return pZERO; 

  pari alpha = extract(gtrans(k3), pow(pTWO,integer(n)) - pONE);
  pari beta = extract(gtrans(k3), pow(pTWO,integer(n)) * (pow(pTWO,integer(n)) - pONE));
  pari mul = matrix(n,n); 
  for (i=0; i<n; i++) {
    if (sc[i].int_value() == 0) continue; 
    mul[i][i] = sc[i+n]/sc[i]; 
  }
  pari sig = signat(alpha * gtrans(beta - alpha * mul)); 
  
  /* signature_Y is signature(inner prod given by above matrix) */ 
  return (sig[0] - sig[1]); 
}


/* this is the signature of a 4-manifold Y obtained from 
   a filled 3-manifold. Y depends on the cusped 3-manifold,
   a choice of meridians and longitudes, and a filling of the 
   3-manifold. */

pari signature_Y(Triangulation *manifold)
{
  int n_cusps = get_num_cusps(manifold); 
  int n_filled = num_filled_cusps(manifold); 
  int n_comp = n_cusps - n_filled; 

  if (n_filled==0) 
    return pZERO; 

  pari mxc = cusp_homology_relations(manifold); 
  pari mxe = edge_homology_relations(manifold); 

  pari sc = surgery_coeffs2(manifold); 

  if (n_filled == n_cusps) 
    return signature_Y(cusp_homology_kernel(mxc,mxe),sc); 

  /* we were NOT given fillings for all of the cusps. 
     instead we compute a "relative signature":
     choose an arbitrary filling of the remaining cusps. 
     compute s1, the signature of the underlying cusped manifold
     with the specified filling of all the cusps. 
     fill just the cusps of the original manifold which were
     specified as filled in the first place to obtain a 
     manifold with fewer cusps. 
     compute s2, the signature of the fewer cusped manifold
     with the same arbitrary fillings of its cusps. 
     
     the theory implies that for large arbitrary fillings s1-s2 will 
     converge to a constant. experimentally, s1-s2 seems to be
     constant for any choice of the arbitrary filling. 

     notice that if we did this computation with specified fillings 
     for all the cusps then s2 would be the signature obtained 
     from a manifold with no cusps which is always zero. (the vector
     space on which the form whose signature we are computing lives, is 
     a subspace of the cusp homology.)

     if we did this computation with none of the cusps filled we
     would obtain the same value twice so s1-s2 would be zero. 
     */

  /* the "arbitrary" filling coeffs */ 
  pari arb = rvector(2 * n_comp); 
  int i; 
  for (i=0; i<n_comp; i++) {
    arb[i] = pONE; 
    arb[i+n_comp] = pTWO; 
  }

  /* all the surgery coeffs together */  
  pari sc_arb = gcopy(sc); 
  int l = 0; 
  for (i=0; i < n_cusps; i++) {
    if (cusp_is_complete(manifold,i)) {
      sc_arb[i] = arb[l]; 
      sc_arb[i+n_cusps] = arb[l+n_comp];
      l++; 
    }
  }
  return signature_Y(cusp_homology_kernel(mxc,mxe),sc_arb) - 
	 signature_Y(cusp_homology_kernel(mxc,mxe,sc),arb);  
}

int I_p_q(int p, int q)
{
  int quot, rem; 

  if (q == 0) return 0; 
  if (p == 0) return 2; /* perhaps this is meaningless?? */ 
  
  if (p < 0) {
    if (q < 0) return I_p_q(-p,-q);
    else return -1 - I_p_q(-p,q);
  } else {
    if (q < 0) return -1 - I_p_q(p,-q); 
  }

  quot = p / q;
  rem = p % q; 

  if (rem == 0) return 2 - quot; 

  return (2 - quot + I_p_q(q, q - rem)); 
}

/* this is 12 p times the classical dedekind sum */ 

int dedekind(int q, int p)
{
  if (p==0) return 0; 

  if (p < 0) return -dedekind(-q,-p); 

  double sum = 0.0;
  double dp = (double)p; 
  int i; 

  for (i=1; i<p; i++) {
    sum += 1.0/(tan((double)i * PI/dp) * tan(PI * (double)((i * q) % (2 * p))/dp));
  }
  return (int)floor(3.0 * sum + .5); 
}

pari eta_rat_part(pari k0, pari sc)
{
  int i, n_cusps = length(sc)/2;
  pari three = integer(3); 
  pari cusp_sum;
  for (i=0; i<n_cusps; i++) {
    if (sc[i]!=pZERO) {
      cusp_sum += (integer(dedekind(sc[i+n_cusps].int_value(),sc[i].int_value())) - sc[i+n_cusps]) / (three * sc[i]);
    } else {
      cusp_sum += -pONE;
    }
  }
  return cusp_sum - signature_Y(gtrans(k0), sc); 
}

pari eta_rat_diff(pari k0, pari sc, pari cusp_maps)
{
  int i, n = length(sc)/2; 

  pari transform = matrix(2*n, 2*n); 
  for (i=0; i<n; i++) {
    transform[i][i] = cusp_maps[i][0][0]; 
    transform[i+n][i] = cusp_maps[i][0][1]; 
    transform[i][i+n] = cusp_maps[i][1][0]; 
    transform[i+n][i+n] = cusp_maps[i][1][1]; 
  }

  return eta_rat_part(k0 * transform, sc * transform) - eta_rat_part(k0,sc); 
}

pari cusp_holonomy_sum(const pari& h, const pari& sc)
{
  pari the_sum; 
  int i, n = length(sc)/3; 
  for (i=0; i<n; i++) {
    if (sc[i+2*n]==pZERO) continue; 

    if (sc[i]!=pZERO)
      the_sum += -gimag(h[i+n])/sc[i];
    else 
      the_sum += gimag(h[i])/sc[i+n]; 
  }
  return the_sum/(integer(6)*pPI); 
}

pari cusp_dedekind_sum(const pari& sc)
{
  pari the_sum, three = integer(3); 
  int i, n = length(sc)/3; 
  for (i=0; i<n; i++) {
    if (sc[i+2*n]==pZERO) continue; 

    if (sc[i]!=pZERO)
      the_sum += integer(dedekind(sc[i+n].int_value(),sc[i].int_value()))/
	(three * sc[i]);
    else 
      the_sum += -pONE;
  }
  return the_sum; 
}

pari cusp_surgery_sum(const pari& sc)
{
  pari the_sum, three = integer(3); 
  int i, n = length(sc)/3; 
  for (i=0; i<n; i++) {
    if (sc[i+2*n]==pZERO) continue; 

    if (sc[i]!=pZERO)
      the_sum += -sc[i+n]/(three * sc[i]);
    else 
      the_sum += sc[i]/(three * sc[i+n]);
  }
  return the_sum; 
}

pari torsion(const pari& h, const pari& sc)
{
  int i, n_cusps = length(h)/2;
  pari cg = rvector(n_cusps); 
  pari u, v; 

  for (i=0; i<n_cusps; i++) {
    if (sc[i+2*n_cusps]==pZERO) continue;

    u = h[i]; v = h[i+n_cusps]; 

    if (abs(greal(v)) < pari(1e-5))
      cg[i] = (greal(u) * gimag(v) - greal(v) * gimag(u)) * 
	      gimag(v) / (- pTWO * pPI * greal(v)); 
    else 
      cg[i] = -gimag(u); 
  }
  return cg; 
}

pari torsion_sum(const pari& h, const pari& sc)
{
  pari the_sum;
  pari tor = torsion(h, sc); 
  int i, n = length(tor); 
  for (i=0; i<n; i++) 
    the_sum += tor[i]; 

  return the_sum/(integer(-6) * pPI); 
}

pari defect_sum(const pari& sc)
{
  pari the_sum; 
  int i, n = length(sc)/3;

  for (i=0; i<n; i++)
    the_sum += I_p_q(sc[i].int_value(),sc[i+n].int_value()); 

  return the_sum; 
}

pari core_sum(const pari& cg)
{
  int i, n_cusps = length(cg); 
  pari core_len_sum; 

  for (i = 0; i < n_cusps; i++)
    core_len_sum += cg[i];

  return core_len_sum; 
}

void snap::check_eta() const
{
  pari shape_sum = -greal(dilog_sum()) / (integer(3) * pPI * pPI); 
  printf("shape sum = ");
  shape_sum.print(); 

  pari h = acc_holonomies(); 
  pari sc = ::surgery_coeffs(manifold); 
  pari holonomy_sum = cusp_holonomy_sum(h,sc); 
  printf("holonomy sum = ");
  holonomy_sum.print(); 

  pari ded_sum = cusp_dedekind_sum(sc); 
  printf("dedekind sum = "); 
  ded_sum.print();

  pari surgery_sum = cusp_surgery_sum(sc); 
  printf("surgery sum = "); 
  surgery_sum.print();

  pari sig = signature_Y(manifold); 
  printf("signature term (to subtract) = "); 
  sig.print(); 

  printf("Eta difference = "); 
  (shape_sum + holonomy_sum + ded_sum + surgery_sum - sig).print();  
/* 
  printf("\n");
  pari tor = torsion_sum(h, sc); 
  printf("torsion sum = "); 
  tor.print(); 

  printf("Eta difference (2nd formula) = "); 
  (shape_sum + tor + ded_sum - sig).print(); 

  printf("\n");
  pari core_len_sum = gimag(core_sum(core_geodesics()))/(integer(-6) * pPI); 
  printf("torsion from core geodesics = "); 
  core_len_sum.print(); 

  pari def = defect_sum(sc); 
  printf("defect sum = "); 
  def.print(); 

  printf("Eta difference (3rd formula) = "); 
  (shape_sum + core_len_sum + def - sig).print(); 
*/
}

bool snap::eta_available() const
{
  return (_eta_invariant_known || _eta_fudge_known); 
}

/* eta_invariant(manifold) = eta_fudge(triangulation) + eta_difference() */ 

pari snap::eta_difference() const
{
  if (!accurate_info_set()) {
    warn("eta_difference() requires accurate shape info\n");
    return pZERO; 
  }
  if (cone_surgery(manifold)) {
    warn("eta_difference probably wrong for cone manifold\n"); 
  }

  pari sc = ::surgery_coeffs(manifold); 

  pari three = integer(3); 
  pari shape_sum = -greal(dilog_sum()) / (three * pPI * pPI); 

  int i; 
  pari cusp_sum;
  pari two_pi = pTWO * pPI; 
  int n_cusps = get_num_cusps(manifold); 
  pari h = acc_holonomies(); 
  for (i=0; i<n_cusps; i++) {
    if (cusp_is_complete(manifold,i)) continue; 
    if (sc[i]!=pZERO) {
      cusp_sum += (-gimag(h[i+n_cusps])/two_pi + 
	integer(dedekind(sc[i+n_cusps].int_value(),sc[i].int_value())) - 
		   sc[i+n_cusps]) / (three * sc[i]);
    } else {
      cusp_sum += gimag(h[i]) / (two_pi * three * sc[i+n_cusps]); /* - pONE */
    }
  }

/* Ouyang's paper says there should be a -pONE here. If I put it in my 
   answers are off by 1. However it is not clear whether the problem is
   really here or in the signature computation. (Nor have I established
   for certain that the answer will be off by exactly 1 for each 0,1 
   filled cusp in every case. This is true for all the redundant cases
   in the surgery relations between standard manifolds */ 

  /* get the sum of the complex lengths of all core geodesics */
  pari sigy = signature_Y(manifold); 

#if 0
  printf("\n"); 
  printf("cusp sum: "); cusp_sum.print(); 

  printf("shape sum: "); shape_sum.print();
  printf("signature: "); sigy.print();
#endif 

  return shape_sum + cusp_sum - sigy; 
}

pari factor_sl2Z_matrix(pari mx)
{
  pari S = p_lisexpr("[0,1;-1,0]"); 
  pari Sinv = p_lisexpr("[0,-1;1,0]"); 
  pari T = p_lisexpr("[1,0;0,1]");
  pari Tinv = p_lisexpr("[1,0;0,1]"); 

  pari t; 

  pari prod = rvector(vecmax(abs(mx[0])).int_value()+5);
  int j, i = 0; 

  while (1) {
    if (mx[0][1]==pZERO) { /* [a,b;0,d] */ 
      if (mx[0][0]>pZERO) { /* a>0 */ 
	if (mx[1][0]==pZERO) break; 
	prod[i++] = mx;
	break; 
      } else { /* a<0 */ 
	prod[i++] = S;
	mx *= Sinv; 
      }
    } else { /* [a,b;c,d] */ 
      if (abs(mx[1][1])<abs(mx[0][1])) { /* |d|<|c| */ 
	if (mx[1][0]>pZERO) { /* swap making a>0 */ 
	  prod[i++] = S; 
	  mx *= Sinv; 
	} else {
	  prod[i++] = Sinv; 
	  mx *= S; 
	} 
      } else { /* |c|<|d|, subtract euc quot from d */ 
	t = gdivround(mx[1][1],mx[0][1]);
	T[1][0] = t; Tinv[1][0] = -t; 
	prod[i++] = gcopy(T);
	mx *= Tinv; 
      }
    }
  }

  pari pr = rvector(i); 
  for (j=0;j<i;j++) 
    pr[j] = prod[i-j-1];
  return pr; 
}

pari cusp_term_adjustment(pari mx, pari pq)
{
  pari l = factor_sl2Z_matrix(mx); 
  pari three = integer(3); 
  pari adj; 

  int i; 
  for (i = length(l)-1; i>=0; i--) {
    if (l[i][0][0]==pZERO) adj -= gsigne(pq[0]*pq[1]); 
    else adj -= (l[i][1][0]/three + 
		 integer(gsigne(pq[0]*pq[1]) - gsigne((pq[0]+l[i][1][0]*pq[1])*pq[1]))); 
    pq = l[i] * pq; 
  }
  return adj; 
}

Boolean re_fill(Triangulation *drilled, Triangulation *original)
{
    const int       num_steps = 5;
    const double    coef[5] = {5.0, 3.0, 2.0, 1.3, 1.0};
    int i, k; 
    int             singularity_index;
    Complex         core_length;

#if 0
    FILE            *core_report;

    core_report = fopen("core_report", "w");
    if (core_report == NULL) {
      warn("couldn't open core report\n");
      return FALSE;
    }
#endif

    for (i = get_num_cusps(drilled); --i >= get_num_cusps(original); ) {
      for (k = 0; k < num_steps; k++) {
	set_cusp_info(drilled, i, FALSE, coef[k], 0.0);
	switch (do_Dehn_filling(drilled)) {
	case geometric_solution:
	case nongeometric_solution:
	  /* OK */
	  break;
	case degenerate_solution:
	  warn("degenerate solution\n");
	  /* fall through to default . . . */
	default:
	  warn("can't find hyperbolic structure\n");
	  // fclose(core_report);
	  return FALSE; 
	}
      }
      
      core_geodesic(drilled,
		    i,
		    &singularity_index,
		    &core_length,
		    NULL);
      if (singularity_index != 1) {
	warn("singularity index is not 1\n");
	// fclose(core_report);
	return FALSE; 
      }
      // fprintf(core_report, "core length %d:  %9.6lf + i %9.6lf\n",
      //   i, core_length.real, core_length.imag);
    }
    // fclose(core_report);

    /*
     *  Check that volume agrees with that of the original manifold.
     */
    if (1e-10 < fabs(volume(drilled, NULL) - volume(original, NULL))) {
      warn("filled drilled volume differs from original volume\n");
      return FALSE; 
    }

    return TRUE; 
}


/* "drilled" and "standard" are isometric manifolds, triangulated differently. 
   Filling "drilled" will give us a manifold isometric to "original". 
   Since we know the CS value of "standard" we can copy that to "drilled". 
   Then we can fill "drilled" to get a manifold isometric to "original" whose
   CS value is known. 
*/ 

Boolean compute_CS_value(
    Triangulation *original, 
    Triangulation *drilled, 
    Triangulation *standard, 
    int orientation_reversed)
{
    Boolean         value_is_known;
    double          the_value;
    int             the_precision;
    Boolean         requires_initialization;

    /*
     *  Set CS value of the drilled manifold from value of standard manifold. 
     *  This also sets the CS fudge which remains invariant under filling.
     */
    get_CS_value(standard,
		 &value_is_known,
		 &the_value,
		 &the_precision,
		 &requires_initialization);
    if (value_is_known == FALSE) {
      printf("no CS value in standard manifold\n");
      return FALSE; 
    }
    if (the_precision < 8) {
      printf("imprecise CS in standard manifold\n");
      return FALSE; 
    }
    if (orientation_reversed)
      the_value = -the_value;
    set_CS_value(drilled, the_value);

    /*
     *  Re-fill the drilled manifold. 
     */
    if (!re_fill(drilled, original)) return FALSE;

    /* 
     * CS value should now be available for manifold. Transfer it to 
     * the original manifold. 
     */
    get_CS_value(drilled,
		 &value_is_known,
		 &the_value,
		 &the_precision,
		 &requires_initialization);
    if (value_is_known == FALSE) {
      warn("no CS for presumed_original\n");
      return FALSE; 
    }
    if (the_precision < 8) {
      warn("imprecise CS for presumed_original\n");
      return FALSE; 
    }
    set_CS_value(original, the_value);

    return TRUE; 
}


bool snap::bootstrap_CS_value(double& value, int report)
{
  /* Drill to a manifold in the standard set. */ 
  Triangulation *drilled; 
  IsometryList *isometry_list = 0; 
  Triangulation *standard;
  int which; 
  if (!drill_to_standard(manifold, &drilled, &isometry_list, &standard, which, report)) 
    return false; 

  /* Check if the isometries are all orientation reversing. */ 
  Boolean or_pres, or_rev; 
  isometry_list_orientations(isometry_list, &or_pres, &or_rev);
  free_isometry_list(isometry_list);

  /* Set the chern simons value for the manifold from 
     that of the standard manifold. */ 
  compute_CS_value(manifold, drilled, standard, !or_pres); 

  int precision; 
  Boolean known, req_init; 
  get_CS_value(manifold, &known, &value, &precision, &req_init);

  return known;
}

bool snap::set_eta_info()
{
  if (!eta_available()) return false; 

  if (!_eta_invariant_known) {
    _eta_invariant = eta_invariant(); 
    _eta_invariant_known = true; 
  } else if (!_eta_fudge_known) {
    if (flattening.type()==t_INT) {
      find_a_flattening(); 
      if (flattening.type()==t_INT) return false; 
    }
    _eta_fudge = eta_fudge(); 
    _eta_fudge_known = true; 
  }
  return true; 
}

void snap::print_eta_info() const
{
  printf("Eta fudge: "); 
  if (_eta_fudge_known) {
    _eta_fudge.print(); 
  } else {
    printf("not known.\n"); 
  }
  printf("Eta invariant: "); 
  if (_eta_invariant_known) {
    _eta_invariant.print(); 
  } else {
    printf("not known.\n"); 
  }
}

void snap::set_eta_invariant(pari x)
{
  _eta_invariant = x; 
  _eta_invariant_known = true; 
  _eta_fudge_known = false; 
}

void snap::set_eta_fudge(pari x)
{
  _eta_fudge = x; 
  _eta_fudge_known = true; 
  _eta_invariant_known = false; 
}

void snap::set_flattening(pari f)
{
  flattening = f; 
}

pari snap::eta_invariant() const
{
  if (_eta_invariant_known) return _eta_invariant; 
  if (_eta_fudge_known) return _eta_fudge + eta_difference();
  return pZERO;
}

pari snap::eta_fudge() const
{
  if (_eta_fudge_known) return _eta_fudge; 
  if (_eta_invariant_known) {
    pari f = _eta_invariant - eta_difference();
    pari fx = ground(integer(18)*f)/integer(18); 
    if (abs(f-fx) > pari(1e-40)) {
      printf("Warning: eta_fudge not a multiple of 1/18.\n"); 
      printf("value found: "); f.print();
      return f; 
    }
    return fx; 
  }
  return pZERO;
}

#if 0
/* dilog varies as follows: its principal branch is cut along
   [1,+inf]. analytically continuing across this in a downward 
   (ie. imaginary part becoming negative) direction we have to 
   add 2 pi i log(z) to the principal branch of the dilog. 
   since other branches of dilog involve this log term they 
   are defined on the plane cut along [-inf,0] and [1,+inf]. 
   for a branch with n*(2 pi i log(z)) added, analytically continuing
   across [-inf,0] in a downward direction we have to add 2 pi i to 
   the principal branch of log and therefore subtract n*(4 pi^2) from 
   the branch of dilog we were on. analytically continuing across
   [1,+inf] simply changes the n in front of the log term. 
   therefore the various branches of the dilog function are given by 
   dilog(z) + n*(2 pi i log(z)) - m*(4 pi^2) for integers n,m. 
*/

pari dilog_with_history(pari z, ShapeInversion *history)
{
  int n=0,m=0,down=1;

  while (history) {
    switch (history->wide_angle) {
    case 0: /* crossing [-inf,0] */ 
      m += down * n; 
      break;
    case 1: /* crossing [1,+inf] */ 
      n += down; 
      break; 
    default:
      break;
    }
    down = -down; 
    history = history->next; 
  }
  return dilog(z) + integer(n)*pTWO_PI*gi*log(z) - integer(m)*pFOUR*pPI*pPI;
}
#endif
