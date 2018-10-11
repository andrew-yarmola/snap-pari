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
#include "field.hh"

#include <algorithm>

extern "C" {
#include "pari_ext.h"
#include <signal.h>
}

// static prt_init s1("field.cc\n"); 
// static say xx("field.cc"); 

inline bool is_real(const pari& r)
{ return numerical_zero(gimag(r)); }

inline pari lindeppart1(pari::tmp x, long bit)
{ return lindeppart1(x.g, bit); }

inline pari lindeppart2(pari::tmp part, pari::tmp x, long bit)
{ return lindeppart2(part.g, x.g, bit); }

const pari_printform pp_short = pari_printform(1,1,'g',16,0); 

const pari_printform pp_tiny = pari_printform(1,0,'f',9,0); 

const pari_printform no_newline = pari_printform(1,0,'g',-1,0);

// result is monic
pari polymod_algdep(const pari& elt)
{
  pari cp = caradj0(elt, 0); 
  pari denom = ggcd(cp, deriv(cp,0));  

  return (cp * denom[lgeff(denom)-1])/denom; 
}

int polymod_lindep(pari vec, pari& dependence)
{
  int i, j; 
  int deg = lgeff(vec[0][0]) - 1; 
  int cols = length(vec); 
  pari mat = matrix(cols, deg+1); 
  for (i = 0; i < cols; i++) {
    for (j = 0; j < lgeff(vec[i][1]); j++) {
      mat[i][j] = vec[i][1][j];
    }
  }

  /* ker finds any linear relations between the columns of mat */ 
  /* matrixqz3 turns a rational matrix into an integer matrix with same Q-span */ 
  pari k = matrixqz3(ker(mat)); 
  if (length(k)==0) return 0; 

  int shortest = deg, shortcol = 0; 
  for (i = 0; i < length(k); i++) {
    /* find length of column i */ 
    for (j = cols-1; j >=0; j--) {
      if (k[i][j] != pZERO) break; 
    }
    if (j < shortest) { 
      shortest = j;
      shortcol = i; 
    }
  }
  dependence = k[shortcol]; 
  return 1; 
}

int polymod_contains(const pari& elt, const pari& moos)
{
  pari dep; 
  return polymod_lindep(concat(moos, -elt), dep); 
}

int polymod_contains(const pari& elt, const pari& moos, const pari& mp, pari& exact)
{
  pari dep; 
  if (!polymod_lindep(concat(moos, -elt), dep)) return 0; 

  int n = length(moos); 
  pari coeffs = extract(dep, pow(pTWO,integer(n))-pONE)/dep[n]; 
  exact = polymod(gtopolyrev(coeffs,0), mp); 

  return 1; 
}

/* gen is a polymod. we find its minimum polynomial and 
   then put powers 0-(degree-1) of gen in moos[1]-moos[degree]. */ 

void set_moos(const pari& gen, pari& moos, pari& mp)
{
  /* look for min poly for first generator */ 
  mp = polymod_algdep(gen); 

  int degree = lgeff(mp)-1;

  /* moos contains powers of primitive element */ 
  moos = cvector(degree);
  moos[0] = polymod(pari(polun[0]),gen[0]);
  if (degree > 1) moos[1] = gen; 
  int c; 
  for(c=2; c<degree; c++)
    moos[c] = moos[c-1] * gen;
}



int find_min_poly(pari elt, pari& poly, int report)
{
  /* do LLL */ 
  poly = algdep2(elt,field::maxdeg,(long)(1+digits*field::fudge));

  if (report) {
    printf("LLL found: "); 
    poly.print();
  }

  /* factorize LLL result */ 

  #ifdef PARI_2_2_OR_LATER
  pari factors = factpol(poly,1); 
  #else
  pari factors = factpol(poly,0,1); 
  #endif

  int n_fact = length(factors[0]); 
  if (report && n_fact > 1) {
    printf("Factors: ");
    factors.print(); 
  }

  /* identify the factor */ 
  int i;
  pari error; 
  for (i = 0; i < n_fact; i++) {
    error = abs(gsubst(factors[0][i],0,elt));
    if (numerical_zero(error)) break; 
  }
  if (i == n_fact) {
    if (report) {
      if (n_fact > 1) 
	printf("Generator does not seem to be a root of any factor!\n"); 
      else 
	printf("Generator does not seem to be a root of this polynomial!\n"); 
    }
    return 0; /* means we didn't find a polynomial */
  }
  poly = factors[0][i];
  if (report) {
    if (n_fact> 1) {
      printf("Generator seems to be a root of factor: ");
      poly.print(); 
    }
    printf("Absolute value of factor(generator): ");
    error.print(pp_short); 
  }

  // Make sure it is monic. 
  poly /= poly[lgeff(poly)-1];

  return 1; 
}

int is_integral(const pari& elt)
{
  return content(caradj0(elt,0))==pONE;
}

#if 0
void cancel_lead_coeff(pari& pol)
{
  int degree = lgeff(pol)-1, i;
  if (pol[degree]==pONE) return; /* perhaps we can get out quick */
  for (i = 0; i < degree; i++) 
    if (gmod(pol[i],pol[degree])!=pZERO) return;
  for (i = 0; i <= degree; i++)
    pol[i] /= pol[degree]; 
}

void make_monic(pari& pol)
{

  /* to get a monic poly defining the same field as
     an x^n + ... + a0
     we can substitute x/an and multiply through by 
     (an)^(n-1). this is equivalent to multiplying the 
     first coeff by an^-1, the second by 1, the third by 
     an, etc. */

  pol /= content(pol); 

  int degree = lgeff(pol)-1, i;
  pari lead = pol[degree]; 
  if (lead==pONE) return; 
  pari mult = pONE/lead; 
  for (i=degree; i>=0; i--) {
    pol[i] *= mult;
    mult *= lead; 
  }
}
#endif

pari poly_denominator(const pari& poly)
{
  int deg = lgeff(poly)-1; 

  // Polynomial should be monic so c <= 1. 
  pari c = content(poly); 
  if (c==pONE) return c; 

  pari primes = factor(c);
  int npr = length(primes[0]); 

  int p, i, j, n; 
  pari coeff_fact; 

  pari req_power; 

  for (i=0; i<npr; i++) primes[1][i] = pZERO; 

  for (p = 0; p < deg; p++) {

    coeff_fact = factor(poly[p]); 
    n = length(coeff_fact[0]); 

    i = 0; 
    for (j=0; j<n; j++) {

      // Skip positive powers. 
      if (coeff_fact[1][j] >= pZERO) continue; 

      // Step i to the next prime in coeff_fact. 
      while (primes[0][i] < coeff_fact[0][j] && i < npr) i++; 
      if (i>=npr) break; 

      // Check that primes[0][i] == coeff_fact[0][j]. 
      if (primes[0][i] > coeff_fact[0][j]) {
	printf("Oops, problem in poly_denominator.\n"); 
	return pONE; 
      }

      // Get the smallest power of this prime required to make the
      // power of this prime positive for this coefficient. 
      req_power = ceil((-coeff_fact[1][j])/integer(deg - p)); 
      if (primes[1][i] < req_power) primes[1][i] = req_power; 

    }
  }

  pari prod = pONE; 

  for (i=0; i<npr; i++) prod *= pow(primes[0][i], primes[1][i]); 

  return prod;
}

pari reduced_poly(const pari& p, pari& newroot_on_old, int report)
{
  int deg = lgeff(p) - 1; 

  pari den = poly_denominator(p); 
  pari poly = pow(den, integer(deg)) * gsubst(p, 0, pari(polx[0])/den); 
  // root of poly is den * root of p. 

  if (report) printf("Finding a reduced polynomial...\n"); 

  pari info = polred0(poly, 2, pZERO); 

  int i; 
  for (i=0; i<length(info[1]); i++) {
    if (lgeff(info[1][i])-1 == deg) break; 
  }

  if (i==length(info[1])) {
    printf("Problem finding reduced polynomial.\n"); 
    newroot_on_old = pari(polx[0]); 
    return poly;
  }

  newroot_on_old = polymod(gsubst(info[0][i],0,den * pari(polx[0])), p); 
  
  return info[1][i];
}


pari canonical_poly(const pari& p, pari& newroot_on_old, int report)
{
  pari poly = p; 
  int deg = lgeff(poly) - 1; 

  pari den = poly_denominator(p); 
  poly = pow(den, integer(deg)) * gsubst(p, 0, pari(polx[0])/den); 
  // root of poly is den * root of p. 

  if (report) printf("Seeking T2 minimal generator...\n"); 

  pari info = polredabs0(poly,5); 
  int num = length(info); 
  // 5 is a flag which says get all the polynomials, and for each one
  // get old_root_on_new. 

  pari fix_sign = rvector(num); 
  pari candidates = rvector(num); 
  pari pol; 

  int i, j; 
  for (j=0; j<num; j++) {
    pol = info[j][0];
    candidates[j] = pol; 

    // Fix the sign of coeffs of powers of x differing from degree by an odd number.
    fix_sign[j] = pONE; 
    for (i = deg-1; i >= 0; i -= 2)
      if (pol[i] != pZERO) break; 

    if (i < 0 || gsigne(pol[i]) < 0) continue;

    // Ok, something to fix. 
    for (; i >= 0; i -= 2) {
      pol[i] = -pol[i]; 
    }

    candidates[j] = pol; 
    fix_sign[j] = -pONE; 
  }

  if (report) {
    printf("T2 minimal equivalents: "); 
    candidates.print(); 
  }


  int which = 0; 
  if (num > 1) {

    pari mat = matrix(num, deg+4); 

    /* Make a matrix whose cols are |disc(f)|,|f[n]|,...,|f[0]|,sigsum for each candidate.
       sigsum is the sum from i=0,...,n of 2^n * sign(f[n]). This ensures a canonical
       choice of polynomial when the |f[i]| are all the same. */ 

    pari pow_two; 

    for (i = 0; i < num; i++) {
      pow_two = pONE; 
      mat[i][0] = abs(discsr(candidates[i]));
      mat[i][deg+3] = integer(i);
      for (j = 0; j < deg+1; j++) {
	mat[i][j+1] = abs(candidates[i][deg-j]);
	mat[i][deg+2] += integer(gsigne(candidates[i][j])) * pow_two;
	pow_two *= pTWO; 
      }
    }
    
    /* Make a vector 1,2,...,deg+3. */ 
    pari vec = rvector(deg+3); 
    for (i = 0; i < deg+3; i++)
      vec[i] = i+1; 
    
    /* Sort lexicographically on first deg+3 rows of mat 
       then extract required poly from first row. */ 

    pari sorted = vecsort0(mat,vec,0); 
    which = sorted[0][deg+3].int_value();

  }

  pari newpoly = candidates[which];
  newroot_on_old = polymodrecip(info[which][1]/den) * fix_sign[which];

  if (report) {
    printf("Canonical polynomial: "); 
    newpoly.print();
  }

  return newpoly; 
}

/* implementation of the field class */ 

// field rationals; 
// field complex_infinity; 

int field::maxdeg = 16;
double field::fudge = 0.6; 

field::field() 
: nf(initalg(pari(polx[0]) - pONE)), root_num(0), _is_canonical(true)
{
  set_root_indices();
  numeric_ibasis = gsubst(nf[6],0,pZERO);
  ibasis_lllmat = lindeppart1(numeric_ibasis,(long)(1+digits*fudge));
}

field::field(const field& f) 
: nf(f.nf), root_num(f.root_num), _is_canonical(f._is_canonical),
  numeric_ibasis(f.numeric_ibasis), ibasis_lllmat(f.ibasis_lllmat)
{
  set_root_indices();
}

field::field(const pari& _min_poly, int _root_num, bool isc) 
: nf(initalg(_min_poly)), root_num(_root_num), _is_canonical(isc)
{
  set_root_indices();
  numeric_ibasis = gsubst(nf[6],0,root());
  ibasis_lllmat = lindeppart1(numeric_ibasis,(long)(1+digits*fudge));
}

field::field(const pari& spec) 
{
  pari mp; 
  if (length(spec)==4) {
    mp = spec[1]; 
    root_num = spec[2].int_value(); 
  } else { // old format
    mp = spec[0]; 
    root_num = spec[1].int_value();
  }
  _is_canonical = (leadingcoeff(mp)==pONE);
  nf = initalg(mp * (_is_canonical ? pONE : -pONE));
  set_root_indices();
  numeric_ibasis = gsubst(nf[6],0,root());
  ibasis_lllmat = lindeppart1(numeric_ibasis,(long)(1+digits*fudge));
}

field::~field()
{}

void field::clear()
{
  nf = initalg(pari(polx[0]) - pONE);
  root_num = 0;
  _is_canonical = true; 
  set_root_indices(); 
  numeric_ibasis = gsubst(nf[6],0,pZERO);
  ibasis_lllmat = lindeppart1(numeric_ibasis,(long)(1+digits*fudge));
}

void field::update_precision()
{
  if (!is_set()) return; 
  nf = initalg(nf[0]);
  numeric_ibasis = gsubst(nf[6],0,root());
  ibasis_lllmat = lindeppart1(numeric_ibasis,(long)(1+digits*fudge));
}

field& field::operator = (const field& f)
{
  nf = f.nf; 
  _is_canonical = f._is_canonical; 
  set_root_indices(); 
  root_num = f.root_num; 
  numeric_ibasis = f.numeric_ibasis; 
  ibasis_lllmat = f.ibasis_lllmat; 

  return *this;
}

int field::is_real() const
{ 
  pari sig = nf[1]; 
  return root_num > 0 && root_num <= sig[0].int_value(); 
}

bool field::amphicheiral() const
{ 
  return contains(gconj(root())); 
}

pari field::spec() const
{
  pari res=rvector(3); 
  res[0] = _is_canonical ? nf[0] : -nf[0];
  res[1] = nf[1]; 
  res[2] = nf[2]; 
  return res; 
}

// We shall sort the roots of a polynomial into a 
// canonical order, using the following order relation. 
// Real roots are smaller than complex roots. 
// If both are real or both are complex, comparison is
// the usual ordering on the real component. 
// For complex roots whose real components coincide, 
// the one whose imaginary part has smaller absolute 
// value is the lesser. (Conjugate roots are considered 
// equal.)
//
bool root_less(const pari& a, const pari& b)
{
  if (is_real(a) != is_real(b)) {
    return is_real(a);
    // a is smaller if a is real and b isn`t. 
  }

  // Here both are real or both are complex. 

  pari diff = greal(a) - greal(b);
  if (!numerical_zero(diff)) {
    return gsigne(diff)==-1; 
  }

  // Here real parts were the same. 
  // If (both) real then return false. (Ie. a is not less than b.)

  if (is_real(a)) return false; 

  // Finally compare abs value of imaginary parts. 
  diff = abs(gimag(a)) - abs(gimag(b)); 
  if (numerical_zero(diff)) return false; 
  return gsigne(diff)==-1; 
}


void field::set_root_indices()
{
  int i, n = length(nf[5]); 
  root_index = vector<int>(n); 
  for (i=0; i<n; i++) root_index[i] = i;
  sort(root_index.begin(), root_index.end(), root_index_less(nf[5]));
}


pari field::root(int r) const
{
  int n = r > 0 ? r : -r; // n = abs(r);
  int r1 = nf[1][0].int_value();
  pari rt = nf[5][root_index[n-1]];
  if (n <= r1) return greal(rt);

  // The == in the next line works as a logical not exclusive or. 
  if ((r > 0) == (gimag(rt) > pZERO))
    return rt; 

  return gconj(rt); 
}

int field::root_number(const pari& r) const
{
  int i; 
  int r1 = nf[1][0].int_value();
  int r2 = nf[1][1].int_value(); 
  pari ri; 
  for (i = 1; i <= r1+r2; i++) {
    ri = nf[5][root_index[i-1]];
    if (numerical_zero(r-ri) ||
	(i > r1 && numerical_zero(r-gconj(ri)))) {
      if (i <= r1) {
	return i; 
      } else { 
	return i * gsigne(gimag(r)); 
      }
    }
  }

  return 0; 
}

pari field::numeric_value(pari element, int n) const
{
  return gsubst(element[1], 0, n ? root(n) : root()); 
}

int field::contains(pari x) const 
{
  if (numerical_zero(x)) return 1; 
  if (x.type()==6) {
    if (numerical_zero(gimag(x))) x = greal(x);
    else if (is_real()) return 0; 
  }

  pari dep = lindeppart2(ibasis_lllmat, -x,(long)(1+digits*fudge));
  return numerical_zero(concat(numeric_ibasis,-x) * gtrans(dep));
}

int field::contains(pari x, pari& exact_x, int report) const 
{
  pari ex; 
  int tx = x.type(); 
  if (tx==17||tx==18) { // row or column vector.
    int i, n = length(x);
    exact_x = tx==17 ? rvector(n):cvector(n); 
    for (i=0; i<n; i++) {
      if (!contains(x[i],ex,report)) break;
      exact_x[i] = ex; 
    }
    return i==n;
  }

  if (numerical_zero(x)) {
    exact_x = polymod(pZERO, nf[0]);
    if (report) {
      printf("Element: 0. Approximate error: "); 
      abs(x).print(pp_short); 
    }
    return 1; 
  }

  if (x.type()==6) {
    if (numerical_zero(gimag(x))) x = greal(x);
    else if (is_real()) {
      if (report) printf("Dependence not found.\n"); 
      return 0; 
    }
  }

  pari dep = lindeppart2(ibasis_lllmat,-x,(long)(1+digits*fudge));
  int n = length(numeric_ibasis); 
  pari ibasis_x = gtrans(extract(dep, pow(pTWO,integer(n))-pONE));

  pari error = abs(concat(numeric_ibasis,-x) * gtrans(dep));

  if (dep[n]==pZERO) {
    if (report) {
      printf("Element: "); 
      ibasis_x.print(no_newline);
      printf("/0.\n"); 
    }
    return 0; 
  } 

  ibasis_x /= dep[n];

  if (report) {
    printf("Element: "); 
    ibasis_x.print(no_newline);
    printf(". Approximate error: "); 
    error.print(pp_short); 
  }

  if (!numerical_zero(error)) return 0; 

  exact_x = basistoalg(nf, ibasis_x);

  return 1;
}

int field::extend(pari const& rt1, pari& exact_rt1, int cp, int report, int gnum)
{
  if (!root_num) 
    root_num = 1; // set it to be the rationals initially

  if (contains(rt1, exact_rt1, is_rational() ? 0 : report)) 
    return 1; // field unchanged

  // No rational elt should get past this point!!

  pari ex; // Dummy. 

  pari rt0 = root();
  int deg = lgeff(nf[0]); 

  if (report && !is_rational()) {
    printf("Field does not appear to contain generator %d\n", gnum); 
    if (2*degree() > maxdeg) return 0; 
    printf("Trying field generated by this instead...\n"); 
  }

  pari poly; 
  if (!find_min_poly(rt1, poly, report)) {
    if (report) printf("Unable to find min poly for generator %d!\n", gnum); 
    return 0; 
  }

  if (lgeff(poly) > deg) {
    if (!set_field(poly, rt1, exact_rt1, cp, report)) return 0; 
    if (contains(rt0, ex, report)) return 2; 
  }

  pari newgen = rt0 + rt1; 
  if (report) 
    printf("Trying field generated by generator %d + previous...\n", gnum); 
  if (!find_min_poly(newgen, poly, report)) {
    if (report) printf("Unable to find min poly for generator %d + previous!\n", gnum); 
    return 0;
  }

  if (lgeff(poly) > deg) {
    if (!set_field(poly, newgen, ex, cp, report)) return 0; 
    if (contains(rt1, exact_rt1, report)) return 2; 
  }

  if (report) printf("Unable to find a field containing generator %d!\n", gnum); 
  return 0; 
}

int field::generated_by(const pari& elts, pari& exact_elts, int cp, int report)
{
  if (!root_num) 
    root_num = 1; // set it to be the rationals initially

  int i, n = length(elts), redo = 0, res; 
  exact_elts = rvector(n); 

  pari ex; 
  for (i = 0; i<n; i++) {

    res = extend(elts[i], ex, cp, report, i);
    if (!res) {
      clear(); 
      return 0; 
    }

    exact_elts[i] = ex; 
    if (res>1) redo = i; // field changed, so must recompute exact_elts[0]..[i-1]
  }

  if (report && redo > 0) 
    printf("Recomputing exact generator(s) 0 - %d\n", redo-1); 
  for (i=0; i<redo; i++) {
    if (!contains(elts[i], ex, report)) { 
      if (report) printf("Oops: unable to find generator %d!\n", i);
      clear(); 
      return 0;
    }
    exact_elts[i] = ex; 
  }

  return 1;
}

int field::set_field(const pari& poly, const pari& rt, int report)
{
  if (report) printf("Computing number field information (initalg)...\n"); 
  nf = initalg(poly); 
  if (report) {
    printf("Basis for ring of integers:  ");
    nf[6].print(); 
  }

  set_root_indices(); 
  root_num = root_number(rt);
  if (!root_num) {
    printf("Oops, set_field was not given a root of poly.\n"); 
    return 0; 
  }

  /* create numeric integer basis */
  numeric_ibasis = gsubst(nf[6],0,root());
  ibasis_lllmat = lindeppart1(numeric_ibasis, (long)(1+digits*fudge));

  return 1; 
}

int field::canonize(pari& newroot_on_old, int report)
{
  if (is_canonical()) {
    newroot_on_old = polymod(pari(polx[0]), min_poly()); 
    return 1; 
  }

  pari rt1_on_old, newroot_on_rt1; 
  pari c_poly = canonical_poly(min_poly(), rt1_on_old, report); 
  if (!set_field(c_poly, gsubst(lift(rt1_on_old),0,root()), report)) return 0; 
  if (!find_canonical_root(newroot_on_rt1, report)) return 0; 
  newroot_on_old = gsubst(lift(newroot_on_rt1),0,rt1_on_old);
  _is_canonical = true; 

  return 1; 
}


/* The following version of set_field finds the canonical polynomial
   and root for the field in question if cp is non-zero, and the
   reduced root otherwise. */

int field::set_field(pari const& poly, pari const& root, pari& exact_root, int cp, int report)
{
  pari rt2_on_root, rt3_on_rt2, c_poly; 

  if (cp) {

    c_poly = canonical_poly(poly, rt2_on_root, report); 
    if (!set_field(c_poly, gsubst(lift(rt2_on_root),0,root), report)) return 0; 
    if (!find_canonical_root(rt3_on_rt2, report)) return 0;
    exact_root = polymodrecip(gsubst(lift(rt3_on_rt2),0,rt2_on_root)); 
    _is_canonical = true; 

  } else {

    c_poly = reduced_poly(poly, rt2_on_root, report); 
    if (!set_field(c_poly, gsubst(lift(rt2_on_root),0,root), report)) return 0; 
    exact_root = polymodrecip(rt2_on_root); 
    _is_canonical = false; 
  }

  return 1; 
}

int field::find_canonical_root(pari& new_root_on_old, int report)
{
  /* find which root (embedding) gives field of the generator */ 
  if (report) printf("Determining the canonical generating root.\n"); 

  int r1 = nf[1][0].int_value();
  int r2 = nf[1][1].int_value();
  int rn, rnmax; 
  if (numerical_zero(gimag(root()))) {
    rn = 1; /* check real roots */ 
    rnmax = r1; 
  } else {
    rn = r1+1; /* complex roots */ 
    rnmax = r1+r2; 
  }

  for (; rn <= rnmax; rn++) {
    if (contains(root(rn), new_root_on_old, report))
      break;
    if (rnmax!=r1 && contains(root(-rn), new_root_on_old, report)) {
      rn = -rn; 
      break;
    }
  }
  if (rn > rnmax) {
    if (report) printf("Failed to find canonical root generating this field!\n"); 
    return 0; 
  }

  if (report)
    printf("Root number: %d\n", rn); 
  if (root_num==rn) return 1; // Can return quickly.

  root_num = rn; 

  numeric_ibasis = gsubst(nf[6],0,root());
  ibasis_lllmat = lindeppart1(numeric_ibasis, (long)(1+digits*fudge));

  return 1; 
}

pari field::exact_value(const pari& x) const
{
  pari value; 
  contains(x, value); 
  return value; 
}

void field::set_exact_value(pari& x) const 
{ 
  x = exact_value(x); 
}


// relation == 0, if there is no known relation between the fields other 
// than a possibility that this is a subfield of f. 
// relation == 1, if we want to give up as soon as a generator is not in f. 
// relation == 2, if we know this contains f.
// relation == 3, if we know they are the same (but still want to find 
// exact expressions and initialize this.)

int field::generated_by(const pari& field_gens, pari& exact_gens, int cp, const field *f, int relation, int report)
{
  int i, n_gens; 
  pari moos, min_pol, old_gen, gen, ex, rt1_on_gen, rt1, rie; 
  pari f_x_gens; 

  if (field_gens[0].type() != 9) { // Then numeric generators.
    if (!f->contains(field_gens, f_x_gens, report)) {
      if (relation==1 || relation==3 || (relation==2 && 2 * f->degree() > maxdeg)) return 0; 
      return generated_by(field_gens, exact_gens, cp, report); // Do numeric version.
    } else if (relation>=2) { // this == f. 
      goto whole_field; 
    } 
  }
  else f_x_gens = field_gens;
  // might be proper subfield, generated by f_x_gens. 

  root_num = 1; // set field to be the rationals initially

  if (report) {
    printf("Finding subfield generated by:\n"); 
    f_x_gens.print(); 
  }

  n_gens = length(field_gens);
  i = 0; 
  while (lgeff(min_pol)<3 && i < n_gens) {
    if (report) printf("Trying field given by generator %d\n", i); 
    set_moos(f_x_gens[i], moos, min_pol);
    gen = f_x_gens[i]; 
    i++; 
  }

  /* check if the field is still rational */ 
  if (lgeff(min_pol)<3) {
    if (report) printf("Field is rational\n"); 
    for (i=0; i<n_gens; i++) {
      polymod_contains(f_x_gens[i], moos, min_pol, ex); 
      exact_gens[i] = ex; 
    }
    return 1;
  }

  /* check if the field is the same as f */ 
  if (lgeff(min_pol)==lgeff(f->min_poly())) 
    goto whole_field;

  /* extend the field by successive generators */ 
  if (i < n_gens && report) 
    printf("Checking if the other generators are in this subfield\n"); 

  for (;i<n_gens;i++) {

    /* check if the field is the whole of f. */ 
    if (lgeff(min_pol)==lgeff(f->min_poly())) 
      goto whole_field;

    /* check to see if other generators are in this subfield */ 
    if (polymod_contains(f_x_gens[i], moos)) continue; 

    if (report) {
      printf("Subfield does not contain generator %d.\n", i); 
      printf("Trying field generated by this instead.\n"); 
    }
    /* try the new element as generator */ 
    old_gen = gen; 
    gen = f_x_gens[i]; 
    set_moos(gen, moos, min_pol);

    if (polymod_contains(old_gen, moos)) continue; 
  
    /* if that failed, try various combinations of them */ 
    gen = old_gen + f_x_gens[i];
    set_moos(gen, moos, min_pol); 
    if (polymod_contains(old_gen, moos)) continue; 

    gen = pTWO*old_gen + integer(3)*f_x_gens[i]; 
    set_moos(gen, moos, min_pol);
    if (polymod_contains(old_gen, moos)) continue; 

    /* if we get here, we got stuck */ 
    if (report) printf("Failed to find a field containing generator %s\n", i);
    return 0;

  }

  if (cp) { 

    // Find canonical polynomial and root. 
    min_pol = canonical_poly(min_pol, rt1_on_gen, report); 
    rt1 = gsubst(lift(rt1_on_gen),0,f->numeric_value(gen)); 
    if (!set_field(min_pol, rt1, report)) return 0;
    if (!find_canonical_root(ex, report)) return 0;
    _is_canonical = true; 

  } else {

    // Find reduced polynomial. 
    min_pol = reduced_poly(min_pol, rt1_on_gen, report); 
    rt1 = gsubst(lift(rt1_on_gen),0,f->numeric_value(gen)); 
    if (!set_field(min_pol, rt1, report)) return 0;
    _is_canonical = false; 
  }

  /* find exact expressions for generators */ 
  if (!f->contains(root(), rie)) return 0; 
  set_moos(rie, moos, min_pol); 
  exact_gens = rvector(n_gens); 
  for (i=0; i<n_gens; i++) {
    polymod_contains(f_x_gens[i], moos, min_pol, ex); 
    exact_gens[i] = ex; 
  }

  return 1;

 whole_field:

  if (report) printf("Subfield is whole field.\n"); 
  *this = *f; // Make this the same as f. 

  exact_gens = f_x_gens; 

  if (cp && !f->is_canonical()) {
    pari newroot_on_old; 
    if (!canonize(newroot_on_old, report)) return 0; 
    exact_gens = gsubst(f_x_gens, 0, polymodrecip(newroot_on_old)); 
  }

  return 1; 
}

// how=0 condensed one-line format, human readable
// how=1 expanded multi-line format, human readable
// how=2 pari nf component, pari readable, omits root information. 
// how=3 minimum essential information, [min-poly, root], pari 
//       readable; min-poly * -1 (i.e. lead coeff==-1) if root is
//       not canonical. 

void field::print(int how) const
{

  if (!is_set() && how < 2) {
    fprintf(outfile, "rationals");
    if (how == 1) printf("\n"); 
    return;
  }

  if (how == 0) {
    nf[0].print(no_newline); // Minimum polynomial. 
    fprintf(outfile, " "); 
    nf[1].print(no_newline); // Signature. 
    fprintf(outfile, " "); 
    nf[2].print(no_newline); // Discriminant. 
    fprintf(outfile, " R(%d) = ", root_num); 
    root().print(pp_tiny); 
  }
  if (how == 1) {
    fprintf(outfile, "%sminumum polynomial: ", _is_canonical ? "canonical " : "");
    nf[0].print(); 
    fprintf(outfile, "%sroot number: %d\n", _is_canonical ? "canonical " : "", root_num);
    fprintf(outfile, "numerical value of root: ");
    root().print(pp_short);
    fprintf(outfile, "signature: ");
    nf[1].print();
    fprintf(outfile, "discriminant: ");
    nf[2].print();
    fprintf(outfile, "factors: ");
    factor(nf[2]).print();
    fprintf(outfile, "index of ring of integers: ");
    nf[3].print();
    fprintf(outfile, "basis for ring of integers: ");
    nf[6].print();
  }
  if (how == 2) {
    nf.print();
  }
  if (how == 3) { 
    // Print so it can be read as a single pari expression. 
    // If not canonical, save -1 times min_poly. 
    fprintf(outfile, "[");
    (nf[0] * (_is_canonical ? pONE : -pONE)).print(no_newline); 
    fprintf(outfile, ", %d]\n", root_num); 
  }
}   

void named_field::print(int how) const 
{
  if (how==0)
    printf("%s: ", _name.c_str()); 
  else if (how==1)
    printf("%s\n", _name.c_str()); 

  field::print(how); 
}
