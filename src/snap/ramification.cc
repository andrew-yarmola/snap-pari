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
#include <sys/times.h>
#include "pariwrap.hh"
#include <iostream>

static pari solve_hilbert_equation(const pari& nf, pari a, pari b, const pari& ideal, int idealpow, int report, int abort_after);

static int ramified(pari hil[2], pari ibasis_hil[2], pari itf, pari prime, int report, int abort_after); 

static int next_ideal_rep(pari& rep, const pari& diag);

static pari lift_solution(pari x, pari y, pari z, int ipow, int to_pow, const pari& ideal, const pari& ideal2, const pari& nf, pari unif_power, const pari& a, const pari& b, int which_one, int e, int report, clock_t abort_time);


pari ramification(const pari hilbert[2], pari itf, int int_traces, int report, int abort_after)
{
  pari fl, pl, n, pr;
  pari hil[2];
  pari ibasis_hil[2]; 
  int i, imax, j; 

  if (hilbert[1]==pZERO) {
    return rvector(0); 
  }

  // int old_gl = garbage_limit; 
  // garbage_limit = 100; 

  hil[0] = hilbert[0];
  hil[1] = hilbert[1];

  if (!int_traces) {
    /* ensure that a,b are integers by multiplying by a suitable square */ 
    for (j=0; j<2; j++) {
      n = pONE; 
      fl = factor(denom(algtobasis(itf,hil[j])));
      imax = matsize(fl)[0].int_value(); 
      for (i=0; i<imax; i++)
	n *= pow(fl[0][i],pTWO*ceil(fl[1][i]/pTWO));
      hil[j] *= n; 
    }

    /* pl[0] is our list of primes; all primes dividing a,b or 2 are candidates
       when the traces are not integral */ 
    pl = idealfactor(itf, pTWO*hil[0]*hil[1]);

  } else {

    /* pl[0] is the list of primes; when traces are integral we only need to 
       consider primes dividing b */ 
    pl = idealfactor(itf, hil[1]);

  }

  /* now convert a,b into expressions on the integer basis */ 
  ibasis_hil[0] = algtobasis(itf,hil[0]); 
  ibasis_hil[1] = algtobasis(itf,hil[1]);

  int n_ramified = 0; 
  imax = length(pl[0]); 
  for (i=0; i<imax; i++) {

    pr = pl[0][i];
    switch (ramified(hil,ibasis_hil,itf,pr,report,abort_after)) {
    case 1:
      n_ramified++; 
      pl[1][i] = pONE; 
      break;
    case -1:
      warn("ramification test aborted\n");
      n_ramified++; 
      pl[1][i] = -pONE; 
      break;
    default:
      pl[1][i] = pZERO; 
    }
  }

  /* copy the primes where ramification occurs out to a single list */ 
  pari ram = rvector(n_ramified); 
  j = 0; 
  for (i=0; i<n_ramified; i++) {
    while (pl[1][j]==pZERO) j++; 
    ram[i] = pl[0][j];
    if (pl[1][j]==-pONE) ram[i][0] = -ram[i][0]; /* flag any aborted cases */ 
    j++;
  }

  // garbage_limit = old_gl;

  return ram; 
}

int ramified(pari hil[2], pari ibasis_hil[2], pari itf, pari prime, int report, int abort_after)
{

  /* First check if p divides 2 or the index of the ring of integers.
     these cases are not handled by the code below so we use
     solve_hilbert_equation instead. If solve_hilbert_equation can`t
     solve ax^2 + by^2 == 1 it returns pZERO which has type 1. this
     indicates that the qt algebra ramifies at this prime. the reason
     that we do not always use solve_hilbert_equation is that it is 
     likely to be slower than what follows. 

     When p divides the index of the ring of integers we can`t reduce
     integers modulo the prime ideal simply by reducing mod p and modulo
     the second generator of the ideal, because the prime may appear
     in the denominators of an algebraic integer expressed wrt the 
     primitive element. */

  if (prime[0]==pTWO || gmod(itf[3],prime[0])==pZERO) {
    pari res = solve_hilbert_equation(itf,ibasis_hil[0],ibasis_hil[1],prime,0,
				      report,abort_after); 
    if (res.type()==1) {
      if (res == pZERO) 
	return 1; 
      else 
	return -1;
    } 
    else return 0; 
  }

  /* get valuations of a and b */ 
  int va = idealval(itf,hil[0],prime);
  int vb = idealval(itf,hil[1],prime);

  /* no ramification if both have even valuation (zero after removing squares) */ 
  if ((va%2==0) && (vb%2==0)) 
    return 0; 

  /* p_inv has valuation -1 at p, >=0 elsewhere */ 
  pari p_inv = polymod((itf[6] * prime[4])/prime[0],itf[0]); 

  /* w = a/p^v(a), b/p^v(b) or -ab/p^v(ab), choosing the one for
     which the power of p in the denominator is even. then the qt
     algebra ramifies iff w is not a quadratic residue in the
     completion of the algebraic integers at p. */
  pari w; 
  if ((va%2==0) && (vb%2!=0)) { 
    w = hil[0] * pow(p_inv, integer(va)); 
  } else if ((va%2!=0) && (vb%2==0)) {
    w = hil[1] * pow(p_inv, integer(vb)); 
  } else {
    w = -hil[0] * hil[1] * pow(p_inv,integer(va+vb));
  }

  /* resf is min poly of the finite residue field. */
  pari resf = itf[6] * gmod(prime[1],prime[0]);
  if (resf==pZERO) resf = itf[0];

  /* by Hensel's Lemma w is a quadratic residue in the 
     completion of the field at p iff it is a quadratic 
     residue in the finite residue field mod p. 
     qt algebra is ramified at p iff w is NOT a quadratic
     residue. */ 

  pari fac, c, qp; 

  if (lgeff(resf)-1 < 2) {

    /* when degree of min poly is 1 we are just in the finite 
       field with prime[0] elements */ 

    w = gsubst(w[1],0,-resf[0]/resf[1]);
    fac = factmod(pow(polx[0],pTWO) - w, prime[0]);

  } else {

    /* change x into y in w and resf: we have to use y instead of x
       because factmod9 requires that the variable in the (polymod)
       coefficients of the polynomial have a lower priority than the
       variable of the polynomial itself. */

    w = gsubst(lift(w),0,polx[1]); 
    c = content(w); 

    qp = ((c.type()==t_FRAC) ? c[1]:pONE)*(pow(polx[0],pTWO)-w);

    resf = gsubst(resf,0,polx[1]); 
#if 0
    std::cout << "about to call factmod9\n";
    std::cout << "poly = "; qp.print(); 
    std::cout << "p = "; prime[0].print(); 
    std::cout << "resf = "; resf.print();
#endif
    fac = factmod9(qp, prime[0], resf); 
  }
  int is_qr = (length(fac[0]) > 1 || fac[1][0] > pONE); 

  return !is_qr; 
}

static clock_t clock_time()
{
  struct tms buf;
  times(&buf); 
  return (buf.tms_utime + buf.tms_stime); 
}

pari solve_hilbert_equation(const pari& nf, pari a, pari b, const pari& ideal2, 
			    int to_pow, int report, int abort_after)
{
  int n = length(ideal2[1]); 
  pari element_zero = cvector(n); 
  pari element_one = cvector(n); 
  element_one[0] = pONE;

  /* reduce valuation of a, b to at most 1 while keeping a, b integral. */
  int val, fix_val, odd_val = 0; 
  val = element_val(nf,a,ideal2); 
  fix_val = val - ((val%2)!=0); /* val - fix_val == 0 or 1 */ 
  if (fix_val) a = element_mul(nf,a,element_pow(nf,ideal2[4]/ideal2[0],integer(fix_val))); 
  odd_val += (val%2); 

  val = element_val(nf,b,ideal2); 
  fix_val = val - ((val%2)!=0); /* val - fix_val == 0 or 1 */ 
  if (fix_val) b = element_mul(nf,b,element_pow(nf,ideal2[4]/ideal2[0],integer(fix_val))); 
  odd_val += (val%2); 
  
  if (report) {
    printf("after removing squares we have\n");
    pari aa = basistoalg(nf,a);
    pari bb = basistoalg(nf,b);

    printf("a: "); aa.print();
    printf("b: "); bb.print();

    printf("factorization of a:\n"); 
    idealfactor(nf, aa).print(); 
    printf("factorization of b:\n"); 
    idealfactor(nf, bb).print(); 
  }

  int dyadic_ramification = 0; 
  if (ideal2[0]==pTWO) {
    dyadic_ramification = ideal2[2].int_value(); 
  }

    /* in general Hensel lifting says that if we can solve an equation
       in the residue field we can lift it to the completion. in the
       case of a dyadic prime we have to solve it in the algebraic
       integers mod the prime to the power 2*e+3 where e is the
       ramification index of the ideal (in the extension of the
       rationals to the field). */

  if (to_pow==0) {
    if (dyadic_ramification) 
      to_pow = 2 * dyadic_ramification + 3; 
    else if (odd_val) 
      to_pow = 3;
    else 
      to_pow = 1; 
  }

  clock_t abort_time = 0; 
  if (abort_after)
    abort_time = clock_time() + abort_after * 60;

  pari ideal = idealhermite(nf, ideal2);    

  pari sol; 
  sol = lift_solution(element_zero, element_zero, element_one, 0, to_pow, ideal, ideal2, 
		      nf, element_one, a, b, 0, dyadic_ramification, report, abort_time); 
  if (sol.type()!=1) return sol; 
  sol = lift_solution(element_one, element_zero, element_zero, 0, to_pow, ideal, ideal2, 
		      nf, element_one, a, b, 1, dyadic_ramification, report, abort_time); 
  if (sol.type()!=1) return sol; 
  sol = lift_solution(element_zero, element_one, element_zero, 0, to_pow, ideal, ideal2, 
		      nf, element_one, a, b, 2, dyadic_ramification, report, abort_time); 
  return sol; 
}

/* unif_pow is a power of a uniformizing element, initially just 1. 
   suppose we have a "solution" of the equation for elements (x,y,z) 
   in a set of representitives for Z/p^ipow. lift_solution looks at 
   representitives for Z/p^(ipow+1) congruent to (x,y,z) mod p^ipow
   and tries to find a solution here. A solution is (x,y,z) such that 
   ax^2 + by^2 - z^2 has valuation ipow+1 if p is not dyadic, a little 
   larger if p is dyadic. */ 

static pari lift_solution(pari x, pari y, pari z, int ipow, int to_pow, const pari& ideal,
			  const pari& ideal2, const pari& nf, pari unif_power, 
			  const pari& a, const pari& b, int which_one, int e, int report,
			  clock_t abort_time)
{
  ipow++;

  int solpow;
  if (e) {

    /* When the ideal is dyadic, the class of ax^2 + by^2 - z^2 
       mod ideal^min(2*ipow,ipow+e) depends only on the class of 
       x,y,z mod ideal^ipow. Thus we can look for solutions modulo
       this larger power of the ideal. */
 
    solpow = (ipow<e) ? 2*ipow : ipow+e;
  } else {
    solpow = ipow; 
  }

  pari x1, y1, z1, u, v, res;

  int n = length(unif_power); 
  u = cvector(n);
  v = cvector(n); 

  switch (which_one) {
  case 0: /* z=1 */ 
    z1 = z; 
    break;
  case 1: /* x=1 */ 
    x1 = x; 
    break;
  default: /* y=1 */ 
    y1 = y; 
    break;
  }


  int i, val; 
  pari elt, e1, sol = rvector(3); 

  if (report) {
    for (i=0; i<ipow; i++) printf(" "); 
    if (ipow > 1) {
      printf("checking classes of Z/p^%d == [above] mod p^%d, for solution mod p^%d\n",
	     ipow,ipow-1,solpow); 
    } else {
      printf("checking classes of Z/p for solution mod p^%d\n", solpow); 
    }
  }

  while (1) {


    /* existing_sol_ideal is the ideal mod which we have an existing solution x,y,z.
       u and v vary congruence mod ideal. thus x1,y1,z1 vary their congruence classes 
       mod this_ideal while staying congruent to x,y,z mod existing_sol_ideal. */ 

    switch (which_one) {
    case 0: /* z=1 */ 
      x1 = element_mul(nf,u,unif_power) + x;
      y1 = element_mul(nf,v,unif_power) + y; 
      break;
    case 1: /* x=1 */ 
      y1 = element_mul(nf,u,unif_power) + y;
      z1 = element_mul(nf,v,unif_power) + z; 
      break;
    default: /* y=1 */  
      x1 = element_mul(nf,u,unif_power) + x;
      z1 = element_mul(nf,v,unif_power) + z; 
      break;
    }

    if (gc_off > 0) {
      std::cout << "GARBAGE collection off!!\n";
      exit(0); 
    }
    if (garbage_count > garbage_limit) 
      pari::collect_garbage(); 

    if (report) {
      for (i=0; i<ipow; i++) printf(" "); 
      printf("trying "); 
      sol[0] = x1; 
      sol[1] = y1; 
      sol[2] = z1; 
      sol.print(); 
    }

    /* Want a solution of the hilbert equation mod p^solpow. */ 

    elt = element_mul(nf,a, element_mul(nf,x1,x1)) + 
      element_mul(nf,b, element_mul(nf,y1,y1)) - 
	element_mul(nf,z1,z1);

    if (gnorml2(elt)!=pZERO)
      val = element_val(nf, elt, ideal2); 
    else 
      val = to_pow; // We have an exact solution. 

    if (val >= to_pow) {
      sol[0] = x1; 
      sol[1] = y1; 
      sol[2] = z1; 
      if (report) {
	for (i=0; i<ipow; i++) printf(" "); 
	printf("solution mod p^%d = ", val); 
	sol.print(); 
      }
      return sol; 
    }

    if (val >= solpow) {

      if (report) {
	sol[0] = x1; 
	sol[1] = y1; 
	sol[2] = z1; 
	for (i=0; i<ipow; i++) printf(" "); 
	printf("this is a solution mod p^%d\n", val); 
      }

      e1 = element_mul(nf, unif_power, ideal2[1]);
      res = lift_solution(x1,y1,z1,ipow,to_pow,ideal,ideal2,nf,e1,
			  a,b,which_one,e,report,abort_time);
      if (res.type()!=1 || res == pONE) return res;
    }

    /* try next congruence mod this_ideal */ 

    if (!next_ideal_rep(u,ideal)) {
      if (!next_ideal_rep(v,ideal)) return pZERO; 
    }

    if (abort_time && abort_time < clock_time()) {
      return pONE; 
    }
  }
}



/* return 0 if rep goes back to [0,0,..,0]~, else 1 */ 

static int next_ideal_rep(pari& rep, const pari& ideal)
{
  int i = 0, n = length(rep); 

  while (1) {
    rep[i] += pONE;
    if (rep[i] < ideal[i][i]) break; 
    rep[i] = pZERO; 
    i++; 
    if (i==n) break; 
  }
  return (i<n); 
}

