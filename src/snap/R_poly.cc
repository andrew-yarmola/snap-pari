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
#include "R_poly.hh"

#ifdef __DECCXX
#define exception math_exeption
#include <math.h>
#undef exception
#else
#include <math.h>
#endif

#include <float.h>

static bool newton_ok; 

// Treat vector doubles as polynomials. 

static double abs(double x)
{
  return x > 0 ? x : -x; 
}

R_poly operator * (R_poly const& a, R_poly const& b)
{
  int amax = a.size(), bmax = b.size(); 
  R_poly res(amax+bmax-1, 0.0); 
  int i, j;

  for (i=0; i<amax; i++)
    for (j=0; j<bmax; j++)
      res[i+j] += a[i]*b[j]; 

  return res; 
}

R_poly sqr(R_poly const& f)
{
  int sz = f.size();
  R_poly res(2*sz-1, 0.0); 
  int i, j;

  for (i=0; i<sz; i++)
    for (j=0; j<sz; j++)
      res[i+j] += f[i]*f[j]; 

  return res; 
}

R_poly times_x_pow(R_poly const& f, int n)
{
  R_poly g;
  if (f.size() + n <= 0) return g; 
  g.resize(f.size() + n); 
  int i; 
  for (i=0; i<n; i++)
    g[i] = 0.0; 
  for (; i<g.size(); i++)
    g[i] = f[i-n]; 
  return g;
}

int x_pow(R_poly const& f)
{
  int i; 
  for (i=0; i<f.size(); i++)
    if (f[i]!=0.0) break; 
  return i; 
}

R_poly sqrt(R_poly const& f, double eps)
{
  R_poly z; 
  if (f.size()%2 != 1) return z; 

  int xp = x_pow(f); 
  if (xp%2 != 0) return z; 

  R_poly g = times_x_pow(f, -xp); 

  if (g[0] < 0.0) return z; 

  R_poly r((g.size()+1)/2);

  r[0] = sqrt(g[0]); 

  // Compute coeffs of r. 
  double sn; 
  int i, n; 
  for (n=1; n<r.size(); n++) {
    sn = 0.0;
    for (i=1; 2*i < n; i++) sn += 2.0*r[i]*r[n-i];
    if (2*i==n) sn += r[i]*r[i];
    r[n] = (g[n] - sn)/(2.0*r[0]);
  }
  // Check if remaining coeffs of g tally with square of r. 
  for (; n<g.size(); n++) {
    sn = 0.0;
    for (i = n-r.size()+1; 2*i < n; i++) sn += 2.0*r[i]*r[n-i];
    if (2*i==n) sn += r[i]*r[i];
    if (abs(g[n] - sn) > eps) return z; 
  }

  if (!xp) return r;
  return times_x_pow(r, xp/2); 
}

R_poly operator + (R_poly const& a, R_poly const& b)
{
  int i, nmin, nmax; 
  if (a.size() < b.size()) {
    nmin = a.size(); 
    nmax = b.size();
  } else {
    nmin = b.size(); 
    nmax = a.size();
  }
  R_poly res(nmax); 
  for (i=0; i<nmin; i++) 
    res[i] = a[i]+b[i]; 
  if (a.size() < b.size())
    for (; i<nmax; i++) res[i] = b[i]; 
  else 
    for (; i<nmax; i++) res[i] = a[i]; 
  return res; 
}

R_poly operator - (R_poly const& a, R_poly const& b)
{
  int i, nmin, nmax; 
  if (a.size() < b.size()) {
    nmin = a.size(); 
    nmax = b.size();
  } else {
    nmin = b.size(); 
    nmax = a.size();
  }
  R_poly res(nmax); 
  for (i=0; i<nmin; i++) 
    res[i] = a[i]-b[i]; 
  if (a.size() < b.size())
    for (; i<nmax; i++) res[i] = -b[i]; 
  else 
    for (; i<nmax; i++) res[i] = a[i]; 
  return res; 
}

R_poly operator * (double r, R_poly const& f)
{
  R_poly res(f.size()); 
  int i;
  for (i=0; i<res.size(); i++)
    res[i]=r*f[i];
  return res; 
}

R_poly operator / (R_poly const& f, double r)
{
  R_poly res(f.size()); 
  int i;
  for (i=0; i<res.size(); i++)
    res[i]=f[i]/r;
  return res; 
}

R_poly& operator *= (R_poly& f, double r)
{
  int i;
  for (i=0; i<f.size(); i++)
    f[i] *= r;
  return f; 
}

R_poly& operator /= (R_poly& f, double r)
{
  int i;
  for (i=0; i<f.size(); i++)
    f[i] /= r;
  return f; 
}

R_poly ediv(R_poly const& n, R_poly d, R_poly& r, double eps)
{
  r = n;
  if (!d.size() || d.size() > n.size()) {
    return R_poly();
  }

  int i, j, dd = d.size()-1; 
  double lead = d[dd];
  d /= lead; 

  R_poly q(r.size()-dd); 
  for (i = q.size()-1; i>=0; i--) {
    q[i] = r[i+dd];
    for (j=0; j<d.size(); j++) r[i+j] -= q[i]*d[j];
    q[i] /= lead; 
  }

  r.resize(d.size()-1);
  chop(r, eps);

  return q; 
}

void rmod(R_poly& n, R_poly d, double eps)
{
  if (!d.size() || d.size() > n.size()) {
    return;
  }

  int i, j, dd = d.size()-1; 
  double lead = d[dd], qi;
  d /= lead; 

  for (i = n.size()-d.size(); i>=0; i--) {
    qi = n[i+dd];
    for (j=0; j<d.size(); j++) n[i+j] -= qi*d[j];
  }

  n.resize(d.size()-1);
  chop(n, eps);
}


R_poly hcf(R_poly f, R_poly g, double eps)
{
  chop(f, eps);
  chop(g, eps); 
  if (g.size() > f.size()) swap(f, g); 

  if (!f.size()) return f; 

  while (g.size()) {

    rmod(f, g, eps);
    swap(f, g);

  }
  
  f /= f[f.size()-1];
  return f; 
}

int order(R_poly f, double x, double eps)
{
  int ord=0;
  R_poly g(2), r;
  g[1] = 1;
  g[0] = -x;
  while (f.size()) {
    f = ediv(f,g,r,eps);
    if (r.size()) break;
    ++ord;
  }
  return ord;
}

double evaluate(R_poly const& p, double x)
{
  if (!p.size()) return 0.0;
  R_poly::const_iterator it = p.end()-1, st = p.begin(); 
  register double sum = *it;
  while (it!=st) {
    it--;
    sum *= x; 
    sum += *it; 
  } 
  return sum; 
}

int value_exponent(R_poly const& p, double x)
{
  int ex, epx = 0, ecp, emax = -1000; 
  if (!p.size()) return emax;
  R_poly::const_iterator e = p.end(), i = p.begin(); 
  if (frexp(x, &ex)==0.) {
    if (frexp(p[0],&ecp)==0.) return emax; 
    return ecp; 
  }
  for (; i!=e; i++) {
    if (frexp(*i, &ecp)!=0. && ecp + epx > emax) 
      emax = ecp + epx; 
    epx += ex; 
  }
  return emax; 
}

R_poly derivative(R_poly const& p)
{
  int i, n = (p.size() > 0) ? p.size()-1 : 0; 
  R_poly res(n); 
  for (i=1; i<=n; i++)
    res[i-1] = double(i) * p[i]; 
  return res; 
}

ostream& operator << (ostream& out, R_poly const& v)
{
  int i, n=v.size(); 
  out << "[";
  for (i=0; i<n; i++) {
    if (i>0) out << ", "; 
    out << v[i];
  }
  return out << "]";
}

double upper_bound(R_poly const& f)
{
  int i, n=f.size();
  double sum = 0.0;
  for (i=0; i<n; i++) sum += abs(f[i]); 
  return sum; 
}

void chop(R_poly& f, double eps)
{
  // set to zero poly if all coeffs are small. 
  double ub = upper_bound(f); 
  if (ub < eps * (f.size()+1)) {
    f.resize(0); 
    return; 
  }

  // find leading nonzero coeff. 
  int i;
  for (i=f.size()-1; i>=0; i--) {
    if (abs(f[i]/ub) > eps) break; 
  }

  // resize f accordingly. 
  if (i+1 < f.size()) f.resize(i+1); 

  // reset small coeffs to zero. 
  for (; i>=0; i--)
    if (abs(f[i]/ub) < eps) f[i] = 0.0; 
}


double newton_find_root(double xlo, double xhi, R_poly const& f, R_poly const& df)
{
  // Polish the roots using newton's method. 
  int k; 
  double adj, approx = (xlo + xhi)/2.0, prev_adj = xhi - xlo; 
  if (xlo >= xhi) return approx;
  bool xlo_neg = (evaluate(f, xlo) < 0.0); 
  bool ok = false; 
  if (xlo < 0.0 && xhi > 0.0 && f[0]==0.0) 
    return 0.0; 
  double xlo_start = xlo, xhi_start = xhi; 
  int do_binary = 0, newton_done = 0; 
  double value;

  int eps_exp = value_exponent(f, (abs(xlo) > abs(xhi) ? xlo : xhi)) - DBL_MANT_DIG + 3;
  double eps = ldexp(1.0, eps_exp); 

  for (k = 100; k > 0; k--) {

    // try one newton step. 
    value = evaluate(f, approx);
    adj = value/evaluate(df, approx); 
    approx -= adj; 
    newton_done++; 

    // see if we're home. 
    if (abs(value) < eps && k > 3) { 
      ok = true; 
      if (adj==0.0) break; 
      k = 3;
    }

    // if we step outside the interval, or convergence is poor, 
    // switch to a binary search instead for a few steps. 
    if (approx < xlo || approx > xhi ||
	!ok && newton_done > 4 && abs(adj/prev_adj) > 0.7) {

      if (approx <= xlo || approx >= xhi)
	approx = (xlo + xhi)/2.0;

      for (do_binary = 10; do_binary > 0; do_binary--) {
	if (xlo_neg == (evaluate(f, approx) < 0.0)) xlo = approx; 
	else xhi = approx; 
	approx = (xlo + xhi)/2.0; 
      }

      newton_done = 0; 
    }

    prev_adj = adj; 
  }

  if (!ok) newton_ok = false; 
#if 0
  if (abs(value) > eps) {
    cout << "inaccurate root " << approx << " in [" << xlo_start << ", " 
	 << xhi_start << "]" << endl; 
    cout << "poly was " << f << endl; 
  }
#endif

  return approx; 
}

R_poly poly_from_roots(vector<double> const& rts)
{
  R_poly lin(2); 
  lin[1] = 1.0; 
  R_poly res(1,1.0); 
  int i, n=rts.size(); 
  for (i=0; i<n; i++) {
    lin[0] = -rts[i]; 
    res = res * lin; 
  }
  return res; 
}

int signof(double x, double eps)
{
  if (x < -eps) return -1; 
  else if (x > eps) return 1; 
  return 0;
}

vector<double> roots_in_monotone_intervals(R_poly const& f, R_poly const& df, vector<double> stat, double eps)
{
  vector<double> roots; 

  int i;
  double x, y, prev_x; 
  int sign, prev_sign = 0; 

  for (i=0; i<stat.size(); i++) {
    x = stat[i]; 
    y = evaluate(f, x);
    sign = signof(y, eps); 

    if (sign==0 && y != 0.) { 
      // This indicates a repeated root. In certain cases we want to 
      // adjust this to the actual sign of y, rather than the sign of 
      // y to within eps. 
      if ((y < 0.) == (prev_sign > 0)) { 
	int ey, ev; 
	frexp(y, &ey); 
	ev = value_exponent(f, x); 
	if (ey > (ev - DBL_MANT_DIG + 4)) 
	  sign = -prev_sign; // ie. sign = sign(y); 
      }
    }

    if (sign==0) { // x is a root. 
      roots.push_back(x); 
    } else if (sign * prev_sign == -1) { // sign change. 
      roots.push_back(newton_find_root(prev_x,x,f,df)); 
    }
    prev_x = x; 
    prev_sign = sign; 
  }
  return roots; 
}

vector<double> roots(R_poly const& f, double eps)
{
  return all_deriv_roots(f, eps)[0]; 
}

vector<vector<double> > all_deriv_roots(R_poly f, double eps)
{
  vector<vector<double> > roots(1); 

  // check if f==0. 
  double ub = upper_bound(f); 
  if (abs(ub) < eps) {
    return roots; 
  }

  // remove any leading zeros. 
  chop(f, eps); 

  // get out quick for constant and linear cases. 
  if (f.size() < 3) { 
    if (f.size()== 2)
      roots[0].push_back(-f[0]/f[1]);
    return roots; 
  }

  // make polynomial monic. 
  int i; 
  double lead = f[f.size()-1]; 
  for (i=0; i<f.size(); i++)
    f[i] /= lead; 

  // compute all derivatives of f down to linear. 
  vector<R_poly> fD(f.size()-1);
  R_poly df = f; 
  for (i=0; i<fD.size(); i++) {
    fD[i] = df; 
    df = derivative(df); 
  }

  roots.resize(f.size()-1);

  i = fD.size()-1; 
  df = fD[i]; 
  roots[i].push_back(-df[0]/df[1]); 
  i--;

  newton_ok = true; 

  int big_neg_sign, big_pos_sign = 1;
  double xlo, xhi; 

  vector<double> intervals; 
  for (; i>=0; i--) {
#if 0
    cout << " roots " << roots[i+1] << endl;
    cout << " poly " << fD[i] << endl;
#endif

    intervals = roots[i+1]; 

    // If polynomial (fD[i]) changes sign before its first stationary
    // value, find a point where it takes its eventual sign on large
    // negative values. Similarly, find a point where it takes its eventual
    // sign for large positive values if this is not the sign at the 
    // last stationary value. 
    
    big_neg_sign = ((fD[i].size()-1) & 1) ? -1 : 1; 
    if (!intervals.size()) { 
      if (big_neg_sign != big_pos_sign) {
	xlo = -1.0;
	while (signof(evaluate(fD[i],xlo), eps)!=big_neg_sign) xlo *= 2; 
	xhi = 1.0;
	while (signof(evaluate(fD[i],xlo), eps)!=big_pos_sign) xhi *= 2; 
	intervals.push_back(xlo);
	intervals.push_back(xhi); 
      }
    } else {
      if (signof(evaluate(fD[i],intervals.front()), eps) * big_neg_sign == -1) {
	xlo = intervals.front()-1.0;
	while (signof(evaluate(fD[i],xlo), eps)!=big_neg_sign) xlo -= (intervals.front()-xlo);
	intervals.insert(intervals.begin(), xlo);
      }
      if (signof(evaluate(fD[i],intervals.back()), eps) * big_pos_sign == -1) {
	xhi = intervals.back()+1.0; 
	while (signof(evaluate(fD[i],xhi), eps)!=big_pos_sign) xhi += (xhi-intervals.back()); 
	intervals.push_back(xhi); 
      }
    }

    roots[i] = roots_in_monotone_intervals(fD[i], fD[i+1], intervals, eps); 
  }

  if (!newton_ok) {
    cout << "Root finding failed for poly: " << f << endl; 
  }

  return roots; 
}

vector<double> roots(R_poly const& f, double xlo, double xhi, double eps)
{
  return all_deriv_roots(f, xlo, xhi, eps)[0];
}

vector<vector<double> > all_deriv_roots(R_poly f, double xlo, double xhi, double eps)
{
  vector<vector<double> > roots(1); 

  // Dispose of trivial cases. 
  if (xlo > xhi) return roots; 

  double ub = upper_bound(f); 
  if (abs(ub) < eps) return roots; // return empty list if f = 0. 

  // remove any leading zeros. 
  chop(f, eps); 

  // Get out quick for constant and linear cases. 
  if (f.size() < 2) return roots; 
  double rt; 
  if (f.size()== 2) {
    rt = -f[0]/f[1];
    if (xlo < rt && rt < xhi)
      roots[0].push_back(rt);
    return roots; 
  }

  // Make polynomial monic. 
  int i; 
  double lead = f[f.size()-1]; 
  for (i=0; i<f.size(); i++)
    f[i] /= lead; 

  // Compute all derivatives of f down to linear. 
  vector<R_poly> fD(f.size()-1);
  R_poly df = f; 
  for (i=0; i<fD.size(); i++) {
    fD[i] = df; 
    df = derivative(df); 
  }

  roots.resize(f.size() - 1); 

  double xlolo = xlo - 0.05*(xhi-xlo); 
  double xhihi = xhi + 0.05*(xhi-xlo); 
  i = fD.size()-1; 
  df = fD[i]; 
  rt = -df[0]/df[1];
  if (xlolo < rt && rt < xhihi)
    roots[i].push_back(rt); 
  i--;

  newton_ok = true; 

  vector<double> intervals; 
  for (; i>=0; i--) {
#if 0
    cout << " roots " << roots[i+1] << endl;
    cout << " poly " << fD[i] << endl;
#endif

    intervals = roots[i+1]; 

    // Add the endpooints of the interval to roots if it
    // does not already contain them. 
    if (intervals.size()) {
      if (xlolo < intervals.front())
	intervals.insert(intervals.begin(), xlolo);
      if (xhihi > intervals.back())
	intervals.push_back(xhihi); 
    } else {
      intervals.push_back(xlolo);
      intervals.push_back(xhihi);
    }

    roots[i] = roots_in_monotone_intervals(fD[i], fD[i+1], intervals, eps); 
  }

  if (!newton_ok) {
    cout << "Root finding failed for poly: " << f << endl; 
  }

  // Remove any roots which fall outside the interval [xlo, xhi]. 
  vector<double>::iterator j; 
  for (i=0; i<roots.size(); i++) {
    while (roots[i].size() && roots[i][0] < xlo) roots[i].erase(roots[i].begin());
    j = roots[i].end(); j--;
    while (roots[i].size() && *j > xhi) roots[i].erase(j--);
  }

  return roots; 
}



#if WMAIN
#include <iomanip>

main()
{
  double eps = 1e-10; 
  R_poly a(3), b(5); 

  a[0] = 1.0;
  a[1] = 0.0;
  a[2] = 1.0;

  b[0] = 1.0;
  b[1] = 2.0;
  b[2] = -1.0;
  b[3] = -2.0;
  b[4] = 1.0;

  cout << "Poly a: " << a << endl; 
  cout << "Poly b: " << b << endl; 
  cout << "Poly a^2: " << sqr(a) << endl; 

  cout << "Poly a*b: " << (a*b) << endl;   
  cout << "Poly 3*a: " << (3.0*a) << endl;   
  cout << "Poly b-a: " << (b-a) << endl;   

  cout << "Deriv b: " << derivative(b) << endl; 
  cout << "Value b(0.5): " << evaluate(b,0.5) << endl << endl; 

  R_poly c = b-a; 

  cout << "Poly c: " << c << endl; 
  R_poly d = sqr(c); 
  cout << "Square of c: " << d << endl; 
  cout << "Square root: " << sqrt(d, 1e-10) << endl << endl; 

  d = d+a;
  cout << "Poly c^2 + a: " << d << endl;
  R_poly q, r;
  q = ediv(d, 2*c, r); 
  cout << "Quotient of div by 2*c: " << q << endl;
  cout << "Remainder of div by c: " << r << endl;

  d = r; 
  q = ediv(d, a, r); 
  cout << "Quotient of div by a: " << q << endl;
  cout << "Remainder of div by a: " << r << endl << endl;

  vector<double> rts(5); 
  rts[0] = -1.3;
  rts[1] = -.7;
  rts[2] = -.4;
  rts[3] = 0.3;
  rts[4] = 1.6;

  a = poly_from_roots(rts); 
  cout << "Poly a:  " << a << endl; 
  cout << "Roots: " << rts << endl; 

  rts.resize(4);
  rts[0] = -1.3;
  rts[1] = -.7;
  rts[2] = -.7;
  rts[3] = 0.5;

  b = poly_from_roots(rts); 
  cout << "Poly b:  " << b << endl; 
  cout << "Roots: " << rts << endl; 

  c = hcf(a, b, eps); 
  cout << "Poly c = hcf(a,b): " << c << endl; 
  cout << "Roots of c: " << roots(c, eps) << endl << endl; 

  rts.resize(5);
  rts[0] = -1.3;
  rts[1] = -.7;
  rts[2] = -.4;
  rts[3] = 0.3;
  rts[4] = 1.6;

  cout << "Roots: " << rts << endl; 
  R_poly f = poly_from_roots(rts); 
  cout << "Poly:  " << f << endl; 
  cout << "Roots: " << roots(f,eps) << endl; 
  cout << "Roots in [-1,1]: " << roots(f, -1.0, 1.0, eps) << endl << endl; 

  rts.resize(4);
  rts[0] = -1.3;
  rts[1] = -.7;
  rts[2] = -.7;
  rts[3] = 0.3;

  cout << "Roots: " << rts << endl; 
  f = poly_from_roots(rts); 
  cout << "Poly:  " << f << endl; 

  cout << "Roots: " << roots(f,eps) << endl; 
  cout << "Roots in [-1, 1]: " << roots(f,-1.,1.,eps) << endl; 

  vector<vector<double> > adr = all_deriv_roots(f,eps);
  cout << "All deriv roots: ";
  int i;
  for (i=0; i<adr.size(); i++) 
    cout << adr[i] << ' ';
  cout << endl; 

  adr = all_deriv_roots(f,-1.,1.,eps);
  cout << "All deriv roots in [-1,1]: ";
  for (i=0; i<adr.size(); i++) 
    cout << adr[i] << ' ';
  cout << endl << endl; 

#if 0
  double fc[] = {34159104.000053406, -1137459455.9998169, 14585236991.999542, 
  -88999319808., 254089981440.00293, -284648774400.00391, 149276275200.00195, 
  -37308771072.000534, 3605575680.0000544};
#endif


  double fc[] = {-859.0155178152163, 142731.28126273831, -6958791.7981830575, 
  -97808783.000373796, 1514544170.8966789, 202368756.08773908, 
  -3016273948.0315461, -104720971.8528229, 1508707712.1058896};

  cout.precision(16); 

  f.resize(9);
  for (i=0; i<9; i++) f[i] = fc[i]; 

  cout << "Poly: " << f << endl; 
  // rts = roots(f,eps);
  // cout << "Roots: " << rts << endl; 
  rts = roots(f,-1.,1.,eps);
  cout << "Roots in [-1,1]: " << rts << endl; 

  for (i=0; i<rts.size(); i++) rts[i] = evaluate(f, rts[i]); 
  cout << "Values: " << rts << endl; 

#if 0
  chop(f, 1e-7); 

  cout << "Poly: " << f << endl; 
  rts = roots(f,eps);
  cout << "Roots: " << rts << endl; 
#endif 


}
#endif

