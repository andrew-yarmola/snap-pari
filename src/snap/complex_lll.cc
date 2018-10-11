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
#include "snappea/complex.h"
#include <stdio.h>

#include "complex_lll.hh"

using namespace std; 

class complex_lll {
  vector<Complex>& b;   // Elements to be lll reduced. 
  vector<Complex>  bst; // Gram-Schmidt orthogonalization of b. 
  vector<double>   B;   // B[i] == dot(bst[i],bst[i]) 
  vector< vector<double> > mu; 
                        // bst[i] = b[i] - sum(k ; mu[i][k]*bst[k]);
  int k;                // where we are up to so far. 
  int kmax;             // how much of bst is valid. 

  void swap_step(); 
  void red_step(int l); 
public:
  static const double eps; 

  complex_lll(vector<Complex>& _b);
  void reduce(); 
};

void lll_reduce(vector<Complex>& b)
{
  complex_lll LLL(b);
  LLL.reduce();
}

const double complex_lll::eps=1e-8;

inline double dot(const Complex& a, const Complex& b)
{
  return a.real * b.real + a.imag * b.imag;
}

void remove_zeros(vector<Complex>& b)
{
  int n = b.size(); 
  int j = n-1; 
  while (j>=0 && dot(b[j],b[j]) > complex_lll::eps) --j; 
  ++j;
  if (j==0) return; 

  int i = 0; 
  for (; j<n; j++) b[i++] = b[j];
  b.resize(i); 
}

complex_lll::complex_lll(vector<Complex>& _b)
  : b(_b), bst(_b.size()), B(_b.size()), mu(_b.size()), k(1), kmax(0)
{
  int i, n = b.size(); 
  for (i=1; i<n; i++) mu[i].resize(i);
}


void complex_lll::reduce()
{
  int j;
  int n = b.size();

  bst[0] = b[0]; 
  B[0] = dot(b[0],b[0]);

  while (k < n) {

    if (k > kmax) {
      // incremental Gram-Schmidt
      kmax = k;
      
      bst[k] = b[k];
      for (j=0; j<k; j++) {
	if (fabs(B[j]) > eps)
	  mu[k][j] = dot(b[k],bst[j])/B[j];
	else 
	  mu[k][j] = 0.;
	bst[k] -= mu[k][j] * bst[j];
      }
      B[k] = dot(bst[k],bst[k]);

    }

    // test LLL condition
    while (1) {
      red_step(k-1);
      
      if (B[k] >= (0.75 - mu[k][k-1] * mu[k][k-1]) * B[k-1]) break;
      
      swap_step();
      k = max(1,k-1);
    }

    for (j = k-2; j>=0; j--)
      red_step(j);
    k++;
  }
}

void complex_lll::swap_step()
{
  // make bst[k] = projection of b[k] orthogonal to <b[0]...b[k-2]>
  // (just as bst[k-1] is projection of b[k-1] already). 

  bst[k] += mu[k][k-1] * bst[k-1]; 

  swap(b[k-1],  b[k]);
  swap(bst[k-1],bst[k]);

  int j;
  for (j=0; j<k-1; j++)
    swap(mu[k-1][j], mu[k][j]);

  double t, m, BB; 

  m = mu[k][k-1];
  BB = B[k] + m * m * B[k-1]; // the new B[k-1]
  if (BB < eps) { // BB==0 (since all B's are >= 0).
    swap(B[k-1],B[k]);
    for (j=k+1; j <= kmax; j++)
      swap(mu[j][k],mu[j][k-1]);

  } else if (B[k] < eps) {
    B[k-1] = BB;
    mu[k][k-1] = 1.0/m;
    for (j=k+1; j <= kmax; j++)
      mu[j][k-1] /= m;

  } else {
    t          = B[k-1] / BB;
    mu[k][k-1] = t * m;
    B[k]       = t * B[k];
    B[k-1]     = BB;

    for (j=k+1; j <= kmax; j++) {
      t = mu[j][k]; 
      mu[j][k] = mu[j][k-1] - m * t; 
      mu[j][k-1] = t + mu[k][k-1] * mu[j][k];
    }
  }

  bst[k] -= mu[k][k-1] * bst[k-1]; 

}

void complex_lll::red_step(int l)
{
  if (fabs(mu[k][l]) <= 0.5) 
    return ;

  double q = floor(0.5 + mu[k][l]);
  b[k] = b[k] - q * b[l];
  mu[k][l] = mu[k][l] - q;
  int i;
  for (i=0; i<l; i++)
    mu[k][i] = mu[k][i] - q * mu[l][i];
}

bool get_dep(vector<Complex> const& b, Complex z, double C[2])
{
  const float eps = 1e-8; 

  // assume 2-dimensional. 
  double det = b[0].real * b[1].imag - b[0].imag * b[1].real; 
  if (fabs(det) < eps) return false; // expect a proper basis here. 
  C[0] = ( b[1].imag * z.real - b[1].real * z.imag)/det; 
  C[1] = (-b[0].imag * z.real + b[0].real * z.imag)/det; 
  return true; 
}


#if 0
void print_vc(const vector<Complex>& vc)
{
  vector<Complex>::const_iterator vcci;
  cout << "[";
  for (vcci = vc.begin(); vcci != vc.end();) {
    print(*vcci); 
    vcci++;
    if (vcci != vc.end()) cout << ", "; 
  }
  cout << "]"; 
}

int main()
{
  vector<Complex> b(3);

  b[0] = Complex(0.,0.);
  b[1] = Complex(0.,0.);
  b[2] = Complex(0.,0.);

  printf("inital numbers: "); 
  print_vc(b);
  printf("\n");

  lll_reduce(b); 

  printf("final numbers : ");
  print_vc(b);
  printf("\n");

  return 0;
}
#endif
