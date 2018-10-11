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

#include "rc_matrix.hh"
#include <vector>
#include <stdio.h>

using namespace std;

double rc_matrix::gauss_eps = 1e-11;

rc_matrix::rc_matrix(int nr, double d)
  : r(nr), c(nr)
{
  M = new double[r*c]; 
  int i, j=0, n=r*c;
  for (i=0; i<n; i++) {
    if (j==0) { M[i]=d; j=r; }
    else { M[i]=0; j--; }
  }
}

void rc_matrix::operator = (rc_matrix const& m)
{
  if (&m == this) return; 

  // ensure M is big enough
  if (r*c < m.r*m.c) {
    delete[] M; 
    M = new double [m.r*m.c]; 
  }

  r = m.r;
  c = m.c; 
  copy(m.M); 
}

void rc_matrix::copy(const double* m)
{
  int i, n = r*c; 
  for (i=0; i<n; i++) M[i] = m[i]; 
}

void rc_matrix::copy(rc_matrix const& m)
{
  int rlim = (m.r < r) ? m.r : r;
  int i, j, n = r*c; 
  if (c==m.c) { // fast easy case
    n = rlim*c; 
    for (i=0; i<n; i++) M[i] = m.M[i]; 
    for (; i<n; i++) M[i] = 0.; 
    return; 
  } 
  int clim = (m.c < c) ? m.c : c; 
  for (i=0; i<rlim; i++) {
    for (j=0; j<clim; j++) (*this)[i][j] = m[i][j];
    for (; j<c; j++) (*this)[i][j] = 0; 
  }
  i *= c; 
  for (; i<n; i++) M[i] = 0.; 
}

n_vector rc_matrix::col(int n) const
{
  n_vector v(r);
  int i;
  for (i=0; i<r; i++) v[i] = (*this)[i][n];
  return v; 
}

void rc_matrix::set_row(int n, n_vector const& v)
{
  int clim = (v.dim < c) ? v.dim : c; 
  int j; 
  for (j=0; j<clim; j++) (*this)[n][j] = v[j];
  for (; j<c; j++) (*this)[n][j] = 0; 
}

double det(rc_matrix m) 
{
  if (m.r != m.c) return 0.; 

  // m is a copy
  int rk = m.gauss(0,0); // make m upper triangular
  if (rk < m.r) return 0.;
  double d = 1.0; 
  int i; 
  for (i=0; i<m.r; i++) d *= m[i][i]; 
  return d; 
}


rc_matrix inverse(rc_matrix m, bool *invertible)
{
  // m is a copy
  if (m.r != m.c) {
    printf("Error, attempt to invert non-square matrix\n"); 
    return rc_matrix();
  }

  rc_matrix inv(m.r, 1.0); // identity matrix. 
  int rk = m.gauss(&inv, 2); 
  if (invertible) *invertible = (rk==m.r); 

  return inv;
}


// third argument default is zero. 

void rc_matrix::swap_rows(int r1, int r2, int sc)
{
  double tmp;
  for (; sc < c; sc++) {
    tmp = (*this)[r1][sc]; 
    (*this)[r1][sc] = (*this)[r2][sc];
    (*this)[r2][sc] = tmp; 
  }
}

// what=0, makes matrix upper triangular
// what=1 gives diagonal matrix
// what=2 gives identity, inverse in m1. 

int rc_matrix::gauss(rc_matrix* m1, int what)
{
  int sr, sc, i, j;

  double maxabs, x;
  int maxrow; 

  if (m1 && m1->r != r) {
    printf("Gauss called for matrices with different numbers of rows\n");
    return -1; 
  }

  sr = 0; sc = 0; 
  while (sr < r && sc < c) {

    // find row with biggest entry in column sc. 
    maxabs = 0.; 
    maxrow = -1; 
    for (i=sr; i<r; i++) {
      x = fabs((*this)[i][sc]); 
      if (x > maxabs + gauss_eps) {
	maxabs = x; 
	maxrow = i; 
      }
    }

    // if all zeros, chop to exact, then go on to next column
    if (maxrow < 0) { 
      for (i=sr; i<r; i++) (*this)[i][sc] = 0.; 
      sc++; 
      continue; 
    }

    // bring it to the top. 
    if (sr!=maxrow) { 
      swap_rows(sr, maxrow, sc); 
      if (m1) m1->swap_rows(sr, maxrow); 
    }

    maxabs = (*this)[sr][sc]; 
    if (what > 1) { // compute full inverse
      for (j=sc+1; j<c; j++) (*this)[sr][j] /= maxabs; 
      if (m1) {
	for (j=0; j<c; j++) {
	  (*m1)[sr][j] /= maxabs; 
	}
      }
      (*this)[sr][sc] = 1.0; 
      maxabs = 1.0; 
    }

    // eliminate column below sr,sc & above if full. 
    i = (what > 0) ? 0 : sr+1; 
    for (; i<r; i++) {
      if (i==sr) continue; 
      x = ((*this)[i][sc])/maxabs; 
      (*this)[i][sc] = 0.; 
      for (j=sc+1; j<c; j++) {
	(*this)[i][j] -= x * (*this)[sr][j]; 
      }
      if (m1) {
	for (j=0; j<c; j++) {
	  (*m1)[i][j] -= x * (*m1)[sr][j]; 
	}
      }
    }

    sr++; sc++; 
  }

  return sr; 
}

// Read in format [a,b;c,d]

rc_matrix::rc_matrix(const char* sc)
  : M(0), r(0), c(0)
{
  const char* cc; 
  while (*sc && *sc=='[' || *sc==' ') sc++; // skip [

  // scan for number of columns. 
  cc = sc; 
  while (*cc && *cc != ']' && *cc!=';') { // immediate close gives c=0; 
    c++; 
    while (*cc && *cc!=',' && *cc!=';' && *cc!=']') cc++; 
    if (*cc==',') cc++; // skip the comma
  }

  // scan for number of rows. 
  cc = sc; 
  while (*cc && *cc != ']') {
    r++; 
    while (*cc && *cc!=';' && *cc!=']') cc++; 
    if (*cc==';') cc++; // skip the semicolon
  }

  if (r*c == 0) return; 

  M = new double[r*c]; 

  // actually read the information. 
  int i=0, j=0; 
  cc = sc; 
  while (*cc && *cc != ']') { // immediate close gives c=0; 
    if (sscanf(cc, "%lf", &M[i*c+j])!=1) break; 
    j++; 
    // find the next punctuation
    while (*cc && *cc!=',' && *cc!=';' && *cc!=']') cc++; 
    if (*cc==',') cc++; // skip the comma
    else if (*cc==';') {
      cc++; // skip the semicolon
      j=0; i++; // start a new row
    }
  }
}


#if 0
int main()
{
  rc_matrix m("[1,2,3;4,5,6]"); 

  cout << m << endl;
  int rk = m.gauss();
  cout << m << " rank " << rk << endl; 

  rc_matrix A("[1,3;2,3]"), B; 

  bool invertible; 
  B = inverse(A, &invertible); 

  cout << A << " invertible=" << invertible << endl;
  cout << B << endl; 

  FILE* fp = fopen("testdata", "r");
  if (!fp) { printf("testdata not found\n"); return 0; }
  char buf[BSIZE];
  vector<n_vector> data; 
  while (fgets(buf, BSIZE, fp)) {
    if (buf[0]=='[') data.push_back(n_vector(buf));
  }
  int i, n = data.size(); 
  printf("read %d vectors\n", n); 
  if (!n) return 0;

  bring_li_to_front(data);
  for (i=0; i<n; i++) cout << data[i] << endl; 

  return 0; 
}
#endif

ostream& operator << (ostream& out, rc_matrix const& a)
{
  int i, j; 
  out << '[';
  for (i=0; i<a.r; i++) {
    if (i>0) out << ';'; 
    for (j=0; j<a.c; j++) {
      if (j>0) { out << ','; }
      out << a[i][j];
    }
  }
  out << ']'; 
  return out; 
}
  
  





