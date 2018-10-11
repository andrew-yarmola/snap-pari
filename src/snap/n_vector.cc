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
#include "n_vector.hh"
#include <cstdio>

using std::sscanf;
using std::cout;
using std::endl;

#define EPSILON 1e-8

// Read a vector in format [1.2,2.0,-1.5]

n_vector::n_vector(const char* cs)
{
  double buf[100];
  int i=0;
  while (*cs && *cs=='[' || *cs==' ') cs++; // skip [
  while (*cs && *cs != ']') {
    if (sscanf(cs, "%lf", &buf[i])!=1) break; 
    i++; 
    while (*cs && *cs!=',' && *cs!=']') cs++; 
    if (*cs) cs++; // skip , or ]
  }
  dim=i;
  V = new double[dim];
  for (i=0; i<dim; i++) V[i]=buf[i]; 
}

n_vector::n_vector(int d, double diag)
  : dim(d)
{ 
  V=new double[d]; 
  int j; 
  for (j=0; j<d; j++) V[j] = diag;
}

void n_vector::operator *= (double x)
{
  int i; 
  for (i=0; i<dim; i++) V[i] *= x;
}

void n_vector::operator /= (double x)
{
  int i; 
  for (i=0; i<dim; i++) V[i] /= x;
}

double n_vector::norm_squared() const
{
  double x=0.;
  int i; 
  for (i=0; i<dim; i++) x += V[i]*V[i]; 
  return x; 
}

ostream& operator << (ostream& out, n_vector const& a)
{
  out << '[';
  int i;
  for (i=0; i<a.dim; i++) {
    if (i>0) out << ',';
    out << a[i]; 
  }
  return out << ']';
}

n_vector operator * (double x, n_vector const& a)
{
  int i; 
  n_vector prod(a.dim); 
  for (i=0; i<a.dim; i++) prod[i] = x*a[i]; 
  return prod;
}

n_vector operator * (n_vector const& a, double x)
{
  int i; 
  n_vector prod(a.dim); 
  for (i=0; i<a.dim; i++) prod[i] = x*a[i]; 
  return prod;
}

n_vector operator / (n_vector const& a, double x)
{
  int i; 
  n_vector res(a.dim); 
  for (i=0; i<a.dim; i++) res[i] = a[i]/x; 
  return res;
}

// different size vectors are OK

n_vector operator + (n_vector const& a, n_vector const& b)
{
  int i, d=a.dim, dm=a.dim; 
  if (b.dim < d) d=b.dim; else dm=b.dim; 
  n_vector sum(dm); 
  for (i=0; i<d; i++) sum[i] = a[i]+b[i]; 
  for (; i<a.dim; i++) sum[i] = a[i]; 
  for (; i<b.dim; i++) sum[i] = b[i]; 
  return sum;
}

n_vector operator - (n_vector const& a, n_vector const& b)
{
  int i, d=a.dim, dm=a.dim; 
  if (b.dim < d) d=b.dim; else dm=b.dim; 
  n_vector diff(dm); 
  for (i=0; i<d; i++) diff[i] = a[i]-b[i]; 
  for (; i<a.dim; i++) diff[i] = a[i]; 
  for (; i<b.dim; i++) diff[i] = -b[i]; 
  return diff;
}

double operator * (n_vector const& a, n_vector const& b)
{
  int i, d=a.dim; 
  double x = 0.; 
  if (b.dim < d) d=b.dim; 
  for (i=0; i<d; i++) x+= a[i]*b[i]; 
  return x;
}

// different length vectors are OK. 
bool operator == (n_vector const& a, n_vector const& b)
{
  int i, d=a.dim; 
  if (b.dim < d) d=b.dim; 
  for (i=0; i<d; i++) 
    if (fabs(a[i]-b[i]) > EPSILON) return false; 
  for (; i<a.dim; i++)
    if (fabs(a[i]) > EPSILON) return false; 
  for (; i<b.dim; i++)
    if (fabs(b[i]) > EPSILON) return false; 
  return true;
}

// lhs dim won't change; argument projected if dim is greater. 

void n_vector::operator += (n_vector const& v)
{
  int i, d=dim; 
  if (v.dim < d) d=v.dim; 
  for (i=0; i<d; i++) V[i] += v[i];
}

void n_vector::operator -= (n_vector const& v)
{
  int i, d=dim; 
  if (v.dim < d) d=v.dim; 
  for (i=0; i<d; i++) V[i] -= v[i];
}

// projects if dim lhs is smaller
// copies projection if dim v is smaller

void n_vector::copy(n_vector const& v)
{
  int i, d=dim; 
  if (v.dim < d) d=v.dim; 
  for (i=0; i<d; i++) V[i]=v[i]; 
  for (; i<dim; i++) V[i]=0;
}

void n_vector::copy(const double* v)
{
  int i;
  for (i=0; i<dim; i++) V[i]=v[i]; 
}

void n_vector::dehomogenized_copy(n_vector const& v)
{
  if (dim != v.dim-1) {
    cout << "dimensions wrong for dehomogenized_copy\n";
    return; 
  }
  double div = v[dim]; 
  if (div < EPSILON) {
    cout << "Warning: attempt to dehomogenize point at infinity\n";
    if (div>=0) div=1e-8; else div=-1e-8; 
  }
  int i; 
  for (i=0; i<dim; i++) V[i]=v[i]/div;
}

// lhs dim may change
void n_vector::operator = (n_vector const& v)
{
  // if (this==&v) return; // self-copy is ok. 
  if (dim!=v.dim) {
    delete[] V;
    dim=v.dim;
    V=new double[dim];
  }
  copy(v.V);
}

// initialize to i'th coordinate vector. 
n_vector::n_vector(int d, int i)
  : dim(d) 
{ 
  V=new double[d]; 
  int j; 
  for (j=0; j<d; j++) V[j] = (i==j) ? 1.0:0.0;
}

void swap(n_vector& a, n_vector& b)
{
  int t;
  double* tmp;

  t = a.dim; a.dim = b.dim; b.dim = t; 
  tmp = a.V; a.V = b.V;     b.V = tmp;
}


#if 0
int main()
{
  n_vector a(2,0), b(2,1);
  n_vector c=a+b;
  cout << "c = " << c << endl; 

  c.normalize(); 
  cout << "c normalized: " << c << endl; 

  cout << "a = " << a << endl; 
  cout << "a.c = " << (a*c) << endl; 
  cout << "a+c = " << (a+c) << endl; 

  n_vector C(3);
  C.copy(c);

  cout << "C = " << C << endl; 

  cout << "C==c is " << (C==c) << endl; 

  cout << "a = " << a << ' ' << endl;

  cout << "a+C = " << (a+C) << endl;

  {
    n_vector z;
    cout << "z = " << z << ' ' << z.V << endl; 
    a.copy(z);
  }

  return 0;
}
#endif
