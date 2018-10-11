#ifndef _R_vector_
#define _R_vector_
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


#include <algorithm>
#include <iostream>
#ifdef __DECCXX
#define exception math_exeption
#include <math.h>
#undef exception
#else 
#include <math.h>
#endif

using std::ostream;
using std::cout;

// double R_vector_epsilon = 1e-9; 
#define R_vector_epsilon 1e-9

template<int n> class R_vector;
template<int n> double operator * (R_vector<n> const& a, R_vector<n> const& b);
template<int n> bool operator == (R_vector<n> const& a, R_vector<n> const& b);
template<int n> bool operator <  (R_vector<n> const& a, R_vector<n> const& b);
template<int n> ostream& operator << (ostream& out, const R_vector<n>& v);

template<int n>
class R_vector {
protected:
  double vec[n]; 

  void copy(const R_vector& v); 
public:
  // static double epsilon;

  R_vector() {}
  R_vector(double x); // All components set to value x. 
  R_vector(int i); // i-th component is 1, others 0.
  R_vector(const double* v);
  R_vector(const R_vector& v) { copy(v); }

  R_vector& operator = (const R_vector& rhs)
    { copy(rhs); return *this; }
  R_vector& operator += (const R_vector& rhs); 
  R_vector& operator -= (const R_vector& rhs); 
  R_vector& operator *= (double r); 
  R_vector& operator /= (double r); 
  R_vector& negate(); 
  R_vector& chop(); // Replace values close to zero with exact zero. 

  bool is_zero() const; 

  // Dot product.
  friend double operator * <n>(R_vector<n> const& a, R_vector<n> const& b); 
  friend bool operator == <n>(R_vector<n> const& a, R_vector<n> const& b);
  friend bool operator <  <n>(R_vector<n> const& a, R_vector<n> const& b);
  friend ostream& operator << <n>(ostream& out, const R_vector<n>& v);

  double norm() const 
    { return sqrt(*this * *this); }
  double norm_squared() const
    { return *this * *this; }
  R_vector& normalize() // Result has norm == 1
    { return *this /= norm(); } 

  double operator [] (int i) const { return vec[i]; }
  double& operator [] (int i) { return vec[i]; }

#ifdef __DECCXX
  // In DEC cxx friends of a non-type-parametrized template must be 
  // defined inline. 
  
  // Euclidean dot
friend double operator * (R_vector<n> const& a, R_vector<n> const& b) 
  {
    int i; 
    double dot = 0.0; 

    for (i=0; i<n; i++)
      dot += a[i] * b[i];
    return dot; 
  }

friend bool operator == (R_vector<n> const& a, R_vector<n> const& b)
  {
    int i; 
    for (i=0; i<n; i++) 
      if (fabs(a.vec[i]-b.vec[i]) > R_vector_epsilon) return false; 
    return true; 
  }

friend bool operator < (R_vector<n> const& a, R_vector<n> const& b)
  {
    int i; 
    for (i=0; i<n; i++) 
      if (fabs(a.vec[i]-b.vec[i]) > R_vector_epsilon) break;
    return i<n && a.vec[i] < b.vec[i]; 
  }

friend ostream& operator << (ostream& out, R_vector<n> const& v)
  {
    int i; 
    out << "[" << v.vec[0];
    for (i=1; i<n; i++) out << ", " << v.vec[i];
    return out << "]";
  }

friend R_vector<n> chop(R_vector<n> m)
{ return m.chop(); }

friend R_vector<n> operator - (const R_vector<n>& a)
{ R_vector<n> neg = a; return neg.negate(); }

friend R_vector<n> operator + (const R_vector<n>& a, const R_vector<n>& b)
{ R_vector<n> sum = a; return sum += b; }

friend R_vector<n> operator - (const R_vector<n>& a, const R_vector<n>& b)
{ R_vector<n> diff = a; return diff -= b; }

friend R_vector<n> operator * (double r, const R_vector<n>& v)
{ R_vector<n> prod = v; return prod *= r; }

friend R_vector<n> operator / (const R_vector<n>& v, double r)
{ R_vector<n> quot = v; return quot /= r; }

#endif 

  void print() const; 
};


template<int n>
R_vector<n>::R_vector(double x)
{
  int i; 
  for (i=0; i<n; i++) vec[i] = x; 
}

template<int n>
R_vector<n>::R_vector(int j)
{
  int i;
  for (i=0; i<n; i++) vec[i] = (i==j) ? 1:0; 
}

template<int n>
R_vector<n>& R_vector<n>::chop()
{
  int i;
  for (i=0; i<n; i++) if (fabs(vec[i]) < R_vector_epsilon) vec[i] = 0.0; 
  return *this;
}

template<int n>
bool R_vector<n>::is_zero() const
{
  int i; 
  for (i=0; i<n; i++) 
    if (fabs(vec[i]) > R_vector_epsilon) return false; 
  return true; 
}

template<int n> 
void R_vector<n>::copy(const R_vector& v)
{
  register int  i;
  for (i = 0; i < n; i++)
    vec[i] = v.vec[i]; 
}

template<int n> 
R_vector<n>::R_vector(const double* v)
{
  register int  i;
  for (i = 0; i < n; i++)
    vec[i] = v[i]; 
}

template<int n> 
R_vector<n>& R_vector<n>::operator += (const R_vector& rhs)
{
  int i;
  for (i=0; i<n; i++)
    vec[i] += rhs.vec[i]; 
  return *this;
}

template<int n> 
R_vector<n>& R_vector<n>::operator -= (const R_vector& rhs)
{
  int i;
  for (i=0; i<n; i++)
    vec[i] -= rhs.vec[i]; 
  return *this;
}

template<int n> 
R_vector<n>& R_vector<n>::operator *= (double r)
{
  int i;
  for (i=0; i<n; i++)
    vec[i] *= r; 
  return *this;
}

template<int n> 
R_vector<n>& R_vector<n>::operator /= (double r)
{
  int i;
  for (i=0; i<n; i++)
    vec[i] /= r; 
  return *this;
}

template<int n> 
R_vector<n>& R_vector<n>::negate()
{
  int i;
  for (i=0; i<n; i++)
    vec[i] = -vec[i]; 
  return *this;
}

template<int n> 
void R_vector<n>::print() const
{
  cout << *this;
}

#ifndef __DECCXX

template<int n> 
double operator * (R_vector<n> const& a, R_vector<n> const& b) // Euclidean dot
{
  int i; 
  double dot = 0.0; 

  for (i=0; i<n; i++)
    dot += a[i] * b[i];
  return dot; 
}

template<int n>
bool operator == (R_vector<n> const& a, R_vector<n> const& b)
{
  int i; 
  for (i=0; i<n; i++) 
    if (fabs(a.vec[i]-b.vec[i]) > R_vector_epsilon) return false; 
  return true; 
}

template<int n> 
bool operator < (R_vector<n> const& a, R_vector<n> const& b)
{
  int i; 
  for (i=0; i<n; i++) 
    if (fabs(a.vec[i]-b.vec[i]) > R_vector_epsilon) break;
  return i<n && a.vec[i] < b.vec[i]; 
}

template<int n> 
ostream& operator << (ostream& out, R_vector<n> const& v)
{
  int i; 
  out << "[" << v.vec[0];
  for (i=1; i<n; i++) out << ", " << v.vec[i];
  return out << "]";
}

template<int n>
inline R_vector<n> chop(R_vector<n> m)
{ return m.chop(); }

template<int n> 
inline R_vector<n> operator - (const R_vector<n>& a)
{ R_vector<n> neg = a; return neg.negate(); }

template<int n> 
inline R_vector<n> operator + (const R_vector<n>& a, const R_vector<n>& b)
{ R_vector<n> sum = a; return sum += b; }

template<int n> 
inline R_vector<n> operator - (const R_vector<n>& a, const R_vector<n>& b)
{ R_vector<n> diff = a; return diff -= b; }

template<int n> 
inline R_vector<n> operator * (double r, const R_vector<n>& v)
{ R_vector<n> prod = v; return prod *= r; }

template<int n> 
inline R_vector<n> operator / (const R_vector<n>& v, double r)
{ R_vector<n> quot = v; return quot /= r; }

#endif // #ifndef __DECCXX

#endif
