#ifndef _n_vector_
#define _n_vector_
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


#include <math.h>
#include <iostream>

using std::ostream;

// real vector of arbitrary dimension. 
class n_vector {
public:
  double *V;
  int dim; 

  void copy(n_vector const& v); 
  void copy(const double* v);

  void dehomogenized_copy(n_vector const& v); 

  n_vector() : V(0), dim(0) {} // will behave like zero vector
  n_vector(int d) : dim(d) { V=new double[d]; }
  n_vector(int d, int i);
  n_vector(int d, double diag); 
  n_vector(n_vector const& v) : dim(v.dim) 
  { if (dim) { V=new double[dim]; copy(v); } else V=0; }
  ~n_vector() { delete[] V; }

  n_vector(const char* cs); 

  double& operator [] (int i) { return V[i]; }
  double const& operator [] (int i) const { return V[i]; }

  void operator = (n_vector const& v); 

  void operator *= (double x); 
  void operator /= (double x); 
  void operator += (n_vector const& v);
  void operator -= (n_vector const& v);

  double norm_squared() const;
  double norm() const { return sqrt(norm_squared()); }
  void normalize() { (*this) /= norm(); }

  friend double operator * (n_vector const& a, n_vector const& b); 
  friend n_vector operator + (n_vector const& a, n_vector const& b); 
  friend n_vector operator - (n_vector const& a, n_vector const& b);
  friend n_vector operator * (double x, n_vector const& a);
  friend n_vector operator * (n_vector const& a, double x);
  friend n_vector operator / (n_vector const& a, double x);

  friend void swap(n_vector& a, n_vector& b); 

  friend bool operator == (n_vector const& a, n_vector const& b);

  friend ostream& operator << (ostream& out, n_vector const& a); 
};

#endif
