#ifndef _rc_matrix_
#define _rc_matrix_
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


#include "n_vector.hh"

// real matrix of arbitrary dimension. 
struct rc_matrix {
  double *M;
  int r, c; 

  static double gauss_eps; 

  void copy(rc_matrix const& m); 
  void copy(const double* m);

  rc_matrix() : M(0), r(0), c(0) {} // will behave like zero matrix
  rc_matrix(int R, int C) : r(R), c(C) { M=new double[r*c]; }
  rc_matrix(int d, double diag); 
  rc_matrix(rc_matrix const& m) : r(m.r), c(m.c)
    { M=new double[r*c]; copy(m); }
  ~rc_matrix() { delete[] M; }

  rc_matrix(const char* cs); 

  double* operator [] (int row) { return &M[row*c]; }
  const double* operator [] (int row) const { return &M[row*c]; }

  n_vector row(int n) const { n_vector R(c); R.copy(&M[n*c]); return R; }
  n_vector col(int n) const; 

  void set_row(int n, n_vector const& v); 

  void operator = (rc_matrix const& m); 

  void operator *= (double x); 
  void operator /= (double x); 
  void operator += (rc_matrix const& m);
  void operator -= (rc_matrix const& m);

  friend double det(rc_matrix m);
  friend rc_matrix inverse(rc_matrix m, bool* invertible);

  void swap_rows(int r1, int r2, int sc=0);
  int gauss(rc_matrix* m=0, int what=0); 

  friend double operator * (rc_matrix const& a, rc_matrix const& b); 
  friend rc_matrix operator + (rc_matrix const& a, rc_matrix const& b); 
  friend rc_matrix operator - (rc_matrix const& a, rc_matrix const& b);
  friend rc_matrix operator * (double x, rc_matrix const& a);
  friend rc_matrix operator * (rc_matrix const& a, double x);
  friend rc_matrix operator / (rc_matrix const& a, double x);

  friend bool operator == (rc_matrix const& a, rc_matrix const& b);

  friend ostream& operator << (ostream& out, rc_matrix const& a); 
};

#endif
