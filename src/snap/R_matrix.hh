#ifndef _R_matrix_
#define _R_matrix_
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


#include "R_vector.hh"
#include <cstdio>

using std::cerr;
using std::swap;

template<int n> class R_matrix; 
template<int n> ostream& operator << (ostream& out, R_matrix<n> const& m);

template<int n>
class R_matrix {
protected:
  double mx[n][n];

  void copy(const R_matrix& m); 
public:

  R_matrix() {}
  R_matrix(double x); // Constructs diagonal matrix. 
  R_matrix(R_vector<n> const& v); // Constructs diagonal matrix.
  R_matrix(const R_matrix& m) { copy(m); }
  R_matrix(const double m[][n]);
  R_matrix(int r, int c); // Constructs elementary r,c matrix. 

  double* operator [] (int r) { return mx[r]; }
  const double* operator [] (int r) const { return mx[r]; }
  double operator () (int r, int c) const { return mx[r][c]; }
  double& operator () (int r, int c) { return mx[r][c]; }

  R_matrix& operator = (const R_matrix& rhs)
    { copy(rhs); return *this; }

  R_matrix& operator *= (double r); // multiply by scalar. 
  R_matrix& operator /= (double r); // divide by scalar. 
  R_matrix& operator *= (const R_matrix& m); // multiply on right: *this * m
  R_matrix& left_mul(const R_matrix& m); // multiply on left: m * *this
  R_matrix& operator += (const R_matrix& m); 
  R_matrix& operator -= (const R_matrix& m); 
  
  R_vector<n> row(int r) const; 
  R_vector<n> column(int c) const; 

  void set_row(int r, R_vector<n> const& v);
  void set_column(int c, R_vector<n> const& v);

  double determinant() const; 
  bool is_zero() const; 

  void upper_triangularize(double& sign);
  void gauss(R_matrix& m);
  void gauss(R_vector<n>& m);
  void invert();
  void transpose();
  void negate();
  void chop(); 

  friend ostream& operator << <n>(ostream& out, R_matrix<n> const& m);
  void print(FILE* fp = stdout) const; 
  int read(FILE* fp = stdin); 
  int from_string(const char* str);

#ifdef __DECCXX

friend R_matrix<n> inverse(R_matrix<n> m)
{
  R_matrix<n> res = 1.0; 
  m.gauss(res); // m is just a copy. 
  return res; 
}

friend R_matrix<n> operator * (R_matrix<n> const& a, R_matrix<n> const& b)
{
  R_matrix<n> res;
  int r,i,c;
  for (r=0; r<n; r++) {
    for (c=0; c<n; c++) {
      res[r][c] = 0.0; 
      for (i=0; i<n; i++) res[r][c] += a[r][i] * b[i][c];
    }
  }
  return res; 
}

friend R_vector<n> operator * (R_matrix<n> const& m, R_vector<n> const& v)
{
  R_vector<n> res; 
  int r, i;
  for (r=0; r<n; r++) {
    res[r] = m[r][0] * v[0];
    for (i=1; i<n; i++) res[r] += m[r][i] * v[i]; 
  }
  return res; 
}

friend R_matrix<n> operator * (double r, R_matrix<n> const& m)
{
  R_matrix<n> res(m);
  res *= r;
  return res; 
}

friend R_matrix<n> operator - (R_matrix<n> const& m)
{
  R_matrix<n> res;
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) res[i][j] = -m[i][j];
  }
  return res; 
}

friend R_matrix<n> chop(R_matrix<n> m)
{
  m.chop(); // a copy. 
  return m; 
}

friend R_matrix<n> operator + (R_matrix<n> const& a, R_matrix<n> const& b)
{
  R_matrix<n> res(a);
  return res += b; 
}

friend R_matrix<n> operator - (R_matrix<n> const& a, R_matrix<n> const& b)
{
  R_matrix<n> res(a);
  return res -= b; 
}

friend bool operator == (R_matrix<n> const& a, R_matrix<n> const& b)
{
  int i, j; 
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      if (fabs(a[i][j]-b[i][j]) > R_vector_epsilon) return false; 
  }
  return true; 
}

friend bool operator < (R_matrix<n> const& a, R_matrix<n> const& b)
{
  int i, j;
  bool eq = true; 
  for (i=0; i<n && eq; i++) {
    for (j=0; j<n && eq; j++)
      if (fabs(a[i][j]-b[i][j]) > R_vector_epsilon) eq = false;
  }
  return (!eq) && a[i][j] < b[i][j]; 
}

friend ostream& operator << (ostream& out, R_matrix<n> const& m)
{
  int i, j; 

  out << "[";
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      out << m.mx[i][j]; 
      if (j<n-1) out << ", "; 
      else if (i<n-1) out << "; "; 
    }
  }
  return out << "]"; 
}
#endif

};

template<int n>
void R_matrix<n>::copy(const R_matrix<n>& m)
{
  int i,j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      mx[i][j] = m.mx[i][j];
}

template<int n>
R_matrix<n>::R_matrix(double x) // Constructs diagonal matrix.
{
  int i,j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      mx[i][j] = (i == j) ? x : 0.0;
}

template<int n>
R_matrix<n>::R_matrix(R_vector<n> const& v) // Constructs diagonal matrix.
{
  int i,j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      mx[i][j] = (i == j) ? v[i] : 0.0;
}

template<int n>
R_matrix<n>::R_matrix(int r, int c) // Constructs matrix with 1 in row r, col c. 
{
  int i,j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      mx[i][j] = 0.0;
  mx[r][c] = 1.0; 
}

template<int n>
R_matrix<n>::R_matrix(const double m[][n])
{
  int i,j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      mx[i][j] = m[i][j];
}

template<int n>
void R_matrix<n>::upper_triangularize(double& sign)
{
  int r=0, c, i, rr, rmax; 
  double absmax, mult;
  sign = 1.0; 

  for (c=0; c<n; c++) {
    // find biggest entry in column c, below row r; let that be row rmax.
    absmax = 0.;
    for (i=r; i<n; i++) {
      if (fabs(mx[i][c]) > absmax) {
	absmax = fabs(mx[i][c]); rmax = i; }
    }
    // if all are zero proceed to next column. 
    if (absmax < R_vector_epsilon) continue; 

    // bring rmax into row r. 
    if (rmax != r) {
      sign = -sign;
      for (i=c; i<n; i++) swap(mx[r][i], mx[rmax][i]);
    }

    // subtract appropriate multiple of row r from each row below, 
    // putting exact zero into entries below row r col c. 
    for (rr=r+1; rr<n; rr++) {
      mult = mx[rr][c]/mx[r][c];
      mx[rr][c] = 0.0; 
      for (i=c+1; i<n; i++) mx[rr][i] -= mult * mx[r][i];
    }
    r++; 
  }
}

template<int n>
void R_matrix<n>::gauss(R_matrix<n>& m)
{
  int r=0, c, i, rr, rmax; 
  double absmax, mult, div;

  for (c=0; c<n; c++) {
    // find biggest entry in column c, below row r; let that be row rmax.
    absmax = 0.;
    for (i=r; i<n; i++) {
      if (fabs(mx[i][c]) > absmax) {
	absmax = fabs(mx[i][c]); rmax = i; }
    }

    // if all are zero give up. 
    if (absmax < R_vector_epsilon) {
      cerr << "Warning: non-invertible matrix in gauss\n";
      return; 
    }

    // bring rmax into row r. 
    if (rmax != r) {
      for (i=c; i<n; i++) swap(mx[r][i], mx[rmax][i]); 
      for (i=0; i<n; i++) swap(m[r][i], m[rmax][i]); // do same to m. 
    }

    // divide row r by its first non-zero entry. 
    div = mx[r][c]; 
    mx[r][c] = 1.0; 
    for (i=c+1; i<n; i++) mx[r][i] /= div; 
    for (i=0; i<n; i++) m[r][i] /= div; 

    // subtract appropriate multiple of row r from each row above and below, 
    // putting exact zero into entries below row r col c. 
    for (rr=0; rr<n; rr++) {
      if (rr == r) continue;
      mult = mx[rr][c];
      mx[rr][c] = 0.0; 
      for (i=c+1; i<n; i++) mx[rr][i] -= mult * mx[r][i];
      for (i=0; i<n; i++) m[rr][i] -= mult * m[r][i];
    }
    r++; 
  }
}

template<int n>
void R_matrix<n>::gauss(R_vector<n>& v)
{
  int r=0, c, i, rr, rmax; 
  double absmax, mult, div;

  for (c=0; c<n; c++) {
    // find biggest entry in column c, below row r; let that be row rmax.
    absmax = 0.;
    for (i=r; i<n; i++) {
      if (fabs(mx[i][c]) > absmax) {
	absmax = fabs(mx[i][c]); rmax = i; }
    }

    // if all are zero give up. 
    if (absmax < R_vector_epsilon) {
      cerr << "Warning: non-invertible matrix in gauss\n";
      return;
    }

    // bring rmax into row r. 
    if (rmax != r) {
      for (i=c; i<n; i++) swap(mx[r][i], mx[rmax][i]); 
      swap(v[r], v[rmax]); 
    }

    // divide row r by its first non-zero entry. 
    div = mx[r][c]; 
    mx[r][c] = 1.0; 
    for (i=c+1; i<n; i++) mx[r][i] /= div; 
    v[r] /= div; 

    // subtract appropriate multiple of row r from each row above and below, 
    // putting exact zero into entries below row r col c. 
    for (rr=0; rr<n; rr++) {
      if (rr == r) continue;
      mult = mx[rr][c];
      mx[rr][c] = 0.0; 
      for (i=c+1; i<n; i++) mx[rr][i] -= mult * mx[r][i];
      v[rr] -= mult * v[r];
    }
    r++; 
  }
}

template<int n>
double R_matrix<n>::determinant() const
{
  R_matrix<n> m(*this); // Copy ourself. 
  double sign, det; 
  m.upper_triangularize(sign); 
  det = m[0][0] * sign; 
  int i; 
  for (i=1; i<n; i++) det *= m[i][i]; 
  return det; 
}

template<int n>
void R_matrix<n>::invert()
{
  R_matrix<n> m = 1.0; // Identity matrix. 
  gauss(m); 
  copy(m);
}

template<int n>
void R_matrix<n>::transpose()
{
  int i, j; 
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) swap(mx[i][j], mx[j][i]); 
}

template<int n>
R_matrix<n>& R_matrix<n>::operator *= (R_matrix<n> const& b)
{
  R_matrix<n> tmp = (*this) * b;
  copy(tmp);
  return *this;
}

template<int n>
R_matrix<n>& R_matrix<n>::left_mul(R_matrix<n> const& a)
{
  R_matrix<n> tmp = a * (*this);
  copy(tmp);
  return *this;
}

template<int n>
R_matrix<n>& R_matrix<n>::operator *= (double r)
{
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) mx[i][j] *= r; 
  }
  return *this;
}

template<int n>
void R_matrix<n>::negate()
{
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) mx[i][j] = -mx[i][j]; 
  }
}

template<int n>
void R_matrix<n>::chop()
{
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) 
      if (fabs(mx[i][j]) < R_vector_epsilon) mx[i][j] = 0.0; 
  }
}

template<int n>
R_matrix<n>& R_matrix<n>::operator /= (double r)
{
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) mx[i][j] /= r; 
  }
  return *this;
}

template<int n>
R_matrix<n>& R_matrix<n>::operator += (R_matrix<n> const& m)
{
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) mx[i][j] += m.mx[i][j]; 
  }
  return *this;
}

template<int n>
R_matrix<n>& R_matrix<n>::operator -= (R_matrix<n> const& m)
{
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) mx[i][j] -= m.mx[i][j]; 
  }
  return *this;
}

template<int n>
R_vector<n> R_matrix<n>::row(int r) const
{
  R_vector<n> res; 
  int i;
  for (i=0; i<n; i++) res[i] = mx[r][i];
  return res;
}

template<int n>
R_vector<n> R_matrix<n>::column(int c) const
{
  R_vector<n> res; 
  int i;
  for (i=0; i<n; i++) res[i] = mx[i][c];
  return res;
}

template<int n>
void R_matrix<n>::set_row(int r, R_vector<n> const& v)
{
  int i;
  for (i=0; i<n; i++) mx[r][i] = v[i];
}

template<int n>
void R_matrix<n>::set_column(int c, R_vector<n> const& v)
{
  int i;
  for (i=0; i<n; i++) mx[i][c] = v[i];
}

template<int n>
bool R_matrix<n>::is_zero() const
{
  int i, j; 
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      if (fabs(mx[i][j]) > R_vector_epsilon) return false; 
  }
  return true; 
}

template<int n> 
void R_matrix<n>::print(FILE* fp) const
{
  int i, j; 

  fprintf(fp, "[");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      fprintf(fp, "%lf", mx[i][j]); 
      if (j<n-1) fprintf(fp, ", "); 
      else if (i<n-1) fprintf(fp, "; "); 
    }
  }
  fprintf(fp, "]"); 
}

template<int n>
int R_matrix<n>::read(FILE* fp)
{
  char c[2]; 
  int i, j; 

  fscanf(fp, " %1s ", c); if (c[0] != '[') return 0; 
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      if (j<n-1) {
	if (fscanf(fp, " %lf , ", &mx[i][j]) != 1) return 0; 
      } else {
	fscanf(fp, " %lf %1s \n", &mx[i][j], c); 
	if ((i<n-1 && c[0] != ';') || (i==n-1 && c[0] != ']')) return 0; 
      }
    }
  }
  return 1; 
}

template<int n>
int R_matrix<n>::from_string(const char* str)
{
  char delim;
  const char* s = str;
  int i, j; 

  while (*s && *s!='[') s++;
  if (!*s) return 0; 
  s++;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      if (sscanf(s, " %lf ", &mx[i][j]) != 1) return 0; 
      delim = (j<n-1) ? ',' : ((i<n-1) ? ';' : ']');
      while (*s && *s!=delim) s++;
      if (!*s) return 0; 
      s++;
    }
  }
  return s - str; 
}

#ifndef __DECCXX

template<int n>
R_matrix<n> inverse(R_matrix<n> m)
{
  R_matrix<n> res = 1.0; 
  m.gauss(res); // m is just a copy. 
  return res; 
}

template<int n>
R_matrix<n> operator * (R_matrix<n> const& a, R_matrix<n> const& b)
{
  R_matrix<n> res;
  int r,i,c;
  for (r=0; r<n; r++) {
    for (c=0; c<n; c++) {
      res[r][c] = 0.0; 
      for (i=0; i<n; i++) res[r][c] += a[r][i] * b[i][c];
    }
  }
  return res; 
}

template<int n>
R_vector<n> operator * (R_matrix<n> const& m, R_vector<n> const& v)
{
  R_vector<n> res; 
  int r, i;
  for (r=0; r<n; r++) {
    res[r] = m[r][0] * v[0];
    for (i=1; i<n; i++) res[r] += m[r][i] * v[i]; 
  }
  return res; 
}

template<int n>
R_matrix<n> operator * (double r, R_matrix<n> const& m)
{
  R_matrix<n> res(m);
  res *= r;
  return res; 
}

template<int n>
R_matrix<n> operator - (R_matrix<n> const& m)
{
  R_matrix<n> res;
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) res[i][j] = -m[i][j];
  }
  return res; 
}

template<int n>
R_matrix<n> chop(R_matrix<n> m)
{
  m.chop(); // a copy. 
  return m; 
}

template<int n>
R_matrix<n> operator + (R_matrix<n> const& a, R_matrix<n> const& b)
{
  R_matrix<n> res(a);
  return res += b; 
}

template<int n>
R_matrix<n> operator - (R_matrix<n> const& a, R_matrix<n> const& b)
{
  R_matrix<n> res(a);
  return res -= b; 
}

template<int n>
bool operator == (R_matrix<n> const& a, R_matrix<n> const& b)
{
  int i, j; 
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      if (fabs(a[i][j]-b[i][j]) > R_vector_epsilon) return false; 
  }
  return true; 
}

template<int n> 
bool operator < (R_matrix<n> const& a, R_matrix<n> const& b)
{
  int i, j;
  bool eq = true; 
  for (i=0; i<n && eq; i++) {
    for (j=0; j<n && eq; j++)
      if (fabs(a[i][j]-b[i][j]) > R_vector_epsilon) eq = false;
  }
  return (!eq) && a[i][j] < b[i][j]; 
}

template<int n> 
ostream& operator << (ostream& out, R_matrix<n> const& m)
{
  int i, j; 

  out << "[";
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      out << m.mx[i][j]; 
      if (j<n-1) out << ", "; 
      else if (i<n-1) out << "; "; 
    }
  }
  return out << "]"; 
}

#endif // __DECCXX

#endif
