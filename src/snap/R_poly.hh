#ifndef _R_poly_
#define _R_poly_
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

#include <vector>
#include <iostream>

using std::vector; 
using std::ostream; 
using std::cout;
using std::endl;
// using std::operator <<;

typedef vector<double> R_poly;

inline int degree(R_poly const& f) 
{ return f.size()-1; }

R_poly operator * (R_poly const& a, R_poly const& b);
R_poly operator + (R_poly const& a, R_poly const& b);
R_poly operator - (R_poly const& a, R_poly const& b);
R_poly operator * (double r, R_poly const& f);
R_poly operator / (R_poly const& f, double r);

R_poly& operator *= (R_poly& f, double r);
R_poly& operator /= (R_poly& f, double r);

R_poly ediv(R_poly const& n, R_poly d, R_poly& r, double eps = 1e-10);
void rmod(R_poly& n, R_poly d, double eps = 1e-10);
R_poly hcf(R_poly f, R_poly g, double eps = 1e-10);

double evaluate(R_poly const& p, double x);

int order(R_poly f, double x, double eps=1e-10);

R_poly sqr(R_poly const& f);
R_poly sqrt(R_poly const& f, double eps = 1e-10);

R_poly derivative(R_poly const& p);
double upper_bound(R_poly const& f);

void chop(R_poly& f, double eps = 1e-10); 
vector<double> roots(R_poly const& f, double eps = 1e-10);
vector<double> roots(R_poly const& f, double xlo, double xhi, double eps = 1e-10);
vector<vector<double> > all_deriv_roots(R_poly f, double eps = 1e-10);
vector<vector<double> > all_deriv_roots(R_poly f, double xlo, double xhi, double eps = 1e-10);
R_poly poly_from_roots(vector<double> const& rts);

ostream& operator << (ostream& out, R_poly const& v);

#endif
