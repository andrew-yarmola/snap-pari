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
#include "O31_line.hh"
#include <algorithm>

double O31_line::epsilon = 1e-4; 

using std::swap; 
using std::cout;
using std::endl;

// Compute a matrix which maps (origin, x-tangent vector, ..) 
// to (pt, dir, ..). 
void O31_line::set_line(O31_vector const& pt, O31_vector const& dir)
{
  // The first column of the matrix will be pt, while the second
  // will be dir. The two remaining columns should be 
  // orthogonal to these wrt. the Minkowski inner product. 
  // Any vector can be made orthogonal to the plane <pt,dir>, 
  // using Gram-Schmidt orthonormalization, as long as it 
  // does not lie in the plane. At least two vectors of the 
  // standard basis for R^4 will not lie in this plane. To get
  // the computation to be numerically as well behaved as possible
  // we choose the two which are already closest to being orthogonal. 
  // But the inner product of the ith standard basis vector with
  // a vector is the ith cpt of the vector. The projection of 
  // the ith standard basis vector into <pt, dir> is therefore
  // pt[i] * pt + dir[i] * dir. Since we are assuming pt, dir to be
  // orthonormal already, this vector has size dir[i]^2 - pt[i]^2. 

  // Find the two standard basis vectors closest to being
  // orthogonal to the plane <pt, dir>. 

  double ip3 = fabs(dir[0] * dir[0] + pt[0] * pt[0]); 
  double ip2 = fabs(dir[1] * dir[1] + pt[1] * pt[1]); 
  double ip; 
  int bvec3 = 0, bvec2 = 1, i;

  for (i=2; i<4; i++) {
    ip = fabs(dir[i] * dir[i] + pt[i] * pt[i]); 
    if (ip < ip3) {
      if (ip3 < ip2) {
	bvec2 = bvec3; ip2 = ip3; 
      }
      ip3 = ip; bvec3 = i;
    } else if (ip < ip2) {
      ip2 = ip; bvec2 = i;
    }
  }

  // Set b2, b3 equal to the two standard basis vectors closest to being
  // orthogonal to the plane <pt, dir>. 

  O31_vector b2, b3; 
  for (i=0; i<4; i++) b2[i] = (bvec2 == i) ? 1.0 : 0.0; 
  for (i=0; i<4; i++) b3[i] = (bvec3 == i) ? 1.0 : 0.0; 

  // Apply Gram-Schmidt orthonormalization to b2 and b3. 

  b2 += (pt * b2) * pt; // Assumes pt*pt == -1. 
  b2 -= (dir * b2) * dir; // Assumes dir*dir == 1. 
  b2 /= sqrt(b2 * b2); // Assumes b2*b2 > 0. Gives b2*b2 == 1. 

  b3 += (pt * b3) * pt; // Assumes pt*pt == -1. 
  b3 -= (dir * b3) * dir; // Assumes dir*dir == 1. 
  b3 -= (b2 * b3) * b2; // Assumes b2*b2 == 1. 
  b3 /= sqrt(b3 * b3); // Assumes b3*b3 > 0. Gives b3*b3 == 1. 

  // Copy the results into the O31_matrix. 

  for (i=0; i<4; i++) trans(i,0) = pt[i];
  for (i=0; i<4; i++) trans(i,1) = dir[i];
  for (i=0; i<4; i++) trans(i,2) = b2[i]; 
  for (i=0; i<4; i++) trans(i,3) = b3[i]; 

  // Make sure determinant is +1. 
  if (determinant(trans) < 0.)
    for (i=0; i<4; i++) swap(trans(i,2),trans(i,3)); 
}

// This constructs a line from two points on the light cone. 
O31_line::O31_line(const O31_vector& a, const O31_vector& b, int)
{
  O31_vector pt(a+b);
  O31_vector dir(b-a); 

  // (a+b), (b-a) are already orthogonal wrt. o31 inner prod. so 
  // all we have to do now is normalize each one individually. 
  pt.normalize(); 
  dir.normalize(); 

  set_line(pt, dir); 
}

O31_line::O31_line(O31_vector p, O31_vector d)
{ 
  p.normalize(); 
  d += (p * d) * p; 
  d.normalize(); 
  set_line(p,d); 
}


O31_line::O31_line(const line& l)
{
  O31_vector a(l.end[0]);
  O31_vector b(l.end[1]);

  O31_vector pt = (a+b);
  O31_vector dir = (b-a);
  pt.normalize(); 
  dir.normalize();

  set_line(pt, dir); 
}

O31_line::operator line () const
{
  O31_vector pt(trans.col(0)), dir(trans.col(1));

  // Ends of the geodesic. 
  O31_vector a(pt - dir), b(pt + dir);

  return line(Complex(a), Complex(b)); 
}

O31_vector O31_line::point(double s, double c) const
{
  O31_vector v;
  int r;
  for (r=0; r<4; r++) 
    v[r] = c * trans(r,0) + s * trans(r,1); 
  return v; 
}

O31_vector O31_line::direction(double s, double c) const
{
  O31_vector v;
  int r;
  for (r=0; r<4; r++) 
    v[r] = s * trans(r,0) + c * trans(r,1); 
  return v; 
}

O31_vector O31_line::end(int i) const
{
  O31_vector v = (i==0) ? 
    (trans.col(0) - trans.col(1)) :
    (trans.col(0) + trans.col(1)); 
  v /= v[0];
  return v;
}

int compare(O31_line const& a, O31_line const& b, double eps)
{
  if (close(a.end(0), b.end(0), eps)) {
    return close(a.end(1), b.end(1), eps) ?  1 : 0;
  } 
  if (close(a.end(0), b.end(1), eps)) {
    return close(a.end(1), b.end(0), eps) ? -1 : 0;
  }
  return 0;
}


void O31_line::adjust_basepoint(O31_vector p)
{
  p.normalize();

  O31_vector d = trans.col(1);
  d += (p*d) * p;
  d.normalize(); 

  int r;
  for (r=0; r<4; r++) {
    trans(r,0) = p[r];
    trans(r,1) = d[r];
  }
}

// We parametrize lines as p(t) = (pt cosh(t) + dir sinh(t)).
// To find the point in the plane F perpendicular to normal, 
// we want to find t such that <p(t), normal> = 0. 
// From this we get tanh(t) = -<pt, normal>/<dir, normal>.
// Putting r = tanh(t), sinh(t) = r/sqrt(1-r^2), and 
// cosh(t) = 1/sqrt(1-r^2). [Clearly these satisfy sinh(t)/cosh(t) = r, 
// and cosh^2(t) - sinh^2(t) = 1.]
//
// Returns 1 if dir and normal point to the same side of F,
// -1 if dir points to the opposite side of F, and 0 if they 
// do not cross (in which case s and c are not set). 

int O31_line::crossing_parameter(O31_vector const& normal, double& s, double& c) const
{
  double dir_dot_normal = direction() * normal;
  if (fabs(dir_dot_normal) < epsilon) {
    return 0; // line and plane meet at infinity. 
  }
  double r = -(point() * normal)/dir_dot_normal; 
  if (r*r >= 1.0) return 0; // line and plane do not meet. 

  double d = sqrt(1.0 - r*r); 

  s = r/d; 
  c = 1.0/d;

  return (dir_dot_normal > 0.0 ? 1 : -1); 
}

bool O31_line::is_in_plane(O31_vector const& normal) const
{
  return (fabs(point() * normal) < epsilon && 
	  fabs(direction() * normal) < epsilon);
}

O31_vector O31_line::projection(O31_vector const& p) const
{ 
  O31_vector pt(trans.col(0)), dir(trans.col(1));
  return -(p * pt) * pt + (p * dir) * dir; 
}


// distance to origin. assumes line is normalized. 
double distance_to_org(const O31_line& l) 
{
  double x = l.trans(0,0) * l.trans(0,0) - l.trans(0,1) * l.trans(0,1);
  if (x <= 1.0) return 0.0;
  return acosh(sqrt(x)); 
}

double distance_to_org(const O31_line& l, const O31_vector& p)
{
  O31_vector q = l.projection(p);
  double x = (p*q)*(p*q)/((p*p)*(q*q));
  if (x <= 1.0) return 0.0;
  return acosh(sqrt(x)); 
}

ostream& operator << (ostream& out, const O31_line& l)
{
  return out << '[' << l.point() << ", " << l.direction() << ']';
}

#ifdef WMAIN

main()
{
  O31_vector pt(1.0, 0.1, 0.3, 0.4), dir(0, -.3, 0, 0.6);
  O31_line l(pt, dir);

  pt = l.point();
  dir = l.direction();

  cout << "Line: " << l << endl;
  cout << "pt * dir: " << (pt*dir) << endl; 

  cout << "Origin: " << O31_origin << endl;

  cout << "distance_to_org(l) = " << distance_to_org(l) << endl; 
  cout << "distance_to_org(l,O31_origin) = " << distance_to_org(l,O31_origin) << endl; 

}

#endif
