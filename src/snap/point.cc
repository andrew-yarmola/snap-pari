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
#include "point.hh"

point operator * (const O31_matrix& m, const point& p)
{
  O31_vector r = m * O31_vector(p); 

  return r; 
}

double hyp_distance(const point& a, const point& b)
{
  O31_vector mink_a = O31_vector(a), mink_b = O31_vector(b);

  mink_a.normalize(); 
  mink_b.normalize();
  double dot = -mink_a * mink_b;
  if (dot < .999) return 1e8; 
  return (dot > 1.0) ? acosh(dot) : 0.0;
}

double euc_distance(const point& a, const point& b)
{
  return (b-a).norm(); 
}

double alt_distance(const point& a, const point& b)
{
  if (a.norm() < .99 && b.norm() < .99)
    return hyp_distance(a,b);
  return 14.1 * euc_distance(a,b); 
}
