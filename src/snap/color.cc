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
#include "color.hh"

#ifdef __DECCXX
#define exception math_exeption
#include <math.h>
#undef exception
#else
#include <math.h>
#endif

const color black; 

int operator == (const color& a, const color& b)
{
  return (a.is_grey == b.is_grey && 
	  (a.is_grey ? a.grey == b.grey : 
	   (a.r==b.r && a.g==b.g && a.b==b.b))); 
}

int operator != (const color& a, const color& b)
{ return !(a==b); }

static double intensity(double hue)
{
  double in = 2.0 - 6.0 * fabs(hue); 
  if (in > 1.0) return 1.0; 
  if (in < 0.0) return 0.0;
  return in; 
}

color hsbcolor(double hue, double sat, double bri)
{
  return color((sat * (intensity(hue) + intensity(hue - 1.0)) + (1.0 - sat)) * bri, 
	       (sat * intensity(hue - 0.3333333333) + (1.0 - sat)) * bri, 
	       (sat * intensity(hue - 0.6666666666) + (1.0 - sat)) * bri);
}

