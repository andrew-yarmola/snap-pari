#ifndef _color_hh_
#define _color_hh_
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


class color {
public:
  double r, g, b, grey; 
  int is_grey;

  color() : grey(0), is_grey(1) {} // black
  color(double greylev) : grey(greylev), is_grey(1) {}
  color(double red, double green, double blue) 
    : r(red), g(green), b(blue), is_grey(0) {}
};

int operator == (const color& a, const color& b);
int operator != (const color& a, const color& b);

color hsbcolor(double hue, double sat = 1.0, double bri = 1.0); 

extern const color black; 

#endif
