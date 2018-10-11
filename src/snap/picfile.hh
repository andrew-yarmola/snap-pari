#ifndef _picfile_
#define _picfile_
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

#include <cstdio>
#include <string>
#include <list>
#include "snappea/complex.h"
#include "color.hh"

class bbox; 

using std::string;
using std::list; 

class picfile {
public:
  virtual ~picfile() {} 

  virtual void set_scale(double sc, Complex tr=Zero) {} // does nothing by default
  virtual void set_scale(bbox const& user) {}
  virtual void set_scale(bbox const& user, int nr, int nc, int r, int c, double margin, double aspect, bool outline) {}
  virtual void print_line(const Complex& a, const Complex& b) = 0; 
  virtual void print_line(const list<Complex>& vertices, Complex const& trans=Zero) = 0; 
  virtual void print_ellipse(const Complex& ctr, const Complex& radii) = 0; 
  virtual void print_polygon(const list<Complex>& vertices, const color& c, int outline, Complex const& trans=Zero)=0;
  virtual void print_point(const Complex& p) = 0;
  virtual void print_text(const Complex& p, const string& text) = 0; 
  virtual void close() = 0; 
};

class ps_picfile : public picfile {
  FILE *fp; 
  color current_color; 
  double scale;

  void set_color(const color& c); 
  void do_point(const Complex& pt); 
public:
  ps_picfile(const char *file);
  virtual ~ps_picfile(); 

  virtual void set_scale(double sc, Complex tr=Zero);
  virtual void set_scale(bbox const& user);
  virtual void set_scale(bbox const& user, int nr, int nc, int r, int c, bool outline=false, double aspect=1., double margin=.1);
  virtual void print_line(const Complex& a, const Complex& b); 
  virtual void print_line(const list<Complex>& vertices, Complex const& trans=Zero);
  virtual void print_ellipse(const Complex& ctr, const Complex& radii); 
  virtual void print_polygon(const list<Complex>& vertices, const color& c, int outline, Complex const& trans);
  virtual void print_point(const Complex& p);
  virtual void print_text(const Complex& p, const string& text); 
  virtual void close(); 
};

class mma_picfile : public picfile {
  FILE *fp; 
  color current_color; 
  int point_count; 

  void set_color(const color& c); 
  void do_point(const Complex& pt); 
public:
  mma_picfile(const char *file);
  virtual ~mma_picfile(); 

  virtual void set_scale(double sc, Complex tr=Zero);
  virtual void print_line(const Complex& a, const Complex& b); 
  virtual void print_line(const list<Complex>& vertices, Complex const& trans);
  virtual void print_ellipse(const Complex& ctr, const Complex& radii); 
  virtual void print_polygon(const list<Complex>& vertices, const color& c, int outline, Complex const& trans);
  virtual void print_point(const Complex& p);
  virtual void print_text(const Complex& p, const string& text); 
  virtual void close(); 
};


#endif
