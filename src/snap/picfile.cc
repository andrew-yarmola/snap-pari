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
#include "picfile.hh"
#include "bbox.hh"

static const bbox A4(Zero,Complex(596,842));

ps_picfile::ps_picfile(const char *file)
  : current_color(black), scale(-1.)
{
  fp = fopen(file, "w");
  if (!fp) printf("Couldn't open %s for writing\n", file); 
  fprintf(fp, "%%!\n"); 
  set_scale(288.0); 
}

void ps_picfile::set_scale(double sc, Complex tr)
{
  if (!fp) return; 

  if (scale > 0.) // already set?
    fprintf(fp, "grestore\n"); 

  scale = sc; 

  if (tr==Zero) tr = A4.mid(); 

  fprintf(fp, "\ngsave\n%f %f translate\n%f %f scale\n", tr.real, tr.imag, sc, sc); 
  fprintf(fp, "1 %f div setlinewidth\n", sc); 
  fprintf(fp,"/Times-Roman findfont 8 %f div scalefont setfont\n", sc); 
}

void ps_picfile::close()
{
  if (!fp) return; 
 
  fprintf(fp, "grestore\nshowpage\n"); 
  fclose(fp); 
  fp = 0; 
}

void ps_picfile::set_scale(bbox const& user)
{
  double sc;
  Complex tr; 

  A4.get_fit(user, sc, tr); 
  set_scale(sc, tr); 
}

static void make_bbox_path(FILE* fp, bbox const& b)
{
  fprintf(fp, "newpath\n");
  fprintf(fp, "%f %f moveto\n", b.left(), b.bottom());
  fprintf(fp, "%f %f lineto\n", b.right(), b.bottom());
  fprintf(fp, "%f %f lineto\n", b.right(), b.top());
  fprintf(fp, "%f %f lineto\n", b.left(), b.top());
  fprintf(fp, "closepath\n");
}

void ps_picfile::set_scale(bbox const& user, int nr, int nc, int r, int c, bool outline, double aspect, double margin)
{
  double sc;
  Complex tr; 

  bbox out = A4.subbox(nr, nc, r, c, aspect, margin); 
  out.get_fit(user, sc, tr);
  set_scale(sc,tr);

  if (outline) {
    out.rescale(1.0/sc, -tr/sc); 
    make_bbox_path(fp, out); 
    fprintf(fp, "0. setgray stroke\n"); 
  }
}

ps_picfile::~ps_picfile()
{
  close(); 
}

void ps_picfile::set_color(const color& c)
{
  if (c == current_color) return; 
  if (c.is_grey)
    fprintf(fp, "%f setgray\n", c.grey);
  else
    fprintf(fp, "%f %f %f setrgbcolor\n", c.r, c.g, c.b); 
  current_color = c; 
}

void ps_picfile::do_point(const Complex& pt)
{
  fprintf(fp, "%f %f ", pt.real, pt.imag); 
}

void ps_picfile::print_polygon(const list<Complex>& vertices, const color& c, int outline, Complex const& trans)
{
  set_color(c); 

  fprintf(fp, "newpath\n"); 

  list<Complex>::const_iterator it; 

  for (it = vertices.begin(); it != vertices.end(); it++) {

    do_point(*it + trans); 
    if (it == vertices.begin()) 
      fprintf(fp," moveto\n");
    else 
      fprintf(fp," lineto\n");
  }

  if (outline) {
    fprintf(fp,"closepath gsave fill grestore\n");
    set_color(black); 
    fprintf(fp,"stroke\n");
  } else {
    fprintf(fp,"closepath fill\n");
  }
}

void ps_picfile::print_line(const Complex& a, const Complex& b)
{
  set_color(black); 

  fprintf(fp, "newpath\n"); 
  do_point(a); 
  fprintf(fp," moveto\n");
  do_point(b); 
  fprintf(fp," lineto\n");
  fprintf(fp,"stroke\n");
}

void ps_picfile::print_line(const list<Complex>& vertices, Complex const& trans)
{
  set_color(black); 

  fprintf(fp, "newpath\n"); 
  list<Complex>::const_iterator it; 
  for (it = vertices.begin(); it != vertices.end(); it++) {

    do_point(*it + trans); 
    if (it == vertices.begin()) 
      fprintf(fp," moveto\n");
    else 
      fprintf(fp," lineto\n");
  }

  fprintf(fp,"stroke\n");
}

void ps_picfile::print_ellipse(const Complex& ctr, const Complex& radii)
{
  set_color(black);

  fprintf(fp, "gsave newpath\n");
  do_point(ctr); 
  fprintf(fp, "translate\n"); 
  do_point(radii); 
  fprintf(fp, "scale\n"); 
  fprintf(fp, "0 0 1.0 0 360 arc stroke grestore\n"); 
}

void ps_picfile::print_point(const Complex& ctr)
{
  set_color(black);

  double point_size = 3.0;

  do_point(ctr); 
  fprintf(fp, "%f %f div\n", point_size, scale); 
  fprintf(fp, "0 360 arc fill\n"); 
}


void ps_picfile::print_text(const Complex& p, const string& text)
{
  set_color(black); 
  do_point(p); 
  fprintf(fp," moveto (%s) show\n", text.c_str());
}


/* Mathematica pictures */ 

mma_picfile::mma_picfile(const char *file)
: current_color(black),
  point_count(0)
{
  fp = fopen(file, "w");
  if (!fp) {
    printf("failed to open a Mathematica picture file\n"); 
    return; 
  }
  fprintf(fp, "Show[Graphics[{\n"); 
}

void mma_picfile::set_scale(double sc, Complex tr)
{
}

void mma_picfile::close()
{
  if (!fp) return; 
  fprintf(fp, "{}}, AspectRatio->Automatic]]\n"); 
  // The trailing empty list is legal and means we can terminate each 
  // graphic we print (including the last) with a comma. 
  fclose(fp); 
  fp = 0; 
}

mma_picfile::~mma_picfile()
{
  close(); 
}

void mma_picfile::set_color(const color& c)
{
  if (c == current_color) return; 
  if (c.is_grey)
    fprintf(fp, "GrayLevel[%f],\n", c.grey);
  else
    fprintf(fp, "RGBColor[%f,%f,%f],\n", c.r, c.g, c.b); 
  current_color = c; 
}

void mma_picfile::do_point(const Complex& pt)
{
  if (point_count > 2) {
    point_count = 0; 
    fprintf(fp, "\n");
  } 
  fprintf(fp, "{%f, %f}", pt.real, pt.imag); 
  point_count++; 
}

void mma_picfile::print_polygon(const list<Complex>& vertices, const color& c, int outline, Complex const& trans)
{
  set_color(c); 

  list<Complex>::const_iterator it; 

  fprintf(fp, "Polygon[{"); 
  for (it = vertices.begin(); it != vertices.end(); it++) {
    if (it != vertices.begin()) fprintf(fp, ", ");
    do_point(*it + trans); 
  }
  fprintf(fp,"}],\n");

  if (outline) {
    set_color(black); 
    fprintf(fp, "Line[{"); 
    for (it = vertices.begin(); it != vertices.end(); it++) {
      do_point(*it + trans); 
      fprintf(fp, ", ");
    }
    do_point(*vertices.begin()); // close the line
    fprintf(fp,"}],\n");
  }
}

void mma_picfile::print_line(const Complex& a, const Complex& b)
{
  set_color(black); 

  fprintf(fp, "Line[{"); 
  do_point(a); 
  fprintf(fp,",");
  do_point(b); 
  fprintf(fp,"}],\n");
}

void mma_picfile::print_line(const list<Complex>& vertices, Complex const& trans)
{
  set_color(black); 

  list<Complex>::const_iterator it; 

  fprintf(fp, "Line[{"); 
  for (it = vertices.begin(); it != vertices.end(); it++) {
    if (it != vertices.begin()) fprintf(fp, ", ");
    do_point(*it + trans); 
  }
  fprintf(fp,"}],\n");
}

void mma_picfile::print_ellipse(const Complex& ctr, const Complex& radii)
{
  set_color(black);

  fprintf(fp, "Circle[");
  do_point(ctr); 
  fprintf(fp, ","); 
  do_point(radii); 
  fprintf(fp, "],\n"); 
}

void mma_picfile::print_point(const Complex& ctr)
{
  set_color(black);

  fprintf(fp, "Point[");
  do_point(ctr); 
  fprintf(fp, "],\n"); 
}

void mma_picfile::print_text(const Complex& p, const string& text)
{
  set_color(black); 
  fprintf(fp,"Text[\"%s\", ", text.c_str());
  do_point(p);
  fprintf(fp,"],\n"); 
}

