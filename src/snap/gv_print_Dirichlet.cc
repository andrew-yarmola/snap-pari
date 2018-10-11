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
#include "gv_print_Dirichlet.hh"
#include "color.hh"
#include "closed_geodesics.hh"

#define DIST_EPSILON 1e-3

// static say xx("gv_print_Dirichlet.cc"); 

static void gv_print_ortholine(FILE* fp, WEPolyhedron* poly,
			O31_matrix const& the_line, Ortholine const& ol, color const& col);
static void gv_print_intervals(FILE* fp, list<interval> const& li, color const& col);

// This is based on fixed_points which does the same thing
// for an SL(2,C) matrix. 

static int real_eigenvectors(O31_matrix const& mx, O31_vector& rep, O31_vector& att)
{
  MoebiusTransformation mt = MoebiusTransformation(mx); 

  line l; 
  int num_fixed = fixed_points(mt.matrix, l);

  if (num_fixed != 2) return num_fixed;

  rep = O31_vector(l.end[0]);
  att = O31_vector(l.end[1]);

  return 2; 
}

static void gv_print(FILE* fp, const color& c) 
{
  if (!c.is_grey) 
    fprintf(fp, "%f %f %f", c.r, c.g, c.b); 
  else 
    fprintf(fp, "%f %f %f", c.grey, c.grey, c.grey);
}

void gv_print(FILE* fp, WEPolyhedron* polyhedron)
{
  WEVertex *vertex;
  WEEdge *edge;
  WEFace *face;

  double x, y, z; 

  fprintf(fp, "{ OFF\n"); 

  fprintf(fp, "%d %d 0\n\n", polyhedron->num_vertices, polyhedron->num_faces);

  int i = 0;

  for (vertex = polyhedron->vertex_list_begin.next;
       vertex != &polyhedron->vertex_list_end;
       vertex = vertex->next) {

    // index the vertices via the zero_order field
    vertex->zero_order = i++;

    // print them too 
    x = vertex->x[1]; 
    y = vertex->x[2]; 
    z = vertex->x[3]; 
    fprintf(fp, "%f %f %f\n", x, y, z); 
  }

  for (face = polyhedron->face_list_begin.next;
       face != &polyhedron->face_list_end;
       face = face->next) {

    // print the number of vertices
    fprintf(fp, "%d ", face->num_sides);


    // print the vertex numbers around the face
    edge = face->some_edge;
    do
    {
	vertex =    (edge->f[left] == face) ?
		    edge->v[tail] :
		    edge->v[tip];

	fprintf(fp, "%d ", vertex->zero_order); 

	edge =  (edge->f[left] == face) ?
		edge->e[tip][left] :
		edge->e[tail][right];

    } while (edge != face->some_edge);

    // color it!
    gv_print(fp, hsbcolor(face->f_class->hue)); 

    fprintf(fp, "\n");

  }

  fprintf(fp, "}\n"); 
}

void gv_print(FILE* fp, WEPolyhedron* polyhedron, O31_matrix const& T)
{
  WEVertex *vertex;
  WEEdge *edge;
  WEFace *face;

  double x, y, z; 

  fprintf(fp, "{ OFF\n"); 

  fprintf(fp, "%d %d 0\n\n", polyhedron->num_vertices, polyhedron->num_faces);

  int i = 0;

  for (vertex = polyhedron->vertex_list_begin.next;
       vertex != &polyhedron->vertex_list_end;
       vertex = vertex->next) {

    // index the vertices via the zero_order field
    vertex->zero_order = i++;

    // print them too 
    x = vertex->x[1]; 
    y = vertex->x[2]; 
    z = vertex->x[3]; 
    fprintf(fp, "%f %f %f\n", x, y, z); 
  }

  for (face = polyhedron->face_list_begin.next;
       face != &polyhedron->face_list_end;
       face = face->next) {

    // print the number of vertices
    fprintf(fp, "%d ", face->num_sides);


    // print the vertex numbers around the face
    edge = face->some_edge;
    do
    {
	vertex =    (edge->f[left] == face) ?
		    edge->v[tail] :
		    edge->v[tip];

	fprintf(fp, "%d ", vertex->zero_order); 

	edge =  (edge->f[left] == face) ?
		edge->e[tip][left] :
		edge->e[tail][right];

    } while (edge != face->some_edge);

    // color it!
    gv_print(fp, hsbcolor(face->f_class->hue)); 

    fprintf(fp, "\n");

  }

  fprintf(fp, "}\n"); 
}

void gv_print_w_axes(FILE* fp, WEPolyhedron* polyhedron)
{
  fprintf(fp, "{ LIST\n"); 

  gv_print(fp, polyhedron); 

  O31_line ln; 
  O31_vector a, b;
  double len; 
  WEFace* face; 
  list<interval> crossings; // Bits where a geodesic meets the polyhedron.

  for (face = polyhedron->face_list_begin.next;
       face != &polyhedron->face_list_end;
       face = face->next) {

    if (real_eigenvectors(*face->group_element, a, b) != 2) continue;

    ln = O31_line(a,b,0); // Get line from it`s endpoints on the light cone. 

    len = complex_length(*face->group_element).real; 

    crossings = list<interval>(); // Empty the list. 

    if (!all_crossing_lifts(polyhedron, ln, 0, crossings, false)) continue; 

    gv_print_intervals(fp, crossings, hsbcolor(face->f_class->hue, 0.7, 1.0)); 
  }

  fprintf(fp, "}\n"); 
}      


void gv_print_w_geodesics(FILE* fp, WEPolyhedron* polyhedron, 
			  vector<GeodesicWord> const& geodesics, 
			  vector<int> const& show)
{
  int num_geodesics = geodesics.size(); 

  fprintf(fp, "{ LIST\n"); 

  gv_print(fp, polyhedron); 

  line l; 
  double len; 
  list<interval> crossings; // Bits where a geodesic meets the polyhedron.

  int i, gn; 
  double cstep = show.size() ? 1.0/double(show.size()) : 0.0; 
  color this_color; 

  for (i=0; i<show.size(); i++) {

    gn = show[i]; 
    if (gn < 0 || gn >= num_geodesics) continue; 

    if (fixed_points(geodesics[gn].mt.matrix, l) != 2) continue; 
    len = geodesics[gn].length.real; 

    crossings.clear();

    if (!all_crossing_lifts(polyhedron, l, 0, crossings, true)) continue; 

    this_color = hsbcolor(cstep * double(i), 0.6, 1.0);
    gv_print_intervals(fp, crossings, this_color); 
  }

  fprintf(fp, "}\n"); 
}

void gv_print_geodesic_w_ortholines(FILE* fp, WEPolyhedron* domain, 
				    GeodesicWord const& gw, 
				    vector<Ortholine> const& ort)
{
  fprintf(fp, "{ LIST\n"); 

  gv_print(fp, domain); 

  line l; 
  if (fixed_points(gw.mt.matrix, l) != 2) {
    fprintf(fp, "}\n"); 
    return; 
  }
  O31_line ln(l); 

  double len = gw.length.real; 

  list<interval> crossings;
  all_crossing_lifts(domain, ln, 0, crossings, false); 

  // Print the geodesic red. 
  gv_print_intervals(fp, crossings, hsbcolor(0.0,0.6,1.0)); 

  color green = hsbcolor(0.333,0.6,1.0); 
  O31_matrix line_mx(ln); 
  int i; 
  for (i=0; i<ort.size(); i++) { 
    gv_print_ortholine(fp, domain, line_mx, ort[i], green); 
  }

  fprintf(fp, "}\n"); 
}

static void gv_print_segment(FILE* fp, O31_vector const& a, O31_vector const& b, const color& col)
{
  if (a[0] < DIST_EPSILON || b[0] < DIST_EPSILON) return; 

  fprintf(fp, "{ VECT 1 2 1 \n2 \n1 \n"); 
  fprintf(fp, "%f %f %f\n", a[1]/a[0], a[2]/a[0], a[3]/a[0]); 
  fprintf(fp, "%f %f %f\n", b[1]/b[0], b[2]/b[0], b[3]/b[0]); 
  gv_print(fp, col); fprintf(fp, " 1.0\n"); // for VECTS it is RGBA (A = opacity)
  fprintf(fp, "}\n"); 
}

static void gv_print_intervals(FILE* fp, list<interval> const& li, color const& col)
{
  list<interval>::const_iterator it;
  for (it = li.begin(); it!=li.end(); it++)
    gv_print_segment(fp, (*it).end(0), (*it).end(1), col); 
}

static const O31Matrix x_to_y_mx = {{1.,0., 0.,0.},
				    {0.,0.,-1.,0.},
				    {0.,1., 0.,0.},
				    {0.,0., 0.,1.}};

static const O31_matrix x_to_y(x_to_y_mx); 

static void gv_print_ortholine(FILE* fp, WEPolyhedron* poly,
			O31_matrix const& the_line, Ortholine const& ol, color const& col)
{
  // First get a segment suitably positioned with respect to the x-axis:
  // Make a segment of length ol.distance() lying along [-1,1], (O31 y-axis).
  // Translate it ol.position() along [0,infinity], (O31 x-axis). 

  // Applying the transformation the_line to it takes the x-axis to 
  // our chosen geodesic, and the segment to an orthogonal segment of
  // that geodesic. 

  // Point on y-axis (=[-1,1]).
  // O31_vector ypt(cosh(ol.distance(), 0., sinh(ol.distance()), 0.));
  // The other point we want is just O31_origin. 

  // O31_matrix mx = the_line * O31_x_trans(ol.position()); 

  // An interval is given as a transformation and a segment 
  // of the x-axis to which it is applied. Currently we have 
  // our ortholine as the image of a segment of the y-axis under mx. 
  // We want a transformation which applies to a segment on the 
  // x-axis instead. Therefore we should multiply on the RIGHT by 
  // a transformation taking the x-axis into the y-axis. 

  O31_matrix mx = the_line * O31_x_trans(ol.position()) * x_to_y; 

  interval orthoseg(mx, 0.0, ol.distance().real, FGWord(), 0);

  // Now chop the interval up into little pieces according
  // to how it looks from the point of view of a fundamental domain. 
  list<interval> crossings;
  all_crossing_lifts(poly, orthoseg, crossings, false); 

  gv_print_intervals(fp, crossings, col); 
}

