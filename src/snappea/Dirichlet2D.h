#ifndef _Dirichlet2D_
#define _Dirichlet2D_

#include "Alg_matrices.h"
#include "polygon.h"

WE2DPolygon *Dirichlet2D_from_generators(AlgMatrix **, int, double, Boolean);
WE2DPolygon *Dirichlet2D(Triangulation *manifold, double vertex_epsilon, Boolean centroid_at_origin, Boolean interactive);

FuncResult Dirichlet2D_extras(WE2DPolygon*);

void print_Dirichlet2D_extras(WE2DPolygon *polygon, FILE* out=stdout, const char* prefix="");

void print_polygon(WE2DPolygon*, FILE* out=stdout);
void print_label(WE2DPolygon*, FILE* out=stdout);
void print_group(WE2DPolygon*, FILE* out=stdout);

void free_Dirichlet2D_domain(WE2DPolygon*);

void set_name(CyclicWord**, int);

#endif

