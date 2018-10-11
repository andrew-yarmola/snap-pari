/* header file for functions which want to make use of the 
   extra structure used by peripheral_curves.c */ 

#ifndef _peripheral_curves_
#define _peripheral_curves_

#include "kernel.h"

typedef struct _PerimeterPiece
{
    Tetrahedron     *tet;
    VertexIndex     vertex;
    FaceIndex       face;
    Orientation     orientation;    /* How the PerimeterPiece sees the tetrahedron  */
    Boolean         checked;
    struct _PerimeterPiece  *mate;  /* the PerimeterPiece this one is glued to . . .    */
    GluingParity    gluing_parity;  /* . . . and how they match up              */
    struct _PerimeterPiece  *next;  /* the neighbor in the counterclockwise direction   */
    struct _PerimeterPiece  *prev;  /* the neighbor in the clockwise        direction   */
} PerimeterPiece;

struct extra
{
    /*
     *  Has this vertex been included in the fundamental domain?
     */
    Boolean                 visited;

    /*
     *  Which vertex of which tetrahedron is its parent in the
     *  tree structure?
     *  (parent_tet == NULL at the root.)
     */
    Tetrahedron             *parent_tet;
    VertexIndex             parent_vertex;

    /*
     *  Which side of this vertex faces the parent vertex?
     *  Which side of the parent vertex faces this vertex? 
     */
    FaceIndex               this_faces_parent,
			    parent_faces_this;

    /*
     *  What is the orientation of this vertex in the
     *  fundamental domain?
     */
    Orientation             orientation;

    /*
     *  Which PerimeterPiece, if any, is associated with
     *  a given edge of the triangle at this vertex?
     *  (As you might expect, its_perimeter_piece[i] refers
     *  to the edge of the triangle contained in face i of
     *  the Tetrahedron.)
     */
    PerimeterPiece          *its_perimeter_piece[4];

    /*
     *  When computing intersection numbers in
     *  adjust_Klein_cusp_orientations() we want to allow for
     *  the possibility that the Triangulation's scratch_curves
     *  are already in use, so we copy them to scratch_curve_backup,
     *  and restore them when we're done.
     */
    int                     scratch_curve_backup[2][2][2][4][4];
};

/*
 *  The following enum lists the six possible gluing
 *  patterns for a torus or Klein bottle.
 */
typedef int GluingPattern;
enum
{
    abAB,   /* square torus                     */
    abcABC, /* hexagonal torus                  */
    abAb,   /* standard square Klein bottle     */
    aabb,   /* P^2 # P^2 square Klein bottle    */
    abcAcb, /* standard hexagonal Klein bottle  */
    aabccB  /* P^2 # P^2 hexagonal Klein bottle */
};

void  attach_extra(Triangulation *manifold);
void  initialize_flags(Triangulation *manifold);
void  pick_base_tet(Triangulation *manifold, Cusp *cusp, Tetrahedron **base_tet, VertexIndex *base_vertex);
void  set_up_perimeter(Tetrahedron *base_tet, VertexIndex base_vertex, PerimeterPiece **perimeter_anchor);
void  expand_perimeter(PerimeterPiece *perimeter_anchor);
void  find_mates(PerimeterPiece *perimeter_anchor);
void  simplify_perimeter(PerimeterPiece **perimeter_anchor);

void  advance_to_next_side(PerimeterPiece **pp);

void  free_perimeter(PerimeterPiece *perimeter_anchor);
void  free_extra(Triangulation *manifold);

#endif
