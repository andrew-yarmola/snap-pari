/*
 *  positioned_tet.h
 *
 *  The PositionedTet data structure records the position in which a
 *  Tetrahedron is currently being considered.  The file positioned_tet.c
 *  contains routines for working with PositionedTets.
 *
 *  Imagine a tetrahedron sitting on a table with one of its faces
 *  facing toward you.  The face on the table is bottom_face and the one
 *  facing toward you is near_face.  The faces facing away from you
 *  are left_face and right_face.
 *
 *  The orientation field specifies whether the numbering on the
 *  tetrhedron is right_handed or left_handed when viewed in this position.
 */

#ifndef _positioned_tet_
#define _positioned_tet_

#include "SnapPea.h"
#include "kernel_typedefs.h"
#include "triangulation.h"

typedef struct
{
    Tetrahedron *tet;
    FaceIndex   near_face,
		left_face,
		right_face,
		bottom_face;
    Orientation orientation;
} PositionedTet;

#endif
