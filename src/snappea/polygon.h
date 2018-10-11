#ifndef _polygon_
#define _polygon_

#include "Alg_matrices.h"

#define INFINITE_RADIUS     1e20
#define INFINITE_DISTANCE   1e20
#define INFINITE_LENGTH     1e20



typedef int WE2DEdgeEnd;
/*enum
{
    tail = 0,
    tip  = 1
}; */

typedef int WE2DEdgeSide;
/*enum
{
    left  = 0,
    right = 1
};*/


typedef struct WE2DVertex           WE2DVertex;
typedef struct WE2DEdge         WE2DEdge;
/*typedef struct WE2DFace           WE2DFace;*/

typedef struct WE2DVertexClass      WE2DVertexClass;
typedef struct WE2DEdgeClass        WE2DEdgeClass;
/*typedef struct WE2DFaceClass      WE2DFaceClass;*/

typedef struct WE2DPolygon      WE2DPolygon;


struct WE2DVertex
{
    /*
     *  The vector x gives the position of the WE2DVertex in the
     *  projective model of hyperbolic 3-space.  The projective
     *  model is vieWE2Dd as a subset of the Minkowski space model,
     *  so x is a 4-element vector with x[0] == 1.0.
     */
    O31_vector      x;

    O31_vector      xx;

    double          dist;

    Boolean         ideal;

    double          vertex_angle;

    WE2DVertexClass         *v_class;

    Boolean         cone;
	Boolean                 reflector;
    double          distance_to_plane;
    int         which_side_of_plane;

    int         zero_order;
    Boolean         nature;


    WE2DVertex      *prev,
		*next;

};


struct WE2DEdge
{

    WE2DVertex      *v[2];

    WE2DEdge            *e[2];

/*  WE2DFace            *f[2];*/


    double          dihedral_angle;

    double          dist_line_to_origin,
		dist_edge_to_origin;
    O31_vector  closest_point_on_line,
		closest_point_on_edge;

    double          length;

    WE2DEdgeClass       *e_class;

    Boolean         visible;

    WE2DEdge            *neighbor;
    Boolean         preserves_sides,
		preserves_direction,
		preserves_orientation;
	/* some new fields are added here */
	WE2DEdge                 *mate;
	AlgMatrix              *group_element;
      Boolean              to_be_removed;
	  Boolean              clean;
	  Boolean              copied;
	  Boolean              matched;
	  Boolean             assigned;
      char  symbol[2];

    WE2DEdge            *prev,
		*next;

};


/*struct WE2DFace
{
    WE2DEdge            *some_edge;

    
    WE2DFace            *mate;

	O31Matrix               *groupelement;
    double          dist;
    O31_vector      closest_point;

    Boolean         to_be_removed;

    Boolean         clean;

    Boolean         copied;

    Boolean         matched;

    Boolean         visible;

    int         num_sides;

    WE2DFaceClass       *f_class;

    WE2DFace            *prev,
		*next;

};*/




struct WE2DVertexClass
{

	int             index;
    double          hue;

    int         num_elements;

    double          solid_angle;

    int         singularity_order;

    Boolean         ideal;

    double          dist;

    double          min_dist,
		max_dist;

    WE2DVertexClass *belongs_to_region;
    Boolean         is_3_ball;


    WE2DVertexClass *prev,
	    *next;
};


struct WE2DEdgeClass
{

    int         index;
    double          hue;

    int         num_elements;

    double          dihedral_angle;

    int         singularity_order;

    double          dist_line_to_origin,
		dist_edge_to_origin;

    double          length;

    Orbifold2       link;

    double          min_line_dist,
		max_line_dist;

    double          min_length,
		max_length;

    Boolean         removed;
       /* some fields are added here */
	 MatrixParity           parity;
    WE2DEdgeClass       *prev,
		*next;
};


/*struct WE2DFaceClass
{
    
    int             index;
    double          hue;

    
    int             num_elements;

    
    double          dist;

    MatrixParity    parity;

    WE2DFaceClass       *prev,
		    *next;
};*/


struct WE2DPolygon
{
    int             num_vertices,
		    num_edges;
		/*  num_faces;*/

    int             num_finite_vertices,
		    num_ideal_vertices;

    int             num_vertex_classes,
		    num_edge_classes;
		    /*num_face_classes;*/

    int             num_finite_vertex_classes,
		    num_ideal_vertex_classes;

    double          approximate_volume;
     
    double          inradius,
		outradius;


    double          spine_radius;

    double          deviation;

    double          geometric_Euler_characteristic;

    double          vertex_epsilon;

    WE2DVertex      vertex_list_begin,
		    vertex_list_end;
    WE2DEdge            edge_list_begin,
		edge_list_end;
/*  WE2DFace            face_list_begin,
		    face_list_end;*/


    WE2DVertexClass vertex_class_begin,
		vertex_class_end;
    WE2DEdgeClass       edge_class_begin,
		edge_class_end;
/*  WE2DFaceClass       face_class_begin,
		    face_class_end;*/

};

#endif



