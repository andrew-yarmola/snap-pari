/*
 *  Dirichlet2D_extras.c
 */
#include "Dirichlet2D_prototypes.h"
#include "fundamental_group.h"

#define DIST_EPSILON    1e-3
#define LENGTH_EPSILON  1e-3
#define IDEAL_EPSILON   4e-7
#define TRACE_ERROR_EPSILON 1e-2
#define PI_EPSILON  1e-1
#define SOLID_ANGLE_EPSILON 1e-4
#define EQUAL    1e-3


static void         edge_classes(WE2DPolygon *polygon);
static void         initialize_edge_classes(WE2DPolygon *polygon);
static void         match_incident_edges(WE2DPolygon *);
static void                     find_edge_mates(WE2DPolygon *);
static void         mirror_edge_classes(WE2DPolygon *polygon, int *count);
static void         make_mirror_edge_class(WE2DPolygon *polygon, WE2DEdge *edge,  int index);
static void         ordinary_edge_classes(WE2DPolygon *polygon, int *count);
static void         make_ordinary_edge_class(WE2DPolygon *polygon, WE2DEdge *edge, int index);
static void         vertex_classes(WE2DPolygon *polygon);
static void         create_vertex_class(WE2DPolygon *polygon, WE2DVertex *vertex);
static void         subdivide_edges_where_necessary(WE2DPolygon *polygon);
/*static void           subdivide_faces_where_necessary(WE2DPolygon *polygon);
static void         cone_face_to_center(WEFace *face, WE2DPolygon *polygon);
static void         bisect_face(WEFace *face, WE2DPolygon *polygon);
static void         delete_face_classes(WE2DPolygon *polygon);
static void         delete_edge_classes(WE2DPolygon *polygon);
static void         delete_vertex_classes(WE2DPolygon *polygon);*/
static void         dihedral_angles(WE2DPolygon *polygon);
/*static void           solid_angles(WE2DPolygon *polygon);*/
static FuncResult           vertex_distances(WE2DPolygon *polygon);
static void         compute_vertex_distance(WE2DVertex *vertex);
static FuncResult           edge_distances(WE2DPolygon *polygon);
static void         compute_edge_distance(WE2DEdge *edge);
/*static void           face_distances(WE2DPolygon *polygon);*/
static FuncResult           edge_lengths(WE2DPolygon *polygon);
static void         compute_edge_length(WE2DEdge *edge);
/*static void           compute_approx_volume(WE2DPolygon *polygon);*/
static void         compute_inradius(WE2DPolygon *polygon);
static void         compute_outradius(WE2DPolygon *polygon);
static void         compute_spine_radius(WE2DPolygon *polygon);
static void         attempt_free_edge_removal(WE2DPolygon *polygon);
static void         compute_deviation(WE2DPolygon *polygon);
static void         compute_geometric_Euler_characteristic(WE2DPolygon *);
Boolean                         o31_vector_equal(O31_vector const&,O31_vector const&);
void                            edge_label(WE2DPolygon*);
void                            find_fixed_points(WE2DPolygon *polygon);
extern void                     make_samedirect(WE2DPolygon *);
extern void simplify_word(CyclicWord  *);
extern void free_cyclic_word(CyclicWord *);

void find_cone_points(WE2DPolygon *);
void find_reflector_points(WE2DPolygon *);

void divide_edge(WE2DEdge *);
extern void count_cells(WE2DPolygon *);
void check_ordinary_angle_sum( WE2DPolygon *polygon);
extern void print_pol_edge(WE2DPolygon *);

void print_orbifold_Euler(WE2DPolygon *, FILE* out, const char* prefix);
void print_cone_points(WE2DPolygon *, FILE* out, const char* prefix);
void print_reflector_points(WE2DPolygon *, FILE* out, const char* prefix);


FuncResult Dirichlet2D_extras(
    WE2DPolygon *polygon)
{
   match_incident_edges(polygon);
   subdivide_edges_where_necessary(polygon);
   edge_classes(polygon);
   vertex_classes(polygon);
  dihedral_angles(polygon);
  if (vertex_distances(polygon) == func_failed)
    {printf("failed in vertex\n");return func_failed;}
  
  if (edge_distances(polygon) == func_failed)
   {printf(" in edge distances"); return func_failed;}
  
  if (edge_lengths(polygon) == func_failed)
    {printf("in edge lengths\n");return func_failed;}
  compute_inradius(polygon);
  compute_outradius(polygon);
  compute_spine_radius(polygon);
  compute_deviation(polygon);
  compute_geometric_Euler_characteristic(polygon);
  edge_label(polygon);
  return func_OK;
}

void print_Dirichlet2D_extras(WE2DPolygon *polygon, FILE* out, const char* prefix)
{
  find_cone_points(polygon);
  find_reflector_points(polygon);
   
  check_ordinary_angle_sum(polygon);

  print_cone_points(polygon, out, prefix);
  print_reflector_points(polygon, out, prefix);

  fprintf(out, "%sEuler characteristic(top) = %lf\n", prefix, polygon->geometric_Euler_characteristic);
  print_orbifold_Euler(polygon, out, prefix);
}

static void edge_classes(
    WE2DPolygon *polygon)
{
   WE2DEdge *edge;
   int      count;
   
   count = 0;
   
  for( edge = polygon->edge_list_begin.next;
      edge != &polygon->edge_list_end;
      edge = edge ->next)
     edge->e_class = NULL;
   mirror_edge_classes(polygon, &count);
   ordinary_edge_classes(polygon, &count);
   polygon->num_edge_classes = count;
}


static void match_incident_edges( WE2DPolygon *polygon)
{
  WE2DEdge  *edge, *mate_edge;
  Boolean     matrix_sign,
	      orientation_preserved,
	      direction_preserved;
  O31_vector tmp;
  for( edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge ->next)
    {
      
      matrix_sign = (determinant((*edge->group_element).matrix) >  ZERO);
      orientation_preserved = matrix_sign;
      direction_preserved = orientation_preserved;
     /* if(edge->mate == edge)
    {
     o31_matrix_times_vector(edge->group_element->matrix,edge->v[tip]->x,tmp) ;
     if(o31_vector_equal(tmp,edge->v[tip]->x))
       printf("points in same direction\n");
     else
       printf("They r iin opposite direction\n");
     printf("matrix_sign = %d \n",matrix_sign);
     print_vector(edge->v[tip]->x);
     printf(" = tip\n");
     print_vector(edge->v[tail]->x);
     printf(" = tail\n");
     print_vector(tmp);
     printf(" = tmp\n");
     exit (0);
       }
      mate_edge = edge->mate;
      edge->neighbor  = mate_edge;
      mate_edge->neighbor = edge;*/
      
      mate_edge = edge->mate;
      edge->preserves_direction = direction_preserved;
      mate_edge->preserves_direction = direction_preserved;
      
      edge->preserves_orientation = orientation_preserved;
      mate_edge->preserves_orientation = orientation_preserved;
    }
}
static void mirror_edge_classes(
		WE2DPolygon *polygon,
		int     *count)
{
    WE2DEdge        *edge;

   for (edge = polygon->edge_list_begin.next;
    edge != &polygon->edge_list_end;
    edge = edge->next)
     if (edge->e_class == NULL  && edge->mate == edge)
     make_mirror_edge_class(polygon, edge, (*count)++);
}

static void make_mirror_edge_class(
    WE2DPolygon *polygon,
    WE2DEdge            *edge,
    int             index)
{
    WE2DEdgeClass   *new_class;
    WE2DEdge        *this_edge,
		*next_edge;



    new_class = NEW_STRUCT(WE2DEdgeClass);
    new_class->index        = index;
    new_class->hue          = index_to_hue(index);
    new_class->num_elements = 0;
    INSERT_BEFORE(new_class, &polygon->edge_class_end);

    this_edge       = edge;


	this_edge->e_class = new_class;
	new_class->num_elements++;
    /*if (edge->preserves_direction == TRUE)
      {
	if (this_edge->preserves_direction[leading_side] == TRUE)
	  new_class->link = orbifold_xnn;
	else
	  new_class->link = orbifold_2xn;
      }
    else
      {
	if (this_edge->preserves_direction[leading_side] == TRUE)
	  new_class->link = orbifold_2xn;
	else
	  new_class->link = orbifold_22n;
      }*/
    
      }


static void ordinary_edge_classes(
    WE2DPolygon *polygon,
    int             *count)
{
    WE2DEdge    *edge;

    
    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

     if (edge->e_class == NULL)

     make_ordinary_edge_class(polygon, edge, (*count)++);
}


static void make_ordinary_edge_class(
    WE2DPolygon *polygon,
    WE2DEdge    *edge,
    int     index)
{
  WE2DEdgeClass *new_class;
  WE2DEdge  *this_edge, *next_edge;
  
  new_class = NEW_STRUCT(WE2DEdgeClass);
  new_class->index      = index;
  new_class->hue            = index_to_hue(index);
  new_class->num_elements   = 0;
  INSERT_BEFORE(new_class, &polygon->edge_class_end);
    
    
    this_edge       = edge;
    
    while (TRUE)
      {
	    if(this_edge  != NULL)
	this_edge->e_class = new_class;
	
	new_class->num_elements++;
	
	next_edge = this_edge->mate;
	
	if (next_edge == edge)
	  {
	
	/*  if (this_edge->preserves_direction[leading_side] == TRUE)
	    new_class->link = orbifold_nn;
	    else
	    new_class->link = orbifold_no;*/
	break;
	  }
	
	this_edge = next_edge;
      }
}


static void vertex_classes(
    WE2DPolygon *polygon)
{
    WE2DVertex  *vertex;

    polygon->num_vertex_classes = 0;


    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)

	vertex->v_class = NULL;

    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)

	if (vertex->v_class == NULL)

	    create_vertex_class(polygon, vertex);
}


static void create_vertex_class(
    WE2DPolygon *polygon,
    WE2DVertex  *vertex)
{
  O31_vector tmp;
  WE2DVertexClass   *new_class;
  Boolean       progress;
  WE2DEdge  *edge;
  WE2DEdgeEnd   which_end;
  WE2DEdge  *mate_edge;
  WE2DEdgeEnd   nbr_end;
   

  new_class         = NEW_STRUCT(WE2DVertexClass);
  new_class->index      = polygon->num_vertex_classes++;
  new_class->hue            = index_to_hue(new_class->index);
  new_class->num_elements   = 0;
  INSERT_BEFORE(new_class, &polygon->vertex_class_end);

  vertex->v_class = new_class;
  new_class->num_elements++;


  do
    {
      progress = FALSE;
      
      for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)
    
    for (which_end = 0; which_end < 2; which_end++) 
	
      if (edge->v[which_end]->v_class == new_class)
	{
	  mate_edge = edge->mate;
	  tmp = mate_edge->group_element->matrix * edge->v[which_end]->x;
	  nbr_end = (o31_vector_equal(tmp,mate_edge->v[which_end]->x))?which_end:!which_end;
	  if (mate_edge->v[nbr_end]->v_class == NULL)
	{
	  mate_edge->v[nbr_end]->v_class = new_class;
	  new_class->num_elements++;
	  progress = TRUE;
	}
	}
      
    } while (progress == TRUE);
}


static void subdivide_edges_where_necessary(
    WE2DPolygon *polygon)
{
    WE2DEdge    *edge;
    O31_vector tmp;
   for (edge = polygon->edge_list_begin.next;
    edge != &polygon->edge_list_end;
    edge = edge->next)
     if((edge == edge->mate))
       {
	 tmp = edge->group_element->matrix * edge->v[tip]->x;
    if(o31_vector_equal(tmp,edge->v[tail]->x))
      divide_edge(edge);
      }
}

/*static void cone_face_to_center(
    WE2DFace            *face,
    WE2DPolygon *polygon)
{
    int         old_num_sides;
    WE2DEdge        **side_edge,
		**radial_edge;
    WE2DFace        **new_face;
    WE2DVertex  *central_vertex;
    O31_vector   fixed_point;
    int         i;

    old_num_sides = face->num_sides;
    if (old_num_sides % 2 != 0)
	uFatalError("cone_face_to_center", "Dirichlet_extras");

    all_edges_counterclockwise(face, TRUE);

    side_edge   = NEW_ARRAY(old_num_sides, WE2DEdge *);
    radial_edge = NEW_ARRAY(old_num_sides, WE2DEdge *);
    new_face    = NEW_ARRAY(old_num_sides, WE2DFace *);

    {
	WE2DEdge    *edge;
	int     count;

	edge = face->some_edge;
	count = 0;
	do
	{
	    side_edge[count++] = edge;
	    edge = edge->e[tip][left];
	} while (edge != face->some_edge);

	if (count != old_num_sides)
	    uFatalError("cone_face_to_center", "Dirichlet_extras");
    }

    for (i = 0; i < old_num_sides; i++)
    {
	radial_edge[i] = NEW_STRUCT(WE2DEdge);
	INSERT_BEFORE(radial_edge[i], &polygon->edge_list_end);
    }


    for (i = 0; i < old_num_sides; i++)
    {
	new_face[i] = NEW_STRUCT(WE2DFace);
	INSERT_BEFORE(new_face[i], face);
    }


    central_vertex = NEW_STRUCT(WE2DVertex);
    INSERT_BEFORE(central_vertex, &polygon->vertex_list_end);

    for (i = 0; i < 4; i++)
	fixed_point[i] = (i == 0 ? ONE : ZERO) + (*face->group_element)[i][0];
    o31_constant_times_vector(ONE/fixed_point[0], fixed_point, central_vertex->x);

    for (i = 0; i < old_num_sides; i++)
    {
	int ip,
	    im,
	    io;

	ip = (i + 1) % old_num_sides;
	im = (i - 1 + old_num_sides) % old_num_sides;
	io = (i + (old_num_sides/2)) % old_num_sides;

	radial_edge[i]->v[tail] = side_edge[i]->v[tip];
	radial_edge[i]->v[tip]  = central_vertex;
	radial_edge[i]->e[tail][left]   = side_edge[i];
	radial_edge[i]->e[tail][right]  = side_edge[ip];
	radial_edge[i]->e[tip ][left]   = radial_edge[im];
	radial_edge[i]->e[tip ][right]  = radial_edge[ip];
	radial_edge[i]->f[left]     = new_face[i];
	radial_edge[i]->f[right]    = new_face[ip];

	side_edge[i]->e[tail][left] = radial_edge[im];
	side_edge[i]->e[tip ][left] = radial_edge[i];
	side_edge[i]->f[left] = new_face[i];

	new_face[i]->some_edge      = side_edge[i];
	new_face[i]->mate           = new_face[io];
	new_face[i]->group_element  = NEW_ARRAY(1,O31_matrix);
	o31_copy(*new_face[i]->group_element, *face->group_element);
	new_face[i]->num_sides      = 3;
    }

    REMOVE_NODE(face);
    my_free(face->group_element);
    my_free(face);
    face = NULL;

    polygon->num_vertices++;
    polygon->num_edges += old_num_sides;
    polygon->num_faces += old_num_sides - 1;
    my_free(side_edge);
    my_free(radial_edge);
    my_free(new_face);
}*/




static void delete_edge_classes(
    WE2DPolygon *polygon)
{
    WE2DEdgeClass   *dead_edge_class;
    WE2DEdge        *edge;

    while (polygon->edge_class_begin.next != &polygon->edge_class_end)
    {
	dead_edge_class = polygon->edge_class_begin.next;
	REMOVE_NODE(dead_edge_class);
	my_free(dead_edge_class);
    }

    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

	edge->e_class = NULL;
}


static void delete_vertex_classes(
    WE2DPolygon *polygon)
{
    WE2DVertexClass *dead_vertex_class;
    WE2DVertex      *vertex;

    while (polygon->vertex_class_begin.next != &polygon->vertex_class_end)
    {
	dead_vertex_class = polygon->vertex_class_begin.next;
	REMOVE_NODE(dead_vertex_class);
	my_free(dead_vertex_class);
    }

    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)

	vertex->v_class = NULL;
}


static void dihedral_angles(
		WE2DPolygon *polygon)
{
  WE2DVertex    *vertex;
  WE2DEdge      *edge,*edge1;
  int           i,
			  j;
  AlgMatrix *m[2];
  O31_vector normal[2];
  double        length,
		 angle_between_normals;
  
  make_samedirect(polygon);
  edge1 = edge = polygon->edge_list_begin.next;
  do{    
    for (i = 0; i < 2; i++)
      {
    if(i==0)
      m[i] = edge->group_element;
    else
      m[i] = edge->e[tip]->group_element;
      
    for (j = 0; j < 4; j++)
      normal[i][j] = m[i]->matrix(j,0);
    normal[i][0] -= ONE;
    length = safe_sqrt(o31_inner_product(normal[i], normal[i]));
    for (j = 0; j < 4; j++)
      normal[i][j] /= length;
      }
    angle_between_normals = safe_acos(o31_inner_product(normal[0], normal[1]));
    edge->v[tip]->vertex_angle =PI - angle_between_normals;
    edge = edge->e[tip];     
     }while (edge !=edge1);
}




static FuncResult vertex_distances(
    WE2DPolygon
    *polygon)
{
    WE2DVertex      *vertex;
    WE2DVertexClass *vertex_class;

    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)

	compute_vertex_distance(vertex);

    for (   vertex_class = polygon->vertex_class_begin.next;
	    vertex_class != &polygon->vertex_class_end;
	    vertex_class = vertex_class->next)
    {
	vertex_class->dist      = ZERO;
	vertex_class->min_dist  = INFINITE_DISTANCE;
	vertex_class->max_dist  = ZERO;
    }


    polygon->num_finite_vertices    = 0;
    polygon->num_ideal_vertices = 0;

    polygon->num_finite_vertex_classes  = 0;
    polygon->num_ideal_vertex_classes   = 0;

    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)
    {
	vertex->v_class->dist += vertex->dist;

	if (vertex->dist < vertex->v_class->min_dist)
	    vertex->v_class->min_dist = vertex->dist;

	if (vertex->dist > vertex->v_class->max_dist)
	    vertex->v_class->max_dist = vertex->dist;

	if (vertex->ideal == FALSE)
	    polygon->num_finite_vertices++;
	else
	    polygon->num_ideal_vertices++;

	vertex->v_class->ideal = vertex->ideal;
    }


    for (   vertex_class = polygon->vertex_class_begin.next;
	    vertex_class != &polygon->vertex_class_end;
	    vertex_class = vertex_class->next)
    {
	vertex_class->dist /= vertex_class->num_elements;

	if (vertex_class->max_dist - vertex_class->min_dist > DIST_EPSILON)
	    return func_failed;

	if (vertex_class->ideal == FALSE)
	    polygon->num_finite_vertex_classes++;
	else
	    polygon->num_ideal_vertex_classes++;
    }


    if (polygon->num_finite_vertex_classes
      + polygon->num_ideal_vertex_classes
     != polygon->num_vertex_classes)

	uFatalError("vertex_distances", "Dirichlet_extras");

    return func_OK;
}


static void compute_vertex_distance(
    WE2DVertex  *vertex)
{

    double      norm_squared;

    norm_squared = o31_inner_product(vertex->x, vertex->x);

    if (norm_squared < - IDEAL_EPSILON)
    {
	vertex->dist    = acosh( safe_sqrt( -ONE / norm_squared ) );
	vertex->ideal   = FALSE;
    }
    else
    {
	vertex->dist    = INFINITE_DISTANCE;
	vertex->ideal   = TRUE;
    }
}


static FuncResult edge_distances(
    WE2DPolygon *polygon)
{
    WE2DEdge        *edge;
    WE2DEdgeClass   *edge_class;

    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

	compute_edge_distance(edge);

    for (   edge_class = polygon->edge_class_begin.next;
	    edge_class != &polygon->edge_class_end;
	    edge_class = edge_class->next)
    {
	edge_class->dist_line_to_origin = ZERO;
	edge_class->dist_edge_to_origin = ZERO;

	edge_class->min_line_dist       = INFINITE_DISTANCE;
	edge_class->max_line_dist       = ZERO;
    }


    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)
    {
	edge->e_class->dist_line_to_origin += edge->dist_line_to_origin;
	edge->e_class->dist_edge_to_origin += edge->dist_edge_to_origin;

	if (edge->dist_line_to_origin < edge->e_class->min_line_dist)
	    edge->e_class->min_line_dist = edge->dist_line_to_origin;

	if (edge->dist_line_to_origin > edge->e_class->max_line_dist)
	    edge->e_class->max_line_dist = edge->dist_line_to_origin;
    }


    for (   edge_class = polygon->edge_class_begin.next;
	    edge_class != &polygon->edge_class_end;
	    edge_class = edge_class->next)
    {
	edge_class->dist_line_to_origin /= edge_class->num_elements;
	edge_class->dist_edge_to_origin /= edge_class->num_elements;

	if (edge_class->max_line_dist - edge_class->min_line_dist > DIST_EPSILON)
	    return func_failed;
    }

    return func_OK;
}


static void compute_edge_distance(
    WE2DEdge    *edge)
{
    O31_vector  p[2],
		v[2],
		w,
		u,
		component;
    double      length,
		projection,
		c[3],
		u_coord,
		p0_coord,
		p1_coord;

    static double   basepoint[4] = {ONE, ZERO, ZERO, ZERO};

    o31_copy_vector(p[0], edge->v[tail]->x);
    o31_copy_vector(p[1], edge->v[tip ]->x);

    o31_vector_sum (p[1], p[0], v[0]);
    o31_vector_diff(p[1], p[0], v[1]);

    length = safe_sqrt( - o31_inner_product(v[0], v[0]) );
    o31_constant_times_vector(ONE/length, v[0], v[0]);

    projection = - o31_inner_product(v[0], v[1]);
    o31_constant_times_vector(projection, v[0], component);
    o31_vector_diff(v[1], component, v[1]);

    length = safe_sqrt(o31_inner_product(v[1], v[1]));
    o31_constant_times_vector(ONE/length, v[1], v[1]);

    o31_copy_vector(w, basepoint);

    c[0] = - o31_inner_product(w, v[0]);
    o31_constant_times_vector(c[0], v[0], component);
    o31_vector_diff(w, component, w);

    c[1] =   o31_inner_product(w, v[1]);
    o31_constant_times_vector(c[1], v[1], component);
    o31_vector_diff(w, component, w);

    c[2] = safe_sqrt(o31_inner_product(w, w));

    o31_vector_diff(basepoint, w, u);
    o31_constant_times_vector(ONE/u[0], u, u);
    o31_copy_vector(edge->closest_point_on_line, u);

    edge->dist_line_to_origin = asinh(c[2]);

    u_coord  = o31_inner_product(v[1], u);
    p0_coord = o31_inner_product(v[1], p[0]);
    p1_coord = o31_inner_product(v[1], p[1]);

    if (p0_coord >= p1_coord)
	uFatalError("compute_edge_distance", "Dirichlet_extras");

    if (u_coord < p0_coord)
    {
	o31_copy_vector(edge->closest_point_on_edge, p[0]);
	edge->dist_edge_to_origin = edge->v[tail]->dist;
    }
    else if (u_coord > p1_coord)
    {
	o31_copy_vector(edge->closest_point_on_edge, p[1]);
	edge->dist_edge_to_origin = edge->v[tip]->dist;
    }
    else
    {
	o31_copy_vector(edge->closest_point_on_edge, edge->closest_point_on_line);
	edge->dist_edge_to_origin = edge->dist_line_to_origin;
    }
}




static FuncResult edge_lengths(
    WE2DPolygon *polygon)
{
    WE2DEdge        *edge;
    WE2DEdgeClass   *edge_class;


    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

	compute_edge_length(edge);


    for (   edge_class = polygon->edge_class_begin.next;
	    edge_class != &polygon->edge_class_end;
	    edge_class = edge_class->next)
    {
	edge_class->length = ZERO;

	edge_class->min_length = INFINITE_LENGTH;
	edge_class->max_length = ZERO;
    }


    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)
    {
	edge->e_class->length += edge->length;

	if (edge->length < edge->e_class->min_length)
	    edge->e_class->min_length = edge->length;

	if (edge->length > edge->e_class->max_length)
	    edge->e_class->max_length = edge->length;
    }



    for (   edge_class = polygon->edge_class_begin.next;
	    edge_class != &polygon->edge_class_end;
	    edge_class = edge_class->next)
    {
	edge_class->length /= edge_class->num_elements;

	if (edge_class->max_length - edge_class->min_length > LENGTH_EPSILON)
	    return func_failed;
    }

    return func_OK;
}


static void compute_edge_length(
    WE2DEdge    *edge)
{
    if (edge->v[tail]->dist == INFINITE_DISTANCE
     || edge->v[tip ]->dist == INFINITE_DISTANCE)

	edge->length = INFINITE_LENGTH;

    else

	edge->length = acosh(
	    -o31_inner_product(edge->v[tail]->x, edge->v[tip]->x)
	    /
	    (
		safe_sqrt(-o31_inner_product(edge->v[tail]->x, edge->v[tail]->x))
	      * safe_sqrt(-o31_inner_product(edge->v[tip ]->x, edge->v[tip ]->x))
	    ));
}



static void compute_inradius(
    WE2DPolygon *polygon)
{

    WE2DEdge    *edge;
    double  min_value;


    min_value = INFINITE_RADIUS;

    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

	if ((edge->group_element)->matrix(0,0) < min_value)

	    min_value = (edge->group_element)->matrix(0,0);

    polygon->inradius = HALF * acosh(min_value);
}


static void compute_outradius(
    WE2DPolygon *polygon)
{

    WE2DVertex  *vertex;
    double      max_projective_distance,
		projective_distance;


    max_projective_distance = ZERO;

    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)
    {
	if (vertex->ideal == TRUE)
	{
	    polygon->outradius = INFINITE_RADIUS;
	    return;
	}


	projective_distance = vertex->x[1] * vertex->x[1]
			    + vertex->x[2] * vertex->x[2]
			    + vertex->x[3] * vertex->x[3];
	if (projective_distance > max_projective_distance)
	    max_projective_distance = projective_distance;
    }


    polygon->outradius = acosh( ONE / safe_sqrt(ONE - max_projective_distance) );
}


static void compute_spine_radius(
    WE2DPolygon *polygon)
{


    WE2DEdgeClass       *edge_class;
    WE2DVertexClass *vertex_class,
		    *vc[2],
		    *region[2];
    double          max_value;
    WE2DEdge            *edge,
		    *max_edge;
    Boolean         union_is_3_ball;

    for (   edge_class = polygon->edge_class_begin.next;
	    edge_class != &polygon->edge_class_end;
	    edge_class = edge_class->next)

	edge_class->removed = FALSE;
  
    for (   vertex_class = polygon->vertex_class_begin.next;
	    vertex_class != &polygon->vertex_class_end;
	    vertex_class = vertex_class->next)
    {
	vertex_class->belongs_to_region = vertex_class;
	vertex_class->is_3_ball
	    = (vertex_class->solid_angle > FOUR*PI - PI_EPSILON);
    }


    while (TRUE)
    {
	max_value = ZERO;

	for (   edge = polygon->edge_list_begin.next;
		edge != &polygon->edge_list_end;
		edge = edge->next)

	    if (edge->e_class->removed == FALSE
	     && edge->e_class->dist_edge_to_origin > max_value)
	    {
		max_edge    = edge;
		max_value   = edge->e_class->dist_edge_to_origin;
	    }

	if (max_value == ZERO)
	    uFatalError("compute_spine_radius", "Dirichlet_extras");

	vc[0] = max_edge->v[0]->v_class;
	vc[1] = max_edge->v[1]->v_class;

	if (vc[0]->belongs_to_region != vc[1]->belongs_to_region
	 && (vc[0]->is_3_ball || vc[1]->is_3_ball))
	{

	    max_edge->e_class->removed = TRUE;

	    region[0] = vc[0]->belongs_to_region;
	    region[1] = vc[1]->belongs_to_region;

	    for (   vertex_class = polygon->vertex_class_begin.next;
		    vertex_class != &polygon->vertex_class_end;
		    vertex_class = vertex_class->next)

		if (vertex_class->belongs_to_region == region[1])

		    vertex_class->belongs_to_region = region[0];

	    union_is_3_ball = (vc[0]->is_3_ball && vc[1]->is_3_ball);

	    for (   vertex_class = polygon->vertex_class_begin.next;
		    vertex_class != &polygon->vertex_class_end;
		    vertex_class = vertex_class->next)

		if (vertex_class->belongs_to_region == region[0])

		    vertex_class->is_3_ball = union_is_3_ball;
	}
	else
	{
	    /*attempt_free_edge_removal(polygon);*/

	    if (max_edge->e_class->removed == TRUE)
	    {
	    }
	    else
	    {
		polygon->spine_radius = max_value;
		return;
	    }
	}
    }

}


/*static void attempt_free_edge_removal(
    WE2DPolygon *polygon)
{
    WE2DFace    *face;
    WE2DEdge    *edge,
	    *remaining_edge;
    int     count;


    for (face = polygon->face_list_begin.next;
	 face != &polygon->face_list_end;
	 face = face->next)
    {

	count = 0;
	remaining_edge = NULL;

	edge = face->some_edge;
	do
	{

	    if (edge->e_class->removed == FALSE)
	    {
		count++;
		remaining_edge = edge;
	    }
	    if (edge->f[left] == face)
		edge = edge->e[tip][left];
	    else
		edge = edge->e[tail][right];

	} while (edge != face->some_edge);
	if (count == 1)
	{

	    if (remaining_edge->e_class->dihedral_angle > TWO*PI - PI_EPSILON)
	    {

		remaining_edge->e_class->removed = TRUE;
		face = &polygon->face_list_begin;
	    }
	}
    }
}*/


static void compute_deviation(
    WE2DPolygon *polygon)
{
    WE2DEdge        *edge;
    double      the_deviation;

    polygon->deviation = ZERO;

    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)
    {
	the_deviation = Alg_deviation(edge->group_element);
	if (the_deviation > polygon->deviation)
	    polygon->deviation = the_deviation;
    }
}


static void compute_geometric_Euler_characteristic(
    WE2DPolygon *polygon)
{
   double c[3];
   int euler;
   int edges = 0;
   WE2DVertex *vertex;
   WE2DEdge   *edge;
   double angle = 0.0;
   
   count_cells(polygon);
   c[0] = ZERO;
   find_cone_points(polygon);
   find_reflector_points(polygon);
   for (vertex = polygon->vertex_list_begin.next;
    vertex != &polygon->vertex_list_end;
    vertex = vertex->next)
     {
    if(vertex->cone == TRUE)
      c[0] += 1.0;
    else 
      if(vertex->reflector == TRUE)
	 c[0] += 1;
      else
	     angle += vertex->vertex_angle;
    }
   c[0] += (angle/(2*PI));

 
   c[1] =polygon->num_edge_classes;

   c[2] = ONE;

   polygon->geometric_Euler_characteristic = c[0] - c[1] + c[2] ;
 }


 void  edge_label(WE2DPolygon *polygon)
{ 
  WE2DEdge  *edge;
  char  c ='a';
  O31_vector trans;
  
  make_samedirect(polygon);
  for(edge = polygon->edge_list_begin.next;
      edge != &polygon->edge_list_end;
      edge = edge->next)
    edge->copied =  FALSE;
  for(edge = polygon->edge_list_begin.next;
      edge  !=&polygon->edge_list_end;
      edge = edge->next)
    {
      if(edge->copied == FALSE)
    {
      edge->symbol[0]=c;
      edge->symbol[1]='\0';
      edge->copied = TRUE;
      if(edge->mate !=edge)
	{
	  trans = edge->group_element->matrix * edge->mate->v[tip]->x;
	 if(o31_vector_equal(edge->v[tip]->x,trans)==TRUE)
	   { edge->mate->symbol[0] = c;
	 edge->mate->symbol[1] ='\0';
	   }
	 else
	   {
	 edge->mate->symbol[0]=(c-'a'+'A');
	 edge->mate->symbol[1]='\0';
	   }
	 edge->mate->copied = TRUE;
       }
      c++;
       }
    }
}

Boolean o31_vector_equal(O31_vector const& a,O31_vector const& b)
{ int i;
  for(i=0;i<4;i++)
    if(fabs(a[i]-b[i])>EQUAL)
      return  FALSE;
  return TRUE;
}


void find_fixed_points(WE2DPolygon *polygon)
{
  WE2DEdge *edge;
  WE2DEdge *edgen;
  WE2DVertex *vertex;
  O31_vector tmp;
  
  for(edge = polygon->edge_list_begin.next;
      edge !=  &polygon->edge_list_end;
      edge = edge->next)
    {
      tmp = edge->group_element->matrix * edge->mate->v[tip]->x;
      vertex = edge->v[tip];
      if(o31_vector_equal(tmp,edge->v[tip]->x))
    { 
      edgen = edge->e[tip];
      tmp = edgen->group_element->matrix * edgen->mate->v[tip]->x; 
      if(o31_vector_equal(tmp,vertex->x))
	vertex->nature = TRUE;
      else
	vertex->nature = FALSE;
    }
      else
    {
      vertex->nature = FALSE;
      //if(edge->mate == edge)
	/*  Fixed point in the midle, need to add that later */
    }
    }
}
     
void print_vector(O31_vector const& x)
{
  printf("{%f,%f}, ",x[1],x[2]);
}

void print_polygon(WE2DPolygon *polygon, FILE* out)
{ 
  WE2DVertex * vertex;

  for(vertex = polygon->vertex_list_begin.next;
      vertex != &polygon->vertex_list_end;
      vertex = vertex->next) {

    fprintf(out, "{%f, %f}", vertex->x[1], vertex->x[2]); 

    if (vertex->next != &polygon->vertex_list_end) fprintf(out, ", "); 
  }
  fprintf(out, "\n");
}


void print_group(WE2DPolygon *polygon, FILE* out)
{
  WE2DEdge *edge,*edge1;
  make_samedirect(polygon);
  edge = edge1=edge = polygon->edge_list_begin.next;

  do {
    simplify_word(edge->group_element->name);
    print_word(edge->group_element->name, out);
    edge = edge->e[tip];
  } while (edge != edge1);
}

void print_label(WE2DPolygon *polygon, FILE* out)
{
  WE2DEdge *edge,*edge1;

  edge = edge1=polygon->edge_list_begin.next;

  do
    {
      fprintf(out, "%c", edge->symbol[0] + 15); // Start from "p" rather than "a". 
      edge = edge->e[tip];
    } while(edge !=edge1);
  fprintf(out, "\n");
}
    

void find_cone_points(WE2DPolygon *polygon)
{
  WE2DEdge *edge,*edge1;
  O31_vector tmp;

  make_samedirect(polygon);
  edge = edge1=polygon->edge_list_begin.next;
  do
    {
      if(edge->mate == edge->e[tip])
    {
      tmp = edge->group_element->matrix * edge->v[tip]->x;
      if(o31_vector_equal(tmp,edge->v[tip]->x) == TRUE)
	edge->v[tip]->cone = TRUE;
       else
	 edge->v[tip]->cone = FALSE;
    }
      else
    edge->v[tip]->cone =  FALSE;
      edge = edge->e[tip];
    } while (edge !=edge1);
}

inline double nint(double x)
{
  return floor(x+.5); 
}
   
void print_cone_points(WE2DPolygon *polygon, FILE* out, const char* prefix)
{
  WE2DVertex *vertex;
  double order;
  int near;
  int done_one = 0; 

  fprintf(out, "%sCone angles:\n%s", prefix, prefix);

  for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next) {

    if (vertex->cone == TRUE) { 

      if (done_one) fprintf(out, ", "); 

      if (vertex->ideal == FALSE) {
	order = (2*PI)/vertex->vertex_angle;
	near = (int)(nint(order));

	if (fabs(order - near) > 0.00001)
	  fprintf(out, "(2*PI / %f)", order);
	else
	  fprintf(out, "(2*PI / %d)", near);
      } else {
	fprintf(out, "(2*PI / infinity)");
      }
      done_one = 1; 
    }
  }

  fprintf(out, "\n");
}

void find_reflector_points(WE2DPolygon *polygon)
{
   WE2DEdge *edge,*edge1;
   
   make_samedirect(polygon);
   edge = edge1 = polygon->edge_list_begin.next;
   do
     {
    if((edge == edge->mate) && (edge->e[tip] == edge->e[tip]->mate))
      edge->v[tip]->reflector = TRUE;
    else
      edge->v[tip]->reflector = FALSE;
    edge = edge->e[tip];
     } while( edge != edge1);
}
   
void print_reflector_points(WE2DPolygon *polygon, FILE* out, const char* prefix)
{
  WE2DVertex *vertex;
  double order;
  int near;
  int done_one = 0; 

  fprintf(out, "%sReflector angles:\n%s", prefix, prefix);

  for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next) {

    if (vertex->reflector == TRUE) { 

      if (done_one) fprintf(out, ", "); 

      if (vertex->ideal == FALSE) {
	order = PI/vertex->vertex_angle;
	near = (int)(nint(order));

	if (fabs(order - near) > 0.00001)
	  fprintf(out, "(PI / %f)", order);
	else
	  fprintf(out, "(PI / %d)", near);
      } else {
	fprintf(out, "(PI / infinity)");
      }
      done_one = 1; 
    }
  }

  fprintf(out, "\n");
}

void check_ordinary_angle_sum(WE2DPolygon *polygon)
{
   WE2DVertex *vertex;
   double angle = 0;
   
   for(vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next) {

     if (vertex->cone == FALSE && vertex->reflector == FALSE)
       angle += vertex->vertex_angle;
   }

   if(fabs(angle - 2* PI) > 0.001)
     if(fabs(angle - 0.0) > 0.001)
       printf("Sum of the ordinary angles is not equal to 2*PI or 0.\n");
}

void divide_edge(WE2DEdge *old_edge)
{
   int i;
   O31_vector cut_point;
   WE2DEdge *new_edge, *tail_edge;
   WE2DVertex *new_vertex;
   
   for(i = 0; i<4 ; i++)
     cut_point[i] = (old_edge->v[tip]->x[i] + old_edge->v[tail]->x[i])/2.0;
   new_edge = NEW_STRUCT(WE2DEdge);
   new_vertex = NEW_STRUCT(WE2DVertex);
   
   new_edge->v[tail] = old_edge->v[tail];
   new_edge->v[tip]  = new_vertex;
   
   new_edge->e[tail] =old_edge->e[tail];
   new_edge->e[tip] = old_edge;
   
   old_edge->v[tail] = new_vertex;
   old_edge->e[tail] = new_edge;
   
   tail_edge = new_edge->e[tail];
   if(tail_edge->v[tip] == new_edge->v[tail])
     tail_edge ->e[tip] = new_edge;
   else 
     if(tail_edge->v[tail] == new_edge->v[tail])
	tail_edge->e[tail] = new_edge;
   else
     uFatalError("Divide_Edge", "Dirichlet2D_extras.c");
   
   new_edge->mate = old_edge;
   old_edge->mate = new_edge;
   
   new_edge->group_element = NEW_ARRAY(1, AlgMatrix);
   Alg_copy(&(new_edge->group_element),old_edge->group_element);
   o31_copy_vector(new_vertex->x,cut_point);
   INSERT_BEFORE(new_edge,old_edge);
   INSERT_BEFORE(new_vertex,new_edge->e[tip]->v[tip]);
}
   
void print_orbifold_Euler(WE2DPolygon *polygon, FILE* out, const char* prefix)
{
  double area, orbi_euler1, orbi_euler2=0.0, angle=0.0, ordinary = 0.0;
  int number;
  WE2DVertex *vertex;
  count_cells(polygon);
  for(vertex = polygon->vertex_list_begin.next;
      vertex != &polygon->vertex_list_end;
      vertex = vertex ->next)
    angle += vertex->vertex_angle;

  number = polygon->num_edges;
  area = (number - 2) * PI - angle ;
  orbi_euler1 = -area/(2*PI);
  
  for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex ->next)
    if (vertex->cone == TRUE || vertex->reflector == TRUE)
      orbi_euler2 += (vertex->vertex_angle)/(2*PI);
    else
      ordinary += vertex->vertex_angle;

  orbi_euler2 += (ordinary / (2*PI)) - ((double)number / 2.0) + 1;

  if (fabs(orbi_euler1 - orbi_euler2)>0.0001)
    printf("Inconsistency in orbifold Euler char computation.\n");
  fprintf(out, "%sOrbifold euler characteristic = %lf\n", prefix, orbi_euler1);
}
