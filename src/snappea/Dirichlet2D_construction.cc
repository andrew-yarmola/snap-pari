
/*
 *  Dirichlet2D_construction.c 
 */

#include "Dirichlet2D_prototypes.h"
#include <stdlib.h>     /* needed for qsort() */

#define SQUARE_SIZE     2.0

#define ERROR_EPSILON       1e-4
#define HYPERIDEAL_EPSILON  1e-4

#define VERIFY_EPSILON      1e-4 

#define DEVIATION_EPSILON   1e-3


static WE2DPolygon  *initial_polygon(AlgMatrixPairList *gen_list, double vertex_epsilon);
static WE2DPolygon  *new_WE2DPolygon(void);
static void     make_square(WE2DPolygon *polygon);
static FuncResult   slice_polygon(WE2DPolygon *polygon, AlgMatrixPairList *gen_list);
static FuncResult   intersect_with_halfplanes(WE2DPolygon *polygon, AlgMatrixPair *matrix_pair);
static Boolean      same_image_of_origin(O31_matrix const& m0, O31_matrix const& m1);
static FuncResult   slice_with_line(WE2DPolygon *polygon, AlgMatrix m, WE2DEdge **new_edge);
static void split_edge(WE2DEdge *old_edge, O31_vector const& cut_point, Boolean set_Dirichlet_construction_fields);
static FuncResult   compute_normal_to_Dirichlet_plane(O31_matrix const& m, O31_vector& normal_vector);
static void     compute_vertex_to_line_distances(WE2DPolygon *polygon, O31_vector const& normal_vector);
static Boolean      positive_vertices_exist(WE2DPolygon *polygon);
static void     cut_edges(WE2DPolygon *polygon);
static FuncResult   create_edge(WE2DPolygon *polygon,WE2DEdge **);
static FuncResult   check_topology_of_cut(WE2DPolygon *polygon);
static void     install_new_edge(WE2DPolygon *polygon, WE2DEdge *new_edge);
static Boolean      edge_still_exists(WE2DPolygon *polygon, WE2DEdge *edge);
static Boolean      has_hyperideal_vertices(WE2DPolygon *polygon);
static void     compute_all_products(WE2DPolygon *polygon, AlgMatrixPairList *product_list);
static void     poly_to_current_list(WE2DPolygon *polygon, AlgMatrixPairList *current_list);
static void     current_list_to_product_tree(AlgMatrixPairList *current_list, AlgMatrixPair  **product_tree);
static Boolean      already_on_product_tree(AlgMatrix *product, AlgMatrixPair *product_tree);
static void     add_to_product_tree(AlgMatrix *product, AlgMatrixPair **product_tree);
static void     product_tree_to_product_list(AlgMatrixPair  *product_tree, AlgMatrixPairList *product_list);
static void     append_tree_to_list(AlgMatrixPair *product_tree, AlgMatrixPair *list_end);
static FuncResult   check_edges(WE2DPolygon *polygon);
static FuncResult   pare_edge(WE2DEdge *edge, WE2DPolygon *polygon, Boolean *edge_was_pared);
static FuncResult   pare_mated_edge(WE2DEdge *edge, WE2DPolygon *polygon, Boolean *edge_was_pared);
static FuncResult   pare_mateless_edge(WE2DEdge *edge, WE2DPolygon *polygon, Boolean *edge_was_pared);
static FuncResult   try_this_alpha(AlgMatrix *alpha, WE2DEdge *edge, WE2DPolygon *polygon, Boolean *edge_was_pared);
void        count_cells(WE2DPolygon *polygon);
static void     sort_edges(WE2DPolygon *polygon);
static int CDECL    compare_edge_distance(const void *ptr1, const void *ptr2);
static FuncResult   verify_group(WE2DPolygon *polygon, AlgMatrixPairList *gen_list);
static void     rewrite_gen_list(WE2DPolygon *polygon, AlgMatrixPairList *gen_list);
void    make_samedirect(WE2DPolygon *polygon);
void redirect2D_edge(WE2DEdge *,Boolean);
void print_pol_edge(WE2DPolygon *);
extern void print_vector(O31_vector const&  );
void print_edge(WE2DEdge *);
extern void simplify_word(CyclicWord *);

WE2DPolygon *compute_Dirichlet2D_domain(
    AlgMatrixPairList   *gen_list,
    double      vertex_epsilon)
{
    WE2DPolygon *polygon;
    WE2DEdge * edge;
    polygon = initial_polygon(gen_list, vertex_epsilon);
    if (polygon == NULL)
	return NULL;
	if (check_edges(polygon) == func_failed)
    {  
	free_Dirichlet2D_domain(polygon);
	return NULL;
    }
    count_cells(polygon);
    sort_edges(polygon);

    if (verify_group(polygon, gen_list) == func_failed)
    {
	free_Dirichlet2D_domain(polygon);
	return NULL;
    }


    rewrite_gen_list(polygon, gen_list);
	return polygon;
}


static WE2DPolygon  *initial_polygon(
    AlgMatrixPairList   *gen_list,
    double      vertex_epsilon)
{
    WE2DPolygon *polygon;
    AlgMatrixPairList   product_list;
    int cnt=0;
    polygon = new_WE2DPolygon();

    polygon->vertex_epsilon = vertex_epsilon;

    make_square(polygon);

	if (slice_polygon(polygon, gen_list) == func_failed)
    {
     
	free_Dirichlet2D_domain(polygon);
	return NULL;
    }

    while (has_hyperideal_vertices(polygon) == TRUE)
    { 
      cnt ++;
      compute_all_products(polygon, &product_list);
      if (slice_polygon(polygon, &product_list) == func_failed)
	{
	  free_matrix_pairs(&product_list);
	  free_Dirichlet2D_domain(polygon);
	  return NULL;
	}
      free_matrix_pairs(&product_list);
      if(cnt > 11)
	{
	  printf("\n** not enough generators for the following orbifold**\n");
/*        print_pol_edge(polygon);*/
	  if(polygon != NULL)
	  free_Dirichlet2D_domain(polygon);
	  return NULL;
	}
    }

    return polygon;
}


static WE2DPolygon *new_WE2DPolygon()
{
    WE2DPolygon *new_polygon;

    new_polygon = NEW_STRUCT(WE2DPolygon);

    new_polygon->num_vertices   = 0;
    new_polygon->num_edges      = 0;

    new_polygon->vertex_list_begin.prev = NULL;
    new_polygon->vertex_list_begin.next = &new_polygon->vertex_list_end;
    new_polygon->vertex_list_end  .prev = &new_polygon->vertex_list_begin;
    new_polygon->vertex_list_end  .next = NULL;

    new_polygon->edge_list_begin.prev   = NULL;
    new_polygon->edge_list_begin.next   = &new_polygon->edge_list_end;
    new_polygon->edge_list_end  .prev   = &new_polygon->edge_list_begin;
    new_polygon->edge_list_end  .next   = NULL;

    new_polygon->vertex_class_begin.prev    = NULL;
    new_polygon->vertex_class_begin.next    = &new_polygon->vertex_class_end;
    new_polygon->vertex_class_end  .prev    = &new_polygon->vertex_class_begin;
    new_polygon->vertex_class_end  .next    = NULL;

    new_polygon->edge_class_begin.prev  = NULL;
    new_polygon->edge_class_begin.next  = &new_polygon->edge_class_end;
    new_polygon->edge_class_end  .prev  = &new_polygon->edge_class_begin;
    new_polygon->edge_class_end  .next  = NULL;

    return new_polygon;
}


static void make_square(
    WE2DPolygon *polygon)
{


    int         i,
		j,
		k;
    WE2DVertex  *initial_vertices[4];
    WE2DEdge        *initial_edges[4];

    const static int evdata[4][2] =
	{
	    {0, 1},
	    {2, 3},
	    {0, 2},
	    {1, 3}
	    };
    const static int eedata[4][2] =
	{
	    { 2,  3}, { 2,  3},
	    { 0, 1}, { 0, 1}
	
	};

    for (i = 0; i < 4; i++)
    {
	initial_vertices[i] = NEW_STRUCT(WE2DVertex);
	INSERT_BEFORE(initial_vertices[i], &polygon->vertex_list_end);
    }

    for (i = 0; i < 4; i++)
    {
	initial_edges[i] = NEW_STRUCT(WE2DEdge);
	INSERT_BEFORE(initial_edges[i], &polygon->edge_list_end);
    }

    for (i = 0; i < 4; i++)
    {
	initial_vertices[i]->x[0] = ONE;
	initial_vertices[i]->x[1] = (i & 2) ? SQUARE_SIZE : -SQUARE_SIZE;
	initial_vertices[i]->x[2] = (i & 1) ? SQUARE_SIZE : -SQUARE_SIZE;
	initial_vertices[i]->x[3] = ZERO;
    }

    for (i = 0; i < 4; i++)
    {
	for (j = 0; j < 2; j++)     /*  j = tail, tip   */
	initial_edges[i]->v[j] = initial_vertices[evdata[i][j]];

	for (j = 0; j < 2; j++)     /*  j = tail, tip   */
	
	initial_edges[i]->e[j] = initial_edges[eedata[i][j]];

	   initial_edges[i]->mate =NULL;
	   initial_edges[i]->group_element = NULL;
	initial_edges[i]->clean = FALSE;
    }

}


static FuncResult slice_polygon(
    WE2DPolygon *polygon,
    AlgMatrixPairList   *gen_list)
{
    AlgMatrixPair   *matrix_pair;
    for (matrix_pair = gen_list->begin.next;
     matrix_pair != &gen_list->end;
     matrix_pair = matrix_pair->next)
    
      if (!matrix_pair->m[0]->matrix.is_identity(MATRIX_EPSILON))
	if (intersect_with_halfplanes(polygon, matrix_pair) == func_failed)
	  return func_failed;
    return func_OK;
      }


static FuncResult intersect_with_halfplanes(
    WE2DPolygon *polygon,
    AlgMatrixPair   *matrix_pair)
{
  int       i;
  WE2DEdge  *new_edge[2];

  if (matrix_pair->m[0]->matrix.deviation() > DEVIATION_EPSILON)
    {printf("too much deviation from o31\n");return func_failed;}
  if (same_image_of_origin(matrix_pair->m[0]->matrix, matrix_pair->m[1]->matrix) == TRUE)
    {
      if (close(matrix_pair->m[0]->matrix, matrix_pair->m[1]->matrix, MATRIX_EPSILON) == FALSE)
	uFatalError("intersect_with_halfplanes", "Dirichlet_construction");

      if (slice_with_line(polygon, *matrix_pair->m[0], &new_edge[0]) == func_failed)
		  return func_failed;

if (new_edge[0] != NULL)
  new_edge[0]->mate = new_edge[0];
return func_OK;
}
  for (i = 0; i < 2; i++)
    if (slice_with_line(polygon, *matrix_pair->m[i], &new_edge[i]) == func_failed)
      return func_failed;
  if (new_edge[0] != NULL && edge_still_exists(polygon, new_edge[0]) == FALSE)
    new_edge[0] = NULL;
  for (i = 0; i < 2; i++)
    if (new_edge[i] != NULL)
      new_edge[i]->mate = new_edge[!i];
  return func_OK;
}


static Boolean same_image_of_origin(
		    O31_matrix const&   m0,
		    O31_matrix const&   m1)
{
  int   i;
  
  for (i = 0; i < 4; i++)
    if (fabs(m0(i,0) - m1(i,0)) > MATRIX_EPSILON)
      return FALSE;
  
  return TRUE;
}



static FuncResult slice_with_line(
		  WE2DPolygon   *polygon,
		  AlgMatrix     m,
		  WE2DEdge  **new_edge)
{       
  O31_vector    normal_vector;
  (*new_edge)= NULL;
  if (compute_normal_to_Dirichlet_plane(m.matrix, normal_vector) == func_failed)
    return func_failed;
  compute_vertex_to_line_distances(polygon, normal_vector);
  if (positive_vertices_exist(polygon) == FALSE)
    return func_OK;
  cut_edges(polygon);
  if (create_edge(polygon,new_edge) == func_failed)
    return func_failed; 
  if (check_topology_of_cut(polygon) == func_failed)
    return func_failed;
  (*new_edge)->mate  = NULL;
  (*new_edge)->clean = FALSE;
  (*new_edge)->group_element = NEW_ARRAY(1,AlgMatrix );
  Alg_copy(&((*new_edge)->group_element), &m);
  install_new_edge(polygon, *new_edge);
  return func_OK;
}


static FuncResult compute_normal_to_Dirichlet_plane(
			    O31_matrix const&   m,
			    O31_vector& normal_vector)
{
  int       i;
  double    max_abs;
  
  for (i = 0; i < 4; i++)
    normal_vector[i] = m(i,0);
  
  normal_vector[0] -= ONE;
  
  max_abs = ZERO;
  for (i = 0; i < 4; i++)
    if (fabs(normal_vector[i]) > max_abs)
      max_abs = fabs(normal_vector[i]);
  
    if (max_abs < ERROR_EPSILON)
     {printf("error in computing normal\n"); return func_failed;}
    for (i = 0; i < 4; i++)
      normal_vector[i] /= max_abs;
  return func_OK;
}


static void compute_vertex_to_line_distances(
			 WE2DPolygon    *polygon,
    O31_vector const&   normal_vector)
{
  WE2DVertex    *vertex;

  
  for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next)
    {
      vertex->distance_to_plane = o31_inner_product(vertex->x, normal_vector);
      
      if (vertex->distance_to_plane > polygon->vertex_epsilon)
    vertex->which_side_of_plane = +1;
      else if (vertex->distance_to_plane < - polygon->vertex_epsilon)
    vertex->which_side_of_plane = -1;
      else
    vertex->which_side_of_plane = 0;
    }
}


static Boolean positive_vertices_exist(
    WE2DPolygon *polygon)
{
    WE2DVertex  *vertex;
    Boolean     positive_vertices_exist,
		negative_vertices_exist;

    positive_vertices_exist = FALSE;
    negative_vertices_exist = FALSE;

    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)
    {
	if (vertex->which_side_of_plane == +1)
	    positive_vertices_exist = TRUE;

	if (vertex->which_side_of_plane == -1)
	    negative_vertices_exist = TRUE;
    }

    if (negative_vertices_exist == FALSE)
	uFatalError("positive_vertices_exist", "Dirichlet2D_construction");

    return positive_vertices_exist;
}


static void cut_edges(
	      WE2DPolygon   *polygon)
{
  WE2DEdge      *edge;
  int       i, j;
  double        t;
  O31_vector cut_point;
  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)
    
    for (i = 0; i < 2; i++)
      
      if (edge->v[ i]->which_side_of_plane == -1
      && edge->v[!i]->which_side_of_plane == +1)
    {
      t = ( ZERO  - edge->v[tail]->distance_to_plane) 
	/(edge->v[tip]->distance_to_plane - edge->v[tail]->distance_to_plane);
      for (j = 0; j < 4; j++)
	cut_point[j] = (ONE - t) * edge->v[tail]->x[j]  +  t * edge->v[tip]->x[j];
      split_edge(edge, cut_point, TRUE);
    }
}


static void split_edge(
	WE2DEdge    *old_edge,
	O31_vector const&   cut_point,
	Boolean     set_Dirichlet_construction_fields)
{
  WE2DEdge  *new_edge,
  *left_neighbor,
  *right_neighbor;
  WE2DVertex    *new_vertex;
   new_edge = NEW_STRUCT(WE2DEdge);
  new_vertex    = NEW_STRUCT(WE2DVertex);
  if(old_edge->v[tail]->which_side_of_plane == +1)
    {
      new_edge->v[tail] = old_edge->v[tail];
      new_edge->v[tip]  = new_vertex;
      
      new_edge->e[tail]  = old_edge->e[tail];
      new_edge->e[tip ]  = old_edge;
      old_edge->v[tail] = new_vertex;
      old_edge->e[tail] = new_edge;
      left_neighbor = new_edge->e[tail];
      if (left_neighbor->v[tip] == new_edge->v[tail])
    left_neighbor->e[tip] = new_edge;
      else
    if (left_neighbor->v[tail] == new_edge->v[tail])
      left_neighbor->e[tail] = new_edge;
    else
      uFatalError("split_edge", "Dirichlet_construction");
    }
  else 
    {
      new_edge->v[tip] = old_edge->v[tip];
      new_edge->v[tail] = new_vertex;
      new_edge->e[tip] = old_edge->e[tip];
      new_edge->e[tail]= old_edge;
      old_edge->v[tip] = new_vertex;
      old_edge->e[tip] = new_edge;
      left_neighbor = new_edge->e[tip];
      if(left_neighbor->v[tail]==new_edge->v[tip])
    left_neighbor->e[tail]=new_edge;
      else
    if(left_neighbor->v[tip]==new_edge->v[tip])
      left_neighbor->e[tip]=new_edge;
    else
      uFatalError("split_edge","Dirichlet_construction");
    }
 
   new_edge->mate =  NULL;
    if(old_edge->mate != NULL)
    {
	    old_edge->mate->clean=FALSE;
    }
 
  new_edge->group_element = NULL;
  
  o31_copy_vector(new_vertex->x, cut_point);
  
  
  new_vertex->distance_to_plane = ZERO;
  new_vertex->which_side_of_plane= 0;

  
  INSERT_BEFORE(new_edge, old_edge);
  
  INSERT_BEFORE(new_vertex, new_edge->e[tip]->v[tip]);

  
  
}


static FuncResult create_edge(
    WE2DPolygon *polygon,
	WE2DEdge **new_edge)
{
    
    int     i,
	    count;
    WE2DEdge    *edge,
	    *edge_before_vertex[2],
	    *edge_after_vertex[2],
	    *temp_edge;

    
    make_samedirect(polygon);
	    
    count = 0;
    for(edge = polygon->edge_list_begin.next;
	edge !=  &polygon->edge_list_end;
	edge = edge->next)
      {
	if (edge->v[tip]->which_side_of_plane == 0)
	  {
	if (count == 2)
	 {printf("problem in create edge\n"); return func_failed;}
	if(edge->v[tail]->which_side_of_plane == -1)
	  edge_before_vertex[0] = edge;
	else
	  edge_before_vertex[1] = edge;
	count++;
	  }
      } 
    
    if (count < 2)
      return func_OK;
    
    for (i = 0; i < 2; i++)
      edge_after_vertex[i] = edge_before_vertex[i]->e[tip];
    
    for (i = 0; i < 2; i++)
      if (edge_after_vertex[i] == edge_before_vertex[!i])
	return func_OK;
    
    
    if (edge_before_vertex[0]->v[tail]->which_side_of_plane == -1
	&& edge_after_vertex [0]->v[tip] ->which_side_of_plane == +1)
      {
	
      }
    else
      if (edge_before_vertex[0]->v[tail]->which_side_of_plane == +1
	  && edge_after_vertex [0]->v[tip] ->which_side_of_plane == -1)
	{
	  
	  temp_edge     = edge_before_vertex[0];
	  edge_before_vertex[0] = edge_before_vertex[1];
	  edge_before_vertex[1] = temp_edge;
	  
	  temp_edge     = edge_after_vertex[0];
	  edge_after_vertex[0]  = edge_after_vertex[1];
	  edge_after_vertex[1]  = temp_edge;
	}
    else 
      return func_failed;
    (*new_edge) = NEW_STRUCT(WE2DEdge);
    if ((*new_edge)== NULL)
      printf("mem not allocated\n\n");
    
    (*new_edge)->v[tail]    = edge_before_vertex[0]->v[tip];
    (*new_edge)->v[tip] = edge_before_vertex[1]->v[tip];
    
    (*new_edge)->e[tail]        = edge_before_vertex[0];
    (*new_edge)->e[tip]         = edge_after_vertex[1];
    
    edge_before_vertex[0]->e[tip]   = (*new_edge);
    edge_after_vertex [1]->e[tail]  = (*new_edge);
    
    INSERT_BEFORE((*new_edge), edge_after_vertex[1]);
    return func_OK;
      }


void make_samedirect(
	     WE2DPolygon    *polygon)
{
  WE2DEdge  *edge,*edge1;
  edge1=  edge = polygon->edge_list_begin.next;
  do
    {
      if (edge->v[tip] != edge->e[tip]->v[tail]) 
    redirect2D_edge(edge->e[tip],FALSE);
      edge = edge->e[tip];
    }
  while (edge1 !=edge);
    
}


void redirect2D_edge(
    WE2DEdge    *edge,
    Boolean redirect_neighbor_fields)
{
  WE2DVertex    *temp_vertex;
  WE2DEdge      *temp_edge;

  temp_vertex       = edge->v[tail];
  edge->v[tail] = edge->v[tip];
  edge->v[tip]  = temp_vertex;


  temp_edge         = edge->e[tail];
  edge->e[tail]         = edge->e[tip];
  edge->e[tip]          = temp_edge;



  /*    if (redirect_neighbor_fields)
    {
    int         i;
    WE2DEdge        *nbr_edge;
    WE2DEdgeSide    side,
		    nbr_side;
    WE2DEdge    *temp_edge;
    Boolean     temp_boolean;

    for (side = 0; side < 2; side++)
    {
    nbr_edge = edge->neighbor[side];
    nbr_side = (edge->preserves_sides[side] ? side : !side);
    nbr_edge->preserves_sides[nbr_side]
		  = ! nbr_edge->preserves_sides[nbr_side];
	nbr_edge->preserves_direction[nbr_side]
		= ! nbr_edge->preserves_direction[nbr_side];
    }



    for (side = 0; side < 2; side++)
    {
    edge->preserves_sides[side]     = ! edge->preserves_sides[side];
    edge->preserves_direction[side] = ! edge->preserves_direction[side];
    }

    temp_edge               = edge->neighbor[left];
    edge->neighbor[left]    = edge->neighbor[right];
    edge->neighbor[right]   = temp_edge;

    temp_boolean                = edge->preserves_sides[left];
    edge->preserves_sides[left]     = edge->preserves_sides[right];
    edge->preserves_sides[right]        = temp_boolean;

    temp_boolean                = edge->preserves_direction[left];
    edge->preserves_direction[left]     = edge->preserves_direction[right];
    edge->preserves_direction[right]    = temp_boolean;

    temp_boolean                = edge->preserves_orientation[left];
    edge->preserves_orientation[left]   = edge->preserves_orientation[right];
    edge->preserves_orientation[right]  = temp_boolean;
    }*/
}


static FuncResult check_topology_of_cut(
    WE2DPolygon *polygon)
{
  int       num_zero_edges,
	count;
  WE2DVertex    *vertex,
	*tip_vertex;
  WE2DEdge  *edge,
	*starting_edge;

  for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next)
    vertex->zero_order = 0;

  num_zero_edges = 0;

  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)

    if (edge->v[tail]->which_side_of_plane == 0
    && edge->v[tip] ->which_side_of_plane == 0)
      {
    edge->v[tail]->zero_order++;
    edge->v[tip] ->zero_order++;
    num_zero_edges++;
      }

 /* for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next)
    
    if (vertex->which_side_of_plane == 0
    && vertex->zero_order != 1)
      return func_failed; */ 

  starting_edge = NULL;

  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)

    if (edge->v[tail]->which_side_of_plane == 0
    && edge->v[tip] ->which_side_of_plane == 0)
      {
    starting_edge = edge;
    break;
      }

  if (starting_edge == NULL)
    uFatalError("check_topology_of_cut", "Dirichlet_construction");


}


static void install_new_edge(
    WE2DPolygon *polygon,
    WE2DEdge        *new_edge)
{
  WE2DEdge      *edge;
  WE2DVertex            *vertex;
  WE2DEdge      *dead_edge;
  WE2DVertex    *dead_vertex;
  WE2DVertex    *tip_vertex;



  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)
    
    if (edge->v[tail]->which_side_of_plane == +1
    || edge->v[tip] ->which_side_of_plane == +1)
      {
    if (edge->v[tail]->which_side_of_plane == 0)
      redirect2D_edge(edge, FALSE);
    if(edge->mate !=NULL)
      {
	my_free_array(edge->group_element);
	edge->group_element = NULL;
	edge->mate->mate = NULL;
      }
    if (edge != NULL)
      if(edge->group_element != NULL)
	{
	  my_free_array(edge->group_element);
	  edge->group_element = NULL;
	}
    dead_edge = edge;
    edge = edge->prev;
    REMOVE_NODE(dead_edge);
    my_free(dead_edge);
      }
    
  for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next)
      
    if (vertex->which_side_of_plane == +1)
      {
    dead_vertex = vertex;
    vertex = vertex->prev;
    REMOVE_NODE(dead_vertex);
	  
    my_free(dead_vertex);
      }

}


static Boolean edge_still_exists(
    WE2DPolygon *polygon,
    WE2DEdge        *edge0)
{
    
  WE2DEdge  *edge;

  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)

    if (edge == edge0)
      return TRUE;
  return FALSE;
}


static Boolean has_hyperideal_vertices(
    WE2DPolygon *polygon)
{
  WE2DVertex    *vertex;

  for (vertex = polygon->vertex_list_begin.next;
       vertex != &polygon->vertex_list_end;
       vertex = vertex->next)

    if (o31_inner_product(vertex->x, vertex->x) > HYPERIDEAL_EPSILON)
      return TRUE;

  return FALSE;
}


static void compute_all_products(
    WE2DPolygon *polygon,
    AlgMatrixPairList   *product_list)
{
    AlgMatrixPairList   current_list;
    AlgMatrixPair   *product_tree;

    poly_to_current_list(polygon, &current_list);
    current_list_to_product_tree(&current_list, &product_tree);
    product_tree_to_product_list(product_tree, product_list);
    free_matrix_pairs(&current_list);
}


static void poly_to_current_list(
    WE2DPolygon *polygon,
    AlgMatrixPairList   *current_list)
{
    WE2DEdge        *edge;
    AlgMatrixPair   *matrix_pair;

    current_list->begin.prev = NULL;
    current_list->begin.next = &current_list->end;
    current_list->end  .prev = &current_list->begin;
    current_list->end  .next = NULL;

    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

	edge->copied = FALSE;



    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

	if (edge->group_element != NULL  &&  edge->copied == FALSE)
	{
	    matrix_pair = NEW_STRUCT(AlgMatrixPair);
	    Alg_copy(&matrix_pair->m[0],edge->group_element);
	    Alg_invert(matrix_pair->m[0], &matrix_pair->m[1]);
	    matrix_pair->height = matrix_pair->m[0]->matrix(0,0);
	    INSERT_BEFORE(matrix_pair, &current_list->end);

	    edge->copied = TRUE;
	    if (edge->mate != NULL)
	    edge->mate->copied = TRUE;
	}
}


static void current_list_to_product_tree(
    AlgMatrixPairList   *current_list,
    AlgMatrixPair   **product_tree)
{
  AlgMatrixPair *matrix_pair_a,
	*matrix_pair_b;
  int       i,
		j;
  AlgMatrix *product;


  *product_tree = NULL;

  for (matrix_pair_a = current_list->begin.next;
       matrix_pair_a != &current_list->end;
       matrix_pair_a = matrix_pair_a->next)
    
    for (matrix_pair_b = matrix_pair_a;
     matrix_pair_b != &current_list->end;
     matrix_pair_b = matrix_pair_b->next)
      
      for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
      {
	precise_Alg_product(matrix_pair_a->m[i], matrix_pair_b->m[j],&product);
	if (already_on_product_tree(product, *product_tree) == FALSE)
	  add_to_product_tree(product, product_tree);
      }
}


static Boolean already_on_product_tree(
    AlgMatrix   *product,
    AlgMatrixPair   *product_tree)
{
  AlgMatrixPair *subtree_stack,
	*subtree;
  int       i;

  subtree_stack = product_tree;
  if (product_tree != NULL)
    product_tree->next_subtree = NULL;

  
  while (subtree_stack != NULL)
    {
      subtree           = subtree_stack;
      subtree_stack     = subtree_stack->next_subtree;
      subtree->next_subtree = NULL;

      if (subtree->left_child != NULL
      && product->matrix(0,0) < subtree->height + MATRIX_EPSILON)
    {
      subtree->left_child->next_subtree = subtree_stack;
      subtree_stack = subtree->left_child;
    }
      if (subtree->right_child != NULL
      && product->matrix(0,0) > subtree->height - MATRIX_EPSILON)
    {
      subtree->right_child->next_subtree = subtree_stack;
      subtree_stack = subtree->right_child;
    }

      for (i = 0; i < 2; i++)
    if (Alg_equal(product, subtree->m[i], MATRIX_EPSILON) == TRUE)
      return TRUE;
    }

  return FALSE;
}


static void add_to_product_tree(
    AlgMatrix   *product,
    AlgMatrixPair   **product_tree)
{
    AlgMatrixPair   **home;
    double      product_height;

    product_height = product->matrix(0,0);

    home = product_tree;

    while (*home != NULL)
    {
	if (product_height < (*home)->height)
	    home = &(*home)->left_child;
	else
	    home = &(*home)->right_child;
    }

    (*home) = NEW_STRUCT(AlgMatrixPair);
    Alg_copy(&(*home)->m[0], product);
    Alg_invert((*home)->m[0],&(*home)->m[1]);
    (*home)->height = (*home)->m[0]->matrix(0,0);
    (*home)->left_child     = NULL;
    (*home)->right_child    = NULL;
    (*home)->next_subtree   = NULL;
    (*home)->prev           = NULL;
    (*home)->next           = NULL;
}


static void product_tree_to_product_list(
    AlgMatrixPair   *product_tree,
    AlgMatrixPairList   *product_list)
{
  product_list->begin.prev = NULL;
  product_list->begin.next = &product_list->end;
  product_list->end  .prev = &product_list->begin;
  product_list->end  .next = NULL;

  append_tree_to_list(product_tree, &product_list->end);
}


static void append_tree_to_list(
    AlgMatrixPair   *product_tree,
    AlgMatrixPair   *list_end)
{
  AlgMatrixPair *subtree_stack,
	*subtree;

  subtree_stack = product_tree;
  if (product_tree != NULL)
    product_tree->next_subtree = NULL;

  while (subtree_stack != NULL)
    {
      subtree       = subtree_stack;
      subtree_stack = subtree_stack->next_subtree;
      subtree->next_subtree = NULL;
    
      if (subtree->left_child == NULL  &&  subtree->right_child == NULL)
    {
      INSERT_BEFORE(subtree, list_end);
    }
      else
    {
      if (subtree->right_child != NULL)
	{
	  subtree->right_child->next_subtree = subtree_stack;
	  subtree_stack = subtree->right_child;
	  subtree->right_child = NULL;
	}

      subtree->next_subtree = subtree_stack;
      subtree_stack = subtree;

      if (subtree->left_child != NULL)
	{
	  subtree->left_child->next_subtree = subtree_stack;
	  subtree_stack = subtree->left_child;
	  subtree->left_child = NULL;
	}
    }
    }
}


static FuncResult check_edges(
    WE2DPolygon *polygon)
{
  WE2DEdge      *edge;
  Boolean       edge_was_pared;


  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)

    edge->clean = FALSE;

  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)
    if (edge->clean == FALSE)
      {
    if (pare_edge(edge, polygon, &edge_was_pared) == func_failed)
      { return func_failed;}

    if (edge_was_pared == TRUE)
      edge = &polygon->edge_list_begin;
      }

  return func_OK;
}


static FuncResult pare_edge(
    WE2DEdge    *edge,
    WE2DPolygon *polygon,
    Boolean     *edge_was_pared)
{
  if (edge->mate != NULL)
    return pare_mated_edge(edge, polygon, edge_was_pared);
  else
    return pare_mateless_edge(edge, polygon, edge_was_pared);
}


static FuncResult pare_mated_edge(
    WE2DEdge    *edge,
    WE2DPolygon *polygon,
    Boolean     *edge_was_pared)
{
	int             i;
    WE2DEdge    *n_edge;
    AlgMatrix   *alpha;

	n_edge = edge->mate;
	for(i=0;i<2;i++)
    {
	if (i==0)
	    alpha = n_edge->e[tail]->group_element;
	else
	    alpha = n_edge->e[tip] ->group_element;

    
	if (try_this_alpha(alpha, edge, polygon, edge_was_pared) == func_failed)
	    return func_failed;

	if (*edge_was_pared == TRUE)
	    return func_OK;

	} 
    edge->clean = TRUE;

    *edge_was_pared = FALSE;
	return func_OK;
}


static FuncResult pare_mateless_edge(
    WE2DEdge            *edge,
    WE2DPolygon *polygon,
    Boolean         *edge_was_pared)
{

    WE2DEdge        *edge1;
    AlgMatrix   *alpha;
	for (edge1 = polygon->edge_list_begin.next;
	 edge1 != &polygon->edge_list_end;
	 edge1 = edge1->next)
    {
	alpha = edge1->group_element;

	if (try_this_alpha(alpha, edge, polygon, edge_was_pared) == func_failed)
	  return func_failed;
			 
	if (*edge_was_pared == TRUE)
	    return func_OK;
    }

    {printf(" could not find mate\n");return func_failed;}
}


static FuncResult try_this_alpha(
    AlgMatrix       *alpha,
    WE2DEdge        *edge,
    WE2DPolygon         *polygon,
    Boolean         *edge_was_pared)
{
    WE2DVertex  *vertex;
    WE2DEdge    edge1;
    AlgMatrix   *beta;
    O31_vector  normal;
    AlgMatrixPair   matrix_pair;
	int             i;

    precise_Alg_product(edge->group_element, alpha, &beta);
    compute_normal_to_Dirichlet_plane(beta->matrix, normal);
	for(i=0;i<2;i++)
    {
    
	if (i == 0)
	    vertex = edge->v[tip];
	else
	    vertex = edge->v[tail];

    
	if (o31_inner_product(vertex->x, normal) > polygon->vertex_epsilon)
	{
	    Alg_copy(&(matrix_pair.m[0]), beta);
	    Alg_invert(matrix_pair.m[0],&(matrix_pair.m[1]));

	    matrix_pair.height  = ZERO;
	    matrix_pair.prev    = NULL;
	    matrix_pair.next    = NULL;

	if (intersect_with_halfplanes(polygon, &matrix_pair) == func_failed)
		  return func_failed;

		*edge_was_pared = TRUE;

	    return func_OK;
	}
    }
    *edge_was_pared = FALSE;
    return func_OK;
}


void count_cells(
    WE2DPolygon *polygon)
{
    WE2DVertex  *vertex;
    WE2DEdge        *edge;


    polygon->num_vertices   = 0;
    polygon->num_edges      = 0;

    for (vertex = polygon->vertex_list_begin.next;
	 vertex != &polygon->vertex_list_end;
	 vertex = vertex->next)

    polygon->num_vertices++;

    for (edge = polygon->edge_list_begin.next;
	 edge != &polygon->edge_list_end;
	 edge = edge->next)

	polygon->num_edges++;

}


static void sort_edges(
    WE2DPolygon *polygon)
{
    WE2DEdge    **array,
	*edge;
    int     i;

    if (polygon->num_edges < 3)
	uFatalError("sort_edgxes", "Dirichlet_construction");

    array = NEW_ARRAY(polygon->num_edges, WE2DEdge *);

    for (edge = polygon->edge_list_begin.next,
	    i = 0;
	 edge != &polygon->edge_list_end;
	 edge = edge->next,
	    i++)

	array[i] = edge;

    if (i != polygon->num_edges)
	uFatalError("sort_edges", "Dirichlet_construction");

    qsort(  array,
	    polygon->num_edges,
	    sizeof(WE2DEdge *),
	    compare_edge_distance);

    polygon->edge_list_begin.next = array[0];
    array[0]->prev = &polygon->edge_list_begin;
    array[0]->next = array[1];

    for (i = 1; i < polygon->num_edges - 1; i++)
    {
	array[i]->prev = array[i-1];
	array[i]->next = array[i+1];
    }

    array[polygon->num_edges - 1]->prev = array[polygon->num_edges - 2];
    array[polygon->num_edges - 1]->next = &polygon->edge_list_end;
    polygon->edge_list_end.prev = array[polygon->num_edges - 1];

    my_free_array(array); 
}


static int CDECL compare_edge_distance(
    const void  *ptr1,
    const void  *ptr2)
{
    double  diff;

    diff = (*(*((WE2DEdge **)ptr1))->group_element).matrix(0,0)
	 - (*(*((WE2DEdge **)ptr2))->group_element).matrix(0,0);

    if (diff < ZERO)
	return -1;
    if (diff > ZERO)
	return +1;
    return 0;
}



static FuncResult verify_group(
    WE2DPolygon *polygon,
    AlgMatrixPairList   *gen_list)
{

  AlgMatrixPair *matrix_pair;
  AlgMatrix *m,
		*candidate;
  Boolean       progress;
  WE2DEdge      *edge;
  double        verify_epsilon;

  for (matrix_pair = gen_list->begin.next;
       matrix_pair != &gen_list->end;
       matrix_pair = matrix_pair->next)
    {
      Alg_copy(&m,matrix_pair->m[0]);

      verify_epsilon = VERIFY_EPSILON;

      while (!m->matrix.is_identity(MATRIX_EPSILON))
    {
      progress = FALSE;

      for (edge = polygon->edge_list_begin.next;
	   edge != &polygon->edge_list_end;
	   edge = edge->next)
	{
	  Alg_product(m, edge->group_element,&candidate);

	  if (m->matrix(0,0) - candidate->matrix(0,0) > verify_epsilon)
	{
	  Alg_copy(&m,candidate);
	  progress  = TRUE;
	  break;
	}
	}

      if (progress == FALSE)
	{
	  if (verify_epsilon > ZERO)
	verify_epsilon = ZERO;
	  else
	{
	  uAcknowledge("Please tell Jeff Weeks that SnapPea seems to have computed a Dirichlet domain for a finite-sheeted cover of the manifold/orbifold.");
	  return func_failed;
	}
	}
    }
    }

  return func_OK;
}


static void rewrite_gen_list(
    WE2DPolygon *polygon,
    AlgMatrixPairList   *gen_list)
{
  WE2DEdge  *edge,
		*mate;
  AlgMatrixPair *new_matrix_pair;

  AlgMatrix    Alg_identity;
  Alg_identity.matrix = O31_identity;
  Alg_identity.identity = TRUE;

  free_matrix_pairs(gen_list); 

  new_matrix_pair = NEW_STRUCT(AlgMatrixPair);
  Alg_copy(&new_matrix_pair->m[0], &Alg_identity);
  Alg_copy(&new_matrix_pair->m[1], &Alg_identity);
  new_matrix_pair->height = ONE;
  INSERT_BEFORE(new_matrix_pair, &gen_list->end);


  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)
    edge->copied = FALSE;

  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)
    
    if (edge->copied == FALSE)
      {
    mate = edge->mate;
    new_matrix_pair = NEW_STRUCT(AlgMatrixPair);
    Alg_copy(&new_matrix_pair->m[0], edge->group_element);
    Alg_copy(&new_matrix_pair->m[1], mate->group_element);
    new_matrix_pair->height = edge->group_element->matrix(0,0);
    INSERT_BEFORE(new_matrix_pair, &gen_list->end);
    edge->copied = TRUE;
    mate->copied = TRUE;
      }
}


void print_pol_edge(WE2DPolygon *polygon)
{
  WE2DEdge  *edge;
  for (edge = polygon->edge_list_begin.next;
       edge != &polygon->edge_list_end;
       edge = edge->next)
    print_edge(edge);
}

void print_edge(WE2DEdge *edge)
  {  
    printf("tip is   ");
     print_vector(edge->v[tip]->x);
     printf("\n");
     print_vector(edge->v[tail]->x);
     printf("\n");
  }
 


















