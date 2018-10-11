/*
 *  Dirichlet2D.c
 */
#include "Dirichlet2D_prototypes.h"

#define SIMPLIFY_EPSILON        1e-2

#define LENGTH_EPSILON          1e-2

#define FIXED_BASEPOINT_EPSILON 1e-6

static WE2DPolygon  *Dirichlet2D_with_displacement(Triangulation *manifold, double displacement[3], double vertex_epsilon, Boolean centroid_at_origin, Boolean interactive);
static WE2DPolygon  *Dirichlet2D_from_generators_with_displacement(AlgMatrix * generators[], int num_generators, double displacement[3], double vertex_epsilon, Boolean interactive);
static void     array_to_matrix_pair_list(AlgMatrix *generators[], int num_generators, AlgMatrixPairList *gen_list);
static Boolean      is_matrix_on_list(AlgMatrix* m, AlgMatrixPairList *gen_list);
static void     insert_matrix_on_list(AlgMatrix *m, AlgMatrixPairList *gen_list);
static void     simplify_generators(AlgMatrixPairList *gen_list);
static Boolean      generator_fixes_basepoint(AlgMatrixPairList *gen_list);
static double       product_height(O31_matrix const& a, O31_matrix const& b);
static void     generators_from_polygon(WE2DPolygon *polygon, AlgMatrix **generators, int *num_generators);
extern void free_cyclic_word(CyclicWord *);
static double   null_displacement[3] = {0.0, 0.0, 0.0};


WE2DPolygon *Dirichlet2D(
    Triangulation   *manifold,
    double          vertex_epsilon,
    Boolean         centroid_at_origin,
    Boolean         interactive)
{
    return(Dirichlet2D_with_displacement(   manifold,
			null_displacement,
			vertex_epsilon,
				centroid_at_origin,
			interactive));
}


WE2DPolygon *Dirichlet2D_from_generators(
    AlgMatrix   *generators[],
    int     num_generators,
    double      vertex_epsilon,
    Boolean     interactive)
{
    return(Dirichlet2D_from_generators_with_displacement(   generators,
							num_generators,
							null_displacement,
							vertex_epsilon,
							interactive));
}

static void Moebius_array_to_Alg_array(MoebiusTransformation mob[], AlgMatrix* alg[], int n)
{
  int i;
  for (i=0; i<n; i++) {
    alg[i] = NEW_STRUCT(AlgMatrix);
    set_name(&(alg[i]->name),(i+1));
    alg[i]->identity =  FALSE;
    alg[i]->matrix = mob[i]; 
  }
}

static void free_AlgMatrix_array(AlgMatrix* alg[], int n)
{
  int i;
  for (i=0; i<n; i++)
    my_free(alg[i]); 

  my_free_array(alg); 
}

static WE2DPolygon *Dirichlet2D_with_displacement(
    Triangulation   *manifold,
    double          displacement[3],
    double          vertex_epsilon,
    Boolean         centroid_at_origin,
    Boolean         interactive)
{
    MoebiusTransformation   *Moebius_generators;
    AlgMatrix           **Alg_generators;
    WE2DPolygon         *polygon;

    choose_generators(manifold, FALSE, FALSE);
    Moebius_generators  = NEW_ARRAY(manifold->num_generators, MoebiusTransformation);
    Alg_generators      = NEW_ARRAY(manifold->num_generators, AlgMatrix*);
    matrix_generators(manifold, Moebius_generators, centroid_at_origin, TRUE);
    Moebius_array_to_Alg_array(Moebius_generators, Alg_generators, manifold->num_generators);

    polygon = Dirichlet2D_from_generators_with_displacement(
			    Alg_generators,
			    manifold->num_generators,
			    displacement,
			    vertex_epsilon,
			    interactive);

    my_free_array(Moebius_generators);
    free_AlgMatrix_array(Alg_generators, manifold->num_generators);

    return polygon;
}


static WE2DPolygon *Dirichlet2D_from_generators_with_displacement(
    AlgMatrix   *generators[],
    int     num_generators,
    double      displacement[3],
    double      vertex_epsilon,
    Boolean     interactive)
{
    AlgMatrixPair  *matrix_pair;
    int a,b,c;
    AlgMatrixPairList   gen_list;
    WE2DPolygon *polygon;
    Boolean     basepoint_moved;
    double      small_displacement[3] = {0.01734, 0.02035, 0.00000};
/*  printf("%d\n",num_generators);
    for(a=0;a<num_generators;a++)
      for(b=0;b<4;b++)
	for(c=0;c<4;c++)
	  printf("%lf  ",generators[a]->matrix[b][c]);*/
	array_to_matrix_pair_list(generators, num_generators, &gen_list);

    precise_generators(&gen_list); 

    conjugate_matrices(&gen_list, displacement);

    simplify_generators(&gen_list);

    if (generator_fixes_basepoint(&gen_list) == TRUE)
	conjugate_matrices(&gen_list, small_displacement);
	while (TRUE)
    { 
	polygon = compute_Dirichlet2D_domain(&gen_list, vertex_epsilon);

	if (polygon == NULL)
	{
	    free_matrix_pairs(&gen_list);
	    return NULL;
	}
	
	maximize_injectivity_radius(&gen_list, &basepoint_moved, interactive);
	
	if (basepoint_moved == FALSE)
	{
	  free_matrix_pairs(&gen_list);
	  return polygon;
	  
	}

	free_Dirichlet2D_domain(polygon);
    }


}


static void array_to_matrix_pair_list(
    AlgMatrix       *generators[],
    int             num_generators,
    AlgMatrixPairList   *gen_list)
{
    int         i;
    AlgMatrix           Alg_identity;

    Alg_identity.matrix = O31_matrix(1);
    Alg_identity.identity = TRUE;
      
    gen_list->begin.prev = NULL;
    gen_list->begin.next = &gen_list->end;
    gen_list->end  .prev = &gen_list->begin;
    gen_list->end  .next = NULL;

    insert_matrix_on_list(&Alg_identity, gen_list);

    for (i = 0; i < num_generators; i++)
	if (is_matrix_on_list(generators[i], gen_list) == FALSE)
	    insert_matrix_on_list(generators[i], gen_list);
}


Boolean is_matrix_on_list(
    AlgMatrix       *m,
    AlgMatrixPairList   *gen_list)
{
    AlgMatrixPair   *matrix_pair;
    int         i;

    for (matrix_pair = gen_list->begin.next;
	 matrix_pair != &gen_list->end;
	 matrix_pair = matrix_pair->next)

	for (i = 0; i < 2; i++)

	    if (close(m->matrix, matrix_pair->m[i]->matrix, MATRIX_EPSILON))

		return TRUE;

    return FALSE;
}


void insert_matrix_on_list(
    AlgMatrix       *m,
    AlgMatrixPairList   *gen_list)
{
    AlgMatrix   *m_inverse;
    AlgMatrixPair   *new_matrix_pair,
		*mp;

    Alg_invert(m, &m_inverse);

    new_matrix_pair = NEW_STRUCT(AlgMatrixPair);
    Alg_copy(&new_matrix_pair->m[0], m);
    Alg_copy(&new_matrix_pair->m[1], m_inverse);
    new_matrix_pair->height = m->matrix(0,0);

    mp = gen_list->begin.next;
    while (mp != &gen_list->end  &&  mp->height < new_matrix_pair->height)
	mp = mp->next;

    INSERT_BEFORE(new_matrix_pair, mp);
}


void free_matrix_pairs(
    AlgMatrixPairList   *gen_list)
{
    AlgMatrixPair   *dead_node;

    while (gen_list->begin.next != &gen_list->end)
    {
	dead_node = gen_list->begin.next;
	my_free(dead_node->m[0]);
	my_free(dead_node->m[1]);
	REMOVE_NODE(dead_node);
	my_free(dead_node);
    }
}


void free_Dirichlet2D_domain(
    WE2DPolygon *polygon)
{
    WE2DVertex      *dead_vertex;
    WE2DEdge            *dead_edge;
    WE2DVertexClass *dead_vertex_class;
    WE2DEdgeClass       *dead_edge_class;


    if (polygon == NULL)
	uFatalError("free_Dirichlet2D_domain", "Dirichlet2D");

    while (polygon->vertex_list_begin.next != &polygon->vertex_list_end)
    {
	dead_vertex = polygon->vertex_list_begin.next;
	REMOVE_NODE(dead_vertex);
	my_free(dead_vertex);
    }

    while (polygon->edge_list_begin.next != &polygon->edge_list_end)
    {
	dead_edge = polygon->edge_list_begin.next;
	if(dead_edge->group_element != NULL)
	 {
	free_cyclic_word(dead_edge->group_element->name);
	   my_free_array(dead_edge->group_element);}
	REMOVE_NODE(dead_edge);
	my_free(dead_edge);
    }


    while (polygon->vertex_class_begin.next != &polygon->vertex_class_end)
    {
	dead_vertex_class = polygon->vertex_class_begin.next;
	REMOVE_NODE(dead_vertex_class);
	my_free(dead_vertex_class);
    }

    while (polygon->edge_class_begin.next != &polygon->edge_class_end)
    {
	dead_edge_class = polygon->edge_class_begin.next;
	REMOVE_NODE(dead_edge_class);
	my_free(dead_edge_class);
    }

    my_free(polygon);
}


static void simplify_generators(
    AlgMatrixPairList   *gen_list)
{

    Boolean     progress;
    AlgMatrixPair   *aA,
	    *bB,
	    *best_aA;
    AlgMatrix   **best_aA_factor,
	    **best_bB_factor;
    double      max_improvement,
	    improvement;
    int         i,
		j;


    while (TRUE)
    {
	max_improvement = ZERO;

	for (aA = gen_list->begin.next;
	     aA != &gen_list->end;
	     aA = aA->next)

	for (bB = gen_list->begin.next;
	 bB != &gen_list->end;
	 bB = bB->next)
	{

	    if (aA == bB)
		    continue;
    
	    if (aA->height < bB->height)
	    continue;

	    for (i = 0; i < 2; i++)
	    for (j = 0; j < 2; j++)
	    {
	improvement = aA->height - product_height(aA->m[i]->matrix, bB->m[j]->matrix);
	    if (improvement > max_improvement)
	    {
	    max_improvement = improvement;
	    best_aA = aA;
	    best_aA_factor = &aA->m[i];
		best_bB_factor = &bB->m[j];
		}
	    }
	}

	
	if (max_improvement < SIMPLIFY_EPSILON)
	    break;

	precise_Alg_product(*best_aA_factor, *best_bB_factor, &best_aA->m[0]);
	Alg_invert(best_aA->m[0], &best_aA->m[1]);
	best_aA->height = best_aA->m[0]->matrix(0,0);

	if (best_aA->m[0]->matrix.is_identity(MATRIX_EPSILON))
	{
	    REMOVE_NODE(best_aA);
	    my_free(best_aA);
	}
    }
}


static Boolean generator_fixes_basepoint(
    AlgMatrixPairList   *gen_list)
{
    AlgMatrixPair   *matrix_pair;

    for (matrix_pair = gen_list->begin.next;
     matrix_pair != &gen_list->end;
     matrix_pair = matrix_pair->next)
    if (matrix_pair->m[0]->matrix(0,0) < ONE + FIXED_BASEPOINT_EPSILON)
      if (!matrix_pair->m[0]->matrix.is_identity(MATRIX_EPSILON))
	return TRUE;
    return FALSE;
}


static double product_height(
    O31_matrix const& a,
    O31_matrix const& b)
{

    double  height;
    int     i;

    height = ZERO;

    for (i = 0; i < 4; i++)
	height += a(0,i) * b(i,0);

    return height;
}


void change_basepoint(
    WE2DPolygon **polygon,
    Triangulation   *manifold,
    AlgMatrix       **generators,
    int         num_generators,
    double          displacement[3],
    double          vertex_epsilon,
    Boolean         centroid_at_origin,
    Boolean         interactive)
{
    AlgMatrix   *gen;
    int     num_gen;

    if (*polygon != NULL)
    {
	generators_from_polygon(*polygon, &gen, &num_gen);

	free_Dirichlet2D_domain(*polygon);

	(*polygon) = Dirichlet2D_from_generators_with_displacement(
	    &gen, num_gen, displacement, vertex_epsilon, interactive);

	my_free_array(gen);
    }
    else
    {
	if (manifold != NULL)
	(*polygon) = Dirichlet2D_with_displacement(
	   manifold, displacement, vertex_epsilon, centroid_at_origin, interactive);
	else if (generators != NULL  &&  num_generators > 0)
	(*polygon) = Dirichlet2D_from_generators_with_displacement(
	generators, num_generators, displacement, vertex_epsilon, interactive);
	else
	    uFatalError("change_basepoint", "Dirichlet");
    }
}


static void generators_from_polygon(
    WE2DPolygon *polygon,
    AlgMatrix   **generators,
    int     *num_generators)
{
    WE2DEdge    *edge;
    int i;

    *num_generators = polygon->num_edges;

    *generators = NEW_ARRAY(*num_generators, AlgMatrix );

    for (edge = polygon->edge_list_begin.next,
	    i = 0;
	 edge != &polygon->edge_list_end;
	 edge = edge->next,
	    i++)

	Alg_copy(&(generators)[i], edge->group_element);

    if (i != *num_generators)
	uFatalError("generators_from_polygon", "Dirichlet");
}







