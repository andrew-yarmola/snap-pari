//	quotients.c
//
//	This is a quick hack to create quotients of cusped hyperbolic
//	3-manifolds.  It's being written in a 24-hour time frame,
//	so I've done some ugly things, like copying the contentes of
//	isometry_cusped.c and modifying it to suit local purposes.

#include "kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRUE	1
#define FALSE	0

#ifdef __powerc
#define DOUBLE_FORMAT	"%19.15lf "
#else
#define DOUBLE_FORMAT	"%22.18lf "
#endif

static void					compute_symmetries(Triangulation *aTriangulation, IsometryList **aSymmetryList);
static Complex				compute_cross_ratio(Complex a[4]);

static FuncResult	attempt_isometry(Triangulation *manifold0, Tetrahedron *tet0, Tetrahedron *tet1, Permutation map0);
static Boolean		is_isometry_plausible(Tetrahedron *initial_tet0, Tetrahedron *initial_tet1, Permutation initial_map);
static FuncResult	set_tet_image(Tetrahedron *tet0, Tetrahedron *tet1, Permutation map);
static void			copy_isometry(Triangulation *manifold0, Triangulation *manifold1, Isometry **new_isometry);
static void			compute_cusp_action(Triangulation *manifold0, Triangulation *manifold1, Isometry *isometry);
static void			compute_cusp_image(Triangulation *manifold0, Isometry *isometry);
static void			compute_cusp_map(Triangulation *manifold0, Triangulation *manifold1, Isometry *isometry);
static void			copy_images_to_scratch(Triangulation *manifold, int which_set, Boolean double_copy_on_tori);
static Boolean		does_isometry_extend_to_link(Isometry *isometry);
static void			make_isometry_array(IsometryList *isometry_list, Isometry *the_linked_list);


/* This computes the symmetry group generators, but still have to
   decide what to do about getting a canonical triangulation. */ 

Boolean get_symmetry_transformations(Boolean canonize_flag, Boolean orientable_only_flag, Triangulation* theTriangulation, int *num_symmetries, MoebiusTransformation **symmetry_list)
{

  int theNumSymmetryGenerators, theNumDesiredSymmetryGenerators;
  MoebiusTransformation *theSymmetryMoebiusGenerators;
  SymmetryList *theSymmetryList;
  Complex (*v)[4];
  Tetrahedron *tet;
  Complex a[4], b[4], cr_a, cr_b, kk, b1k, normalization;
  int i, j, k;


  if (get_filled_solution_type(theTriangulation) != geometric_solution
      && get_filled_solution_type(theTriangulation) != nongeometric_solution)
    {
      printf("manifold isn't hyperbolic\n");
      return FALSE;
    }

  if (all_cusps_are_complete(theTriangulation) == FALSE)
    {
      printf("no Dehn fillings are allowed\n");
      return FALSE;
    }

  if (orientable_only_flag == TRUE
      && get_orientability(theTriangulation) == nonorientable_manifold)
    {
      printf("It makes no sense to restrict to orientation preserving symmetries of a nonorientable manifold.\n");
      return FALSE;
    }

  if (canonize_flag == TRUE)
    {
      proto_canonize(theTriangulation);

      if (is_canonical_triangulation(theTriangulation) == FALSE)
	{
	  printf("*** Canonical cell decomposition is not a triangulation. ***\n");
	  return FALSE;
	}
    }

  choose_generators(theTriangulation, TRUE, FALSE);
  /* counts the generators so we can allocate the arrays */

  //	To compute the symmetries we'll use our own personalized
  //	version of the SnapPea kernel's compute_cusped_isometries()
  //	function.  We want to insure that we use theTriangulation
  //	itself, not a canonized copy.
  compute_symmetries(theTriangulation, &theSymmetryList);
  theNumSymmetryGenerators = theSymmetryList->num_isometries;
  theSymmetryMoebiusGenerators = NEW_ARRAY(theNumSymmetryGenerators, MoebiusTransformation);

  typedef Complex C4[4];

  v = NEW_ARRAY(theTriangulation->num_tetrahedra, C4); 

  for (tet = theTriangulation->tet_list_begin.next, i = 0;
       tet != &theTriangulation->tet_list_end;
       tet = tet->next, i++)
    for (j = 0; j < 4; j++)
      v[i][j] = tet->corner[j];
  
  for (i = 0; i < theNumSymmetryGenerators; i++)
    {
      for (j = 0; j < 4; j++)
	{
	  // base tetrahedron corners
	  a[j] = v[0][j];

	  // image tetrahedron corners
	  b[j] = v[theSymmetryList->isometry[i]->tet_image[0]]
	    [EVALUATE(theSymmetryList->isometry[i]->tet_map[0], j)];
	}

      //	compute crossratios to decide whether the isometry
      //	is orientation preserving
      cr_a = compute_cross_ratio(a);
      cr_b = compute_cross_ratio(b);
      theSymmetryMoebiusGenerators[i].parity =
	((cr_a.imag > 0.0) == (cr_b.imag > 0.0)) ?
	orientation_preserving :
	orientation_reversing;

      //	If the MoebiusTransformation is orientation_reversing,
      //	we want to compute a function of z-bar, as explained
      //	in the documentation accompanying the definition of a
      //	MoebiusTransformation in SnapPea.h.
      if (theSymmetryMoebiusGenerators[i].parity == orientation_reversing)
	for (j = 0; j < 4; j++)
	  a[j] = complex_conjugate(a[j]);
      
      //	The formula for the Moebius transformation is copied
      //	from compute_one_generator() in matrix_generators.c
      kk = complex_div(
		       complex_mult(complex_minus(b[2],b[0]), complex_minus(a[2],a[1])),
		       complex_mult(complex_minus(b[2],b[1]), complex_minus(a[2],a[0]))
		       );
      b1k = complex_mult(b[1], kk);
      normalization = complex_sqrt(
				   complex_div(
					       One, 
					       complex_mult(kk,
							    complex_mult(
									 complex_minus(a[1],a[0]),
									 complex_minus(b[1],b[0])
									 )
							    )
					       )
				   );
      
      theSymmetryMoebiusGenerators[i].matrix[0][0] = complex_mult(
								  normalization,
								  complex_minus(b1k, b[0])
								  );
      theSymmetryMoebiusGenerators[i].matrix[0][1] = complex_mult(
								  normalization,
								  complex_minus(
										complex_mult(b[0], a[1]),
										complex_mult(b1k,  a[0])
										)
								  );
      theSymmetryMoebiusGenerators[i].matrix[1][0] = complex_mult(
								  normalization,
								  complex_minus(kk, One)
								  );
      theSymmetryMoebiusGenerators[i].matrix[1][1] = complex_mult(
								  normalization,
								  complex_minus(
										a[1],
										complex_mult(kk,a[0])
										)
								  );
    }

  my_free_array(v);
  free_isometry_list(theSymmetryList);

  if (orientable_only_flag == TRUE)
    {
      for (i = 0, theNumDesiredSymmetryGenerators = 0; i < theNumSymmetryGenerators; i++)
	if (theSymmetryMoebiusGenerators[i].parity == orientation_preserving)
	  {
	    theSymmetryMoebiusGenerators[theNumDesiredSymmetryGenerators] = theSymmetryMoebiusGenerators[i];
	    theNumDesiredSymmetryGenerators++;
	  }
    }
  else
    theNumDesiredSymmetryGenerators = theNumSymmetryGenerators;

  *num_symmetries = theNumDesiredSymmetryGenerators;
  *symmetry_list = theSymmetryMoebiusGenerators;
  
  return TRUE; 
}


static void compute_symmetries(
	Triangulation	*aTriangulation,
	IsometryList	**aSymmetryList)
{
	//	For documentation, please see compute_cusped_isometries()
	//	in isometry_cusped.c.

	Isometry		*partial_isometry_list,
					*new_isometry;
	Tetrahedron		*tet0,
					*tet1;
	int				i;

	*aSymmetryList = NEW_STRUCT(IsometryList);
	(*aSymmetryList)->num_isometries	= 0;
	(*aSymmetryList)->isometry			= NULL;

	partial_isometry_list = NULL;

	number_the_tetrahedra(aTriangulation);

	tet0 = aTriangulation->tet_list_begin.next;

	for (tet1 = aTriangulation->tet_list_begin.next;
		 tet1 != &aTriangulation->tet_list_end;
		 tet1 = tet1->next)

		for (i = 0; i < 24; i++)

			if (attempt_isometry(aTriangulation, tet0, tet1, permutation_by_index[i]) == func_OK)
			{
				copy_isometry(aTriangulation, aTriangulation, &new_isometry);

				new_isometry->next = partial_isometry_list;
				partial_isometry_list = new_isometry;

				(*aSymmetryList)->num_isometries++;
			}

	make_isometry_array(*aSymmetryList, partial_isometry_list);
}


/*
 *	attempt_isometry() checks whether sending tet0 (of manifold0)
 *	to tet1 (of manifold1, which may or may equal manifold0)
 *	defines an isometry between the two manifolds.  If it does,
 *	it leaves the isometry in the "image" and "map" fields of
 *	the Tetrahedra of manifold0 (see the documentation at the top
 *	of isometry.h for details), and returns func_OK.  If it doesn't
 *	define an isometry, it leaves garbage in the "image" and "map"
 *	fields and returns func_failed.
 *
 *	Technical note:  attempt_isometry() could be written using a simple
 *	recursive function to set the image and map fields, but I'm trying
 *	to avoid recursions, because on some machines (e.g. older Macs)
 *	a stack-heap collision is a possibility if one has a deep recursion
 *	with a lot of nonstatic local variables, whereas in the future static
 *	local variables are likely to cause problems in a multithreaded
 *	environment.  So . . . better to avoid the recursion.  The actual
 *	algorithm is documented in the following code.
 */

static FuncResult attempt_isometry(
	Triangulation	*manifold0,
	Tetrahedron		*initial_tet0,
	Tetrahedron		*initial_tet1,
	Permutation		initial_map)
{
	Tetrahedron	*tet0,
				*tet1,
				*nbr0,
				*nbr1,
				**queue;
	int			first,
				last;
	FaceIndex	face0,
				face1;
	Permutation	gluing0,
				gluing1,
				nbr0_map;

	/*
	 *	initial_tet1 and initial_map are arbitrary, so
	 *	the vast majority of calls to attempt_isometry()
	 *	will fail.  Therefore it's worth include a quick
	 *	plausibility check at the beginning, to make the
	 *	algorithm run faster.
	 */
	if (is_isometry_plausible(initial_tet0, initial_tet1, initial_map) == FALSE)
		return func_failed;

	/*
	 *	Initialize all the image fields of manifold0
	 *	to NULL to show they haven't been set.
	 */
	for (	tet0 = manifold0->tet_list_begin.next;
		 	tet0 != &manifold0->tet_list_end;
			tet0 = tet0->next)
		tet0->image = NULL;

	/*
	 *	Allocate space for a queue which is large enough
	 *	to hold pointers to all the Tetrahedra in manifold0.
	 */
	queue = NEW_ARRAY(manifold0->num_tetrahedra, Tetrahedron *);

	/*
	 *	At all times, the Tetrahedra on the queue will be those which
	 *
	 *		(1) have set their image and map fields, but
	 *
	 *		(2)	have not checked their neighbors.
	 */

	/*
	 *	Set the image and map fields for initial_tet0.
	 */
	initial_tet0->image	= initial_tet1;
	initial_tet0->map	= initial_map;

	/*
	 *	Put initial_tet0 on the queue.
	 */
	first = 0;
	last  = 0;
	queue[first] = initial_tet0;

	/*
	 *	While there are Tetrahedra on the queue . . .
	 */
	while (last >= first)
	{
		/*
		 *	Pull the first Tetrahedron off the queue and call it tet0.
		 */
		tet0 = queue[first++];

		/*
		 *	tet0 maps to some Tetrahedron tet1 in manifold1.
		 */
		tet1 = tet0->image;

		/*
		 *	For each face of tet0 . . .
		 */
		for (face0 = 0; face0 < 4; face0++)
		{
			/*
			 *	Let nbr0 be the Tetrahedron which meets tet0 at face0.
			 */
			nbr0 = tet0->neighbor[face0];

			/*
			 *	tet0->map takes face0 of tet0 to face1 of tet1.
			 */
			face1 = EVALUATE(tet0->map, face0);

			/*
			 *	Let nbr1 be the Tetrahedron which meets tet1 at face1.
			 */
			nbr1 = tet1->neighbor[face1];

			/*
			 *	Let gluing0 be the gluing which identifies face0 of
			 *	tet0 to nbr0, and similarly for gluing1.
			 */
			gluing0 = tet0->gluing[face0];
			gluing1 = tet1->gluing[face1];

			/*
			 *						 gluing0
			 *				   tet0  ------>  nbr0
			 *					|				|
			 *		  tet0->map |				| nbr0->map
			 *					|				|
			 *					V    gluing1	V
			 *				   tet1  ------>  nbr1
			 *
			 *	We want to ensure that tet1 and nbr1 enjoy the same
			 *	relationship to each other in manifold1 that tet0 and
			 *	nbr0 do in manifold0.  The conditions
			 *
			 *					nbr0->image == nbr1
			 *	and
			 *		nbr0->map == gluing1 o tet0->map o gluing0^-1
			 *
			 *	are necessary and sufficient to insure that we have a
			 *	combinatorial equivalence between the Triangulations.
			 *	(The proof relies on the fact that we've already checked
			 *	(near the beginning of compute_cusped_isometries() above)
			 *	that the Triangulations have the same number of Tetrahedra;
			 *	otherwise one Triangulation could be a (possibly branched)
			 *	covering of the other.)
			 */

			/*
			 *	Compute the required value for nbr0->map.
			 */
			nbr0_map = compose_permutations(
							compose_permutations(
								gluing1,
								tet0->map
							),
							inverse_permutation[gluing0]
						);

			/*
			 *	If nbr0->image and nbr0->map have already been set,
			 *	check that they satisfy the above conditions.
			 */
			if (nbr0->image != NULL)
			{
				if (nbr0->image != nbr1  ||  nbr0->map != nbr0_map)
				{
					/*
					 *	This isn't an isometry.
					 */
					my_free(queue);
					return func_failed;
				}
			}
			/*
			 *	else . . .
			 *	nbr0->image and nbr0->map haven't been set.
			 *	Set them, and put nbr0 on the queue.
			 */
			else
			{
				nbr0->image	= nbr1;
				nbr0->map	= nbr0_map;
				queue[++last] = nbr0;
			}
		}
	}

	/*
	 *	A quick error check.
	 *	Is it plausible that each Tetrahedron
	 *	has been on the queue exactly once?
	 */
	if (first != manifold0->num_tetrahedra
	 || last  != manifold0->num_tetrahedra - 1)
		uFatalError("attempt_isometry", "isometry");

	/*
	 *	Free the queue, and report success.
	 */
	my_free(queue);
	return func_OK;
}


static Boolean is_isometry_plausible(
	Tetrahedron		*initial_tet0,
	Tetrahedron		*initial_tet1,
	Permutation		initial_map)
{
	/*
	 *	To check whether an Isometry taking initial_tet0 to
	 *	initial_tet1 via initial_map is even plausible, let's
	 *	check whether their EdgeClass orders match up.
	 */

	int	i,
		j;

	for (i = 0; i < 4; i++)

		for (j = i + 1; j < 4; j++)

			if (initial_tet0->edge_class[edge_between_vertices[i][j]]->order
			 != initial_tet1->edge_class[edge_between_vertices[EVALUATE(initial_map, i)][EVALUATE(initial_map, j)]]->order)

				return FALSE;

	return TRUE;
}


/*
 *	copy_isometry() assumes the "image" and "map" fields of
 *	manifold0 describe an isometry to a (not necessarily
 *	distinct) manifold1.  It also assumes that the Tetrahedra
 *	in both manifolds have been numbered as described in
 *	symmetry.h.  It allocates an Isometry data structure,
 *	writes the isometry into it, and sets *new_isometry to
 *	point to it.
 */

static void copy_isometry(
	Triangulation	*manifold0,
	Triangulation	*manifold1,
	Isometry		**new_isometry)
{
	Tetrahedron	*tet0;
	int			i;

	/*
	 *	Allocate the Isometry.
	 */

	*new_isometry				= NEW_STRUCT(Isometry);
	(*new_isometry)->tet_image	= NEW_ARRAY(manifold0->num_tetrahedra, int);
	(*new_isometry)->tet_map	= NEW_ARRAY(manifold0->num_tetrahedra, Permutation);
	(*new_isometry)->cusp_image	= NEW_ARRAY(manifold0->num_cusps, int);
	(*new_isometry)->cusp_map	= NEW_ARRAY(manifold0->num_cusps, MatrixInt22);

	/*
	 *	Set the num_tetrahedra and num_cusps fields.
	 */

	(*new_isometry)->num_tetrahedra	= manifold0->num_tetrahedra;
	(*new_isometry)->num_cusps		= manifold0->num_cusps;

	/*
	 *	Copy the isometry from the Triangulation
	 *	to the Isometry data structure.
	 */

	for (tet0 = manifold0->tet_list_begin.next, i = 0;
		 tet0 != &manifold0->tet_list_end;
		 tet0 = tet0->next, i++)
	{
		(*new_isometry)->tet_image[i]	= tet0->image->index;
		(*new_isometry)->tet_map[i]		= tet0->map;
	}

	/*
	 *	How does this Isometry act on the Cusps?
	 */

	compute_cusp_action(manifold0, manifold1, *new_isometry);

	/*
	 *	We don't use the "next" field.
	 */

	(*new_isometry)->next = NULL;
}


/*
 *	Given a Triangulation *manifold0 whose "image" and "map" fields
 *	describe an isometry to some (not necessarily distinct) manifold,
 *	compute_cusp_action() computes the action on the Cusps.  Only real
 *	cusps are considered -- finite vertices are ignored.
 */

static void compute_cusp_action(
	Triangulation	*manifold0,
	Triangulation	*manifold1,
	Isometry		*isometry)
{
	compute_cusp_image(manifold0, isometry);
	compute_cusp_map(manifold0, manifold1, isometry);
	isometry->extends_to_link = does_isometry_extend_to_link(isometry);
}


static void compute_cusp_image(
	Triangulation	*manifold0,
	Isometry		*isometry)
{
	Tetrahedron	*tet;
	VertexIndex	v;

	/*
	 *	Examine each vertex of each Tetrahedron to see which
	 *	cusp is being taken to which.
	 *
	 *	There's a tremendous amount of redundancy here -- each
	 *	cusp image is set over and over -- but it hardly matters.
	 *
	 *	Ignore negatively indexed Cusps -- they're finite vertices.
	 */

	for (tet = manifold0->tet_list_begin.next;
		 tet != &manifold0->tet_list_end;
		 tet = tet->next)

		for (v = 0; v < 4; v++)

			if (tet->cusp[v]->index >= 0)

				isometry->cusp_image[tet->cusp[v]->index]
					= tet->image->cusp[EVALUATE(tet->map, v)]->index;
}


static void compute_cusp_map(
	Triangulation	*manifold0,
	Triangulation	*manifold1,
	Isometry		*isometry)
{
	Tetrahedron	*tet;
	VertexIndex	v;
	int			i;

	/*
	 *	Copy the manifold1's peripheral curves into
	 *	scratch_curves[0], and copy the images of manifold0's
	 *	peripheral curves into scratch_curves[1].
	 *
	 *	When the manifold is orientable and the Isometry is
	 *	orientation-reversing, and sometimes when the manifold
	 *	is nonorientable, the images of the peripheral curves
	 *	of a torus will lie on the "wrong" sheet of the Cusp's
	 *	orientation double cover.  (See peripheral_curves.c for
	 *	background material.)  Therefore we copy the images of
	 *	the peripheral curves of torus cusps to both sheets of
	 *	the orientation double cover, to guarantee that the
	 *	intersection numbers come out right.
	 */

	copy_curves_to_scratch(manifold1, 0, FALSE);
	copy_images_to_scratch(manifold0, 1, TRUE);

	/*
	 *	Compute the intersection numbers of the images of manifold0's
	 *	peripheral curves with manifold1's peripheral curves..
	 */

	compute_intersection_numbers(manifold1);

	/*
	 *	Now extract the cusp_maps from the linking numbers.
	 *
	 *	There's a lot of redundancy in this loop, but a trivial
	 *	computation so the redundancy hardly matters.
	 *
	 *	Ignore negatively indexed Cusps -- they're finite vertices.
	 */

	for (tet = manifold0->tet_list_begin.next;
		 tet != &manifold0->tet_list_end;
		 tet = tet->next)

		for (v = 0; v < 4; v++)

			if (tet->cusp[v]->index >= 0)

				for (i = 0; i < 2; i++)		/* i = M, L */
				{
					isometry->cusp_map[tet->cusp[v]->index][M][i]
						= + tet->image->cusp[EVALUATE(tet->map, v)]->intersection_number[L][i];

					isometry->cusp_map[tet->cusp[v]->index][L][i]
						= - tet->image->cusp[EVALUATE(tet->map, v)]->intersection_number[M][i];
				}
}


static void copy_images_to_scratch(
	Triangulation	*manifold0,
	int				which_set,
	Boolean			double_copy_on_tori)
{
	Tetrahedron	*tet;
	int			i,
				j,
				jj,
				k,
				kk,
				l,
				ll;

	/*
	 *	This function is modelled on copy_curves_to_scratch()
	 *	in intersection_numbers.c.
	 */

	/*
	 *	Note that even though we are passed manifold0 as an
	 *	explicit function parameter, we are ultimately writing
	 *	the image curves in manifold1.
	 */

	for (tet = manifold0->tet_list_begin.next;
		 tet != &manifold0->tet_list_end;
		 tet = tet->next)

		for (i = 0; i < 2; i++)

			for (k = 0; k < 4; k++)
			{
				kk = EVALUATE(tet->map, k);

				for (l = 0; l < 4; l++)
				{
					ll = EVALUATE(tet->map, l);

					if (tet->cusp[k]->topology == torus_cusp
					 && double_copy_on_tori == TRUE)

						tet->image->scratch_curve[which_set][i][right_handed][kk][ll] =
						tet->image->scratch_curve[which_set][i][ left_handed][kk][ll] =
							  tet->curve[i][right_handed][k][l]
							+ tet->curve[i][ left_handed][k][l];

					else
						/*
						 *		tet->cusp[k]->topology == Klein_cusp
						 *	 || double_copy_on_tori == FALSE
						 */

						for (j = 0; j < 2; j++)
						{
							/*
							 *	parities can be tricky.
							 *
							 *	When discussing gluings from a face of one Tetrahedron
							 *	to a face of another, an even parity corresponds to an
							 *	orientation-reversing gluing.
							 *
							 *	When discussing a map from one Tetrahedron onto
							 *	another, an even parity is an orientation-preserving
							 *	map.
							 */

							/*
							 *	If tet->map has even parity (i.e. if it's an
							 *	orientation-preserving map) it will send the
							 *	right-handed vertex cross section to a right-
							 *	handed image.
							 *
							 *	If tet->map has odd parity (i.e. if it's an
							 *	orientation-reversing map) it will send the
							 *	right-handed vertex cross section to a left-
							 *	handed image.
							 */
							jj = (parity[tet->map] == 0) ? j : !j;

							tet->image->scratch_curve[which_set][i][jj][kk][ll]
								   = tet->curve[i][j][k][l];
						}
				}
			}
}


static Boolean does_isometry_extend_to_link(
	Isometry	*isometry)
{
	/*
	 *	This function assumes the cusp_maps have been
	 *	computed, and checks whether they all take
	 *	meridians to meridians.
	 */

	int	i;

	for (i = 0; i < isometry->num_cusps; i++)

		if (isometry->cusp_map[i][L][M] != 0)

			return FALSE;

	return TRUE;
}


static void make_isometry_array(
	IsometryList	*isometry_list,
	Isometry		*the_linked_list)
{
	int			i;
	Isometry	*an_isometry;

	if (isometry_list->num_isometries == 0)

		isometry_list->isometry = NULL;

	else
	{
		isometry_list->isometry = NEW_ARRAY(isometry_list->num_isometries, Isometry *);

		for (	an_isometry = the_linked_list, i = 0;
				an_isometry != NULL && i < isometry_list->num_isometries;
				an_isometry = an_isometry->next, i++)

			isometry_list->isometry[i] = an_isometry;

		/*
		 *	A quick error check.
		 */
		if (an_isometry != NULL || i != isometry_list->num_isometries)
			uFatalError("make_isometry_array", "isometry");
	}
}


static Complex compute_cross_ratio(
	Complex	a[4])
{
	return complex_div
		(
			complex_mult
			(
				complex_minus(a[3], a[1]),
				complex_minus(a[2], a[0])
			),
			complex_mult
			(
				complex_minus(a[2], a[1]),
				complex_minus(a[3], a[0])
			)
		);
}
