/*
 *  filling.c
 *
 *  This file contains the functions
 *
 *      Triangulation   *fill_cusps(Triangulation   *manifold,
 *                                  Boolean         fill_cusp[],
 *                                  char            *new_name);
 *
 *      Triangulation   *fill_reasonable_cusps(Triangulation *manifold);
 *
 *      Boolean         cusp_is_fillable(Cusp *cusp);
 *
 *  which the kernel provides to the UI.
 *
 *  fill_cusps() permanently fills k of the cusps of an n-cusp
 *  manifold (k < n).  It returns an ideal Triangulation of the
 *  resulting (n - k)-cusp manifold.
 *
 *  Arguments:
 *
 *      manifold    is the original manifold.  The Dehn filling coefficients
 *                  cusp->m and cusp->l specify how each cusp is to
 *                  be filled.
 *
 *      fill_cusp   says which cusps are to be filled.  The cusp of index i
 *                  will be filled iff fill_cusp[i] is TRUE.
 *
 *      new_name    provides the name for the new Triangulation.
 *
 *  The UI should decide how to present fill_cusps() to the user.
 *  Should all currently Dehn filled cusps be filled at once?
 *  Should the user be presented with a list of check boxes to
 *  specify which cusps to fill?  Should cusps be filled one at a time?
 *  My hope is that fill_cusps() is sufficiently general to support
 *  whatever approach the UI developer prefers.
 *
 *  Having said that, let me now mention fill_reasonable_cusps(), which
 *  makes a decision about which cusps to fill, and then makes a call
 *  to fill_cusp().  fill_reasonable_cusps() will fill all cusps which
 *  have relatively prime Dehn filling coefficients, unless this would
 *  leave no unfilled cusps, in which case it leaves cusp 0 unfilled.
 *  It copies the name from the manifold being filled.
 *
 *  cusp_is_fillable() determines whether an individual cusp is fillable.
 *
 *  The original manifold is always left unaltered.
 *
 *  The files subdivide.c, close_cusps.c, and remove_finite_vertices.c
 *  document the algorithm in detail.
 */

#include "kernel.h"

static void     check_fill_cusp_array(Triangulation *manifold, Boolean fill_cusp[]);
static Boolean  cusp_is_fillable_x(Cusp *cusp);
static Boolean  no_cusps_to_be_filled(int num_cusps, Boolean fill_cusp[]);


Triangulation *fill_cusps(
    Triangulation   *manifold,
    Boolean         fill_cusp[],
    char            *new_name)
{
    Triangulation   *new_triangulation;

    /*
     *  95/10/1  JRW
     *  The following algorithm works correctly even if no cusps are
     *  to be filled, but we can speed it up a bit by simply copying
     *  the Triangulation.
     */
    if (no_cusps_to_be_filled(manifold->num_cusps, fill_cusp) == TRUE)
    {
	copy_triangulation(manifold, &new_triangulation);
	return new_triangulation;
    }

    /*
     *  First let's do a little error checking on the fill_cusp[] array.
     */
    check_fill_cusp_array(manifold, fill_cusp);

    /*
     *  Subdivide the triangulation, introducing finite vertices.
     *  Note that the original triangulation is left unharmed.
     */
    new_triangulation = subdivide(manifold, new_name);

    /*
     *  Close the Cusps specified in the fill_cusp[] array.
     */
    close_cusps(new_triangulation, fill_cusp);

    /*
     *  Retriangulate with no finite vertices.
     */
    remove_finite_vertices(new_triangulation);

    /*
     *  If the old manifold had a hyperbolic structure,
     *  try to find one for the new_triangulation as well.
     */
    if (manifold->solution_type[complete] != not_attempted)
    {
	find_complete_hyperbolic_structure(new_triangulation);
	do_Dehn_filling(new_triangulation);

	/*
	 *  If the old manifold had a known Chern-Simons invariant,
	 *  pass it to the new_triangulation.
	 */
	if (manifold->CS_value_is_known == TRUE)
	{
	    new_triangulation->CS_value_is_known        = manifold->CS_value_is_known;
	    new_triangulation->CS_value[ultimate]       = manifold->CS_value[ultimate];
	    new_triangulation->CS_value[penultimate]    = manifold->CS_value[penultimate];

	    /*
	     *  The solution_type may or may not be good enough to compute
	     *  the fudge factor, but we'll let compute_CS_fudge_from_value()
	     *  worry about that.
	     */
	    compute_CS_fudge_from_value(new_triangulation);
	}
    }

    return new_triangulation;
}


Triangulation *fill_reasonable_cusps(
    Triangulation   *manifold)
{
    Boolean         *fill_cusp;
    Cusp            *cusp;
    int             i;
    Boolean         all_cusps_are_fillable;
    Triangulation   *new_triangulation;

    /*
     *  Allocate the fill_cusp[] array.
     */

    fill_cusp = NEW_ARRAY(manifold->num_cusps, Boolean);

    /*
     *  See which cusps are fillable.
     */

    for (cusp = manifold->cusp_list_begin.next;
	 cusp != &manifold->cusp_list_end;
	 cusp = cusp->next)

	fill_cusp[cusp->index] = cusp_is_fillable_x(cusp);

    /*
     *  If all the cusps are fillable, leave cusp 0 unfilled.
     */

    all_cusps_are_fillable = TRUE;

    for (i = 0; i < manifold->num_cusps; i++)
	if (fill_cusp[i] == FALSE)
	    all_cusps_are_fillable = FALSE;

    if (all_cusps_are_fillable == TRUE)
	fill_cusp[0] = FALSE;

    /*
     *  Call fill_cusps().
     */

    new_triangulation = fill_cusps(manifold, fill_cusp, manifold->name);

    /*
     *  Free the fill_cusp[] array.
     */

    my_free_array(fill_cusp);

    /*
     *  Done.
     */

    return new_triangulation;
}


static void check_fill_cusp_array(
    Triangulation   *manifold,
    Boolean         fill_cusp[])
{
    Boolean at_least_one_cusp_is_left;
    Cusp    *cusp;

    at_least_one_cusp_is_left = FALSE;

    for (cusp = manifold->cusp_list_begin.next;
	 cusp != &manifold->cusp_list_end;
	 cusp = cusp->next)

	if (fill_cusp[cusp->index])
	{
	    if (cusp_is_fillable_x(cusp) == FALSE)
		uFatalError("check_fill_cusp_array", "filling");
	}
	else
	    at_least_one_cusp_is_left = TRUE;

    if (at_least_one_cusp_is_left == FALSE)
	uFatalError("check_fill_cusp_array", "filling");
}


Boolean cusp_is_fillable(               /* For external use */
    Triangulation   *manifold,
    int             cusp_index)
{
    return cusp_is_fillable_x(find_cusp(manifold, cusp_index));
}


static Boolean cusp_is_fillable_x(      /* For internal use */
    Cusp    *cusp)
{
    return( cusp->is_complete == FALSE
	 && Dehn_coefficients_are_relatively_prime_integers(cusp) == TRUE);
}


static Boolean no_cusps_to_be_filled(
    int     num_cusps,
    Boolean fill_cusp[])
{
    int i;

    for (i = 0; i < num_cusps; i++) 

	if (fill_cusp[i] == TRUE)

	    return FALSE;

    return TRUE;
}
