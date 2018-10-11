/*
 *  Dirichlet_basepoint2D.c
 *
 */

#include "Dirichlet2D_prototypes.h"
#include <stdlib.h>     /* needed for qsort() */


#define BASEPOINT_EPSILON       (1e4 * DBL_EPSILON)


#define BIG_BASEPOINT_EPSILON   1e-5

#define MAX_TOTAL_DISTANCE      0.25

#define DERIVATIVE_EPSILON      1e-6

#define MAX_STEP_SIZE           1.0

#define IDENTITY_EPSILON        1e-6

#define CONSTRAINT_EPSILON      1e-6


#define SADDLE_EPSILON          1e-6


#define ZERO_DERIV_EPSILON      1e-6

#define MIN_PIVOT               (1e5 * DBL_EPSILON)


#define EVALUATE_EQN(eqn, pt)   \
    (eqn[0]*pt[0] + eqn[1]*pt[1] + eqn[2]*pt[2] + eqn[3])


typedef double ObjectiveFunction[4];

typedef double Constraint[4];

typedef double Solution[3];

static int          count_matrix_pairs(AlgMatrixPairList *gen_list);
static void         verify_gen_list(AlgMatrixPairList *gen_list, int num_matrix_pairs);
static FuncResult   set_objective_function(ObjectiveFunction objective_function, AlgMatrixPair *matrix_pair);
static void         step_size_constraints(Constraint *constraints, ObjectiveFunction objective_function);
static void         regular_constraints(Constraint *constraints, AlgMatrixPairList *gen_list, ObjectiveFunction objective_function, Boolean *may_be_saddle_point);
static void         linear_programming(ObjectiveFunction objective_function, int num_constraints, Constraint *constraints, Solution solution);
static Boolean      apex_is_higher(double height1, double height2, Solution apex1, Solution apex2);
static FuncResult   solve_two_equations(Constraint *equations[3], Solution solution);
static void         initialize_t2(Solution solution, O31_matrix& t2);
static void         sort_gen_list(AlgMatrixPairList *gen_list, int num_matrix_pairs);
static int CDECL    compare_image_height(const void *ptr1, const void *ptr2);
static double       length3(double v[3]);
static double       inner3(double u[3], double v[3]);
static void         copy3(Solution dest, const Solution source);


void maximize_injectivity_radius(
    AlgMatrixPairList   *gen_list,
    Boolean         *basepoint_moved,
	Boolean         interactive)
{
  int           num_matrix_pairs;
  double        distance_moved,
	    prev_distance_moved,
	    total_distance_moved;
  Boolean           keep_going;
  ObjectiveFunction objective_function;
  int           num_constraints;
  Constraint        *constraints;
  Solution      solution;
  Boolean       may_be_saddle_point,
			 saddle_query_given;
  int           choice;
  
  static const Solution zero_solution = {ZERO, ZERO, ZERO},
  small_displacement = {0.001734, 0.002035, 0.000000};
  
  static char       *saddle_message = "The basepoint may be at a saddle point of the injectivity radius function.";
  static char       num_saddle_responses = 2;
  static char       *saddle_responses[2] = {
						 "Continue On",
			 "Stop Here and See Dirichlet Domain"};
  static int        saddle_default = 1;
  
  static char       *zero_deriv_message = "The derivative of the distance to the closest translate of the basepoint is zero.";
  static char       num_zero_deriv_responses = 2;
  static char       *zero_deriv_responses[2] = {
						     "Displace Basepoint and Continue On",
			     "Stop Here and See Dirichlet Domain"};
  static int        zero_deriv_default = 1;
  
  
  num_matrix_pairs = count_matrix_pairs(gen_list);

  verify_gen_list(gen_list, num_matrix_pairs);

  *basepoint_moved = FALSE;

  total_distance_moved = ZERO;

  prev_distance_moved = DBL_MAX;

  saddle_query_given = FALSE;


  do
    {

      
      if (set_objective_function(objective_function, gen_list->begin.next->next) == func_OK)
    {
      num_constraints = (num_matrix_pairs - 2) + 2;
      constraints = NEW_ARRAY(num_constraints, Constraint);
      
      step_size_constraints(constraints, objective_function);
	       
      regular_constraints(constraints, gen_list, objective_function, &may_be_saddle_point);

      if (may_be_saddle_point == FALSE)
	linear_programming(objective_function, num_constraints, constraints, solution);
	    
      else
	{
	  if (interactive == TRUE)
	{
	  if (saddle_query_given == FALSE)
	    {
	      choice = uQuery(  saddle_message,
		      num_saddle_responses,
		      saddle_responses,
		      saddle_default);
	      saddle_query_given = TRUE;
	    }
	  else
	    choice = 0; /* continue on */
	}
	  else
	
	choice = 0; /* continue on */
	  
	  switch (choice)
	{
	case 0:
	  
	  linear_programming(objective_function, num_constraints, constraints, solution);
	  break;
	  
	case 1:
	  copy3(solution, zero_solution);
	  *basepoint_moved = FALSE;
	  break;
	}
	}
      
      my_free_array(constraints);
    }
      
	else
	  {
	    
	    switch (    interactive == TRUE ?
		
		uQuery( zero_deriv_message,
		       num_zero_deriv_responses,
		       zero_deriv_responses,
		       zero_deriv_default) :

		
			0   /* continue on */ )
	      {
	      case 0:
		copy3(solution, small_displacement);
		break;
		
	      case 1:
		copy3(solution, zero_solution);
		*basepoint_moved = FALSE;
		break;
	      }
	      }

      conjugate_matrices(gen_list, solution);

      sort_gen_list(gen_list, num_matrix_pairs);
      
      distance_moved = length3(solution);
      
      total_distance_moved += distance_moved;


      if (distance_moved > BASEPOINT_EPSILON)
    {
      *basepoint_moved = TRUE;
      keep_going = (total_distance_moved < MAX_TOTAL_DISTANCE);
    }
      else
	    keep_going = FALSE;
      
      
      if (prev_distance_moved < BIG_BASEPOINT_EPSILON
	 &&      distance_moved < BIG_BASEPOINT_EPSILON)
    {

      keep_going = FALSE;

      if (total_distance_moved < BIG_BASEPOINT_EPSILON)
	*basepoint_moved = FALSE;
    }
      prev_distance_moved = distance_moved;

    } while (keep_going == TRUE);
}


static int count_matrix_pairs(
    AlgMatrixPairList   *gen_list)
{
    int         num_matrix_pairs;
    AlgMatrixPair   *matrix_pair;

    num_matrix_pairs = 0;

    for (   matrix_pair = gen_list->begin.next;
	    matrix_pair != &gen_list->end;
	    matrix_pair = matrix_pair->next)

	num_matrix_pairs++;

    return num_matrix_pairs;
}



static void verify_gen_list(
    AlgMatrixPairList   *gen_list,
    int             num_matrix_pairs)
{
    AlgMatrixPair   *matrix_pair;


    if (num_matrix_pairs < 2)
	uFatalError("verify_gen_list", "Dirichlet_basepoint");



    if (gen_list->begin.next->height > ONE + IDENTITY_EPSILON)
	uFatalError("verify_gen_list", "Dirichlet_basepoint");



    for (   matrix_pair = gen_list->begin.next->next;
	    matrix_pair != &gen_list->end;
	    matrix_pair = matrix_pair->next)

	if (matrix_pair->height < matrix_pair->prev->height)

	    uFatalError("verify_gen_list", "Dirichlet_basepoint");
}


static FuncResult set_objective_function(
    ObjectiveFunction   objective_function,
    AlgMatrixPair           *matrix_pair)
{
  int       i;

  for (i = 0; i < 2; i++)
    objective_function[i] = matrix_pair->m[0]->matrix(0,i+1) - matrix_pair->m[0]->matrix(i+1,0);
  objective_function[2] = 0;
  objective_function[3] = matrix_pair->m[0]->matrix(0,0);


  if (length3(objective_function) > DERIVATIVE_EPSILON)
    return func_OK;
  else
    return func_failed;
}


static void step_size_constraints(
    Constraint          *constraints,
    ObjectiveFunction   objective_function)
{
  int       i,
		j,
		i0;
  double    v[2][3],
		 w[2][3],
		 max_abs,
	    length;

  length = length3(objective_function);
  for (i = 0; i < 3; i++)
    v[0][i] = objective_function[i] / length;

/*  max_abs = ZERO;
    for (i = 0; i < 3; i++)
	if (fabs(v[0][i]) > max_abs)
	{
	    i0 = i;
	    max_abs = fabs(v[0][i]);
	}
*/
    
    v[1][0]     = -v[0][1] ;
    v[1][1] =  v[0][0];
    v[1][2] = ZERO;

/*  length = length3(v[1]);
    for (i = 0; i < 3; i++)
	v[1][i] /= length;
*/
    
/*  for (i = 0; i < 3; i++)
	v[2][i] = v[0][(i+1)%3] * v[1][(i+2)%3]  -  v[0][(i+2)%3] * v[1][(i+1)%3];*/

    
    for (j = 0; j < 3; j++)
    {
	w[0][j] = v[0][j] + (HALF*v[1][j]);
	w[1][j] = v[0][j] + (-HALF*v[1][j]);
/*      w[2][j] = v[0][j] + (-HALF*v[1][j] - ROOT_3_OVER_2*v[2][j]);*/
    }


    for (i = 0; i < 2; i++)
    {
	for (j = 0; j < 3; j++)
	    constraints[i][j] = w[i][j];
	constraints[i][3] = -MAX_STEP_SIZE;
    }
}


static void regular_constraints(
    Constraint          *constraints,
    AlgMatrixPairList       *gen_list,
    ObjectiveFunction   objective_function,
    Boolean             *may_be_saddle_point)
{

    int         i;
    AlgMatrixPair   *matrix_pair;
    Constraint  *constraint;
    double      h[4],
		c;

    
    *may_be_saddle_point = FALSE;

    

    for (   matrix_pair = gen_list->begin.next->next->next,
		constraint = constraints + 2;

	    matrix_pair != &gen_list->end;

	    matrix_pair = matrix_pair->next,
		constraint++)
    {
    

      for (i = 0; i < 2; i++)
	h[i] = matrix_pair->m[0]->matrix(0,i+1) - matrix_pair->m[0]->matrix(i+1,0);
      h[2] =  0;
      h[3] = matrix_pair->m[0]->matrix(0,0);

    
      for (i = 0; i < 4; i++)
	(*constraint)[i] = objective_function[i] - h[i];

    

      if ((*constraint)[3] > - CONSTRAINT_EPSILON)
	{
	  
	  
	  if (length3(*constraint) > ZERO_DERIV_EPSILON)
	{
	  
	  c = inner3(objective_function, *constraint) /
	    (length3(objective_function) * length3(*constraint));

	    
	  
	  if (fabs(c) > ONE - SADDLE_EPSILON)
	    *may_be_saddle_point = TRUE;

	}
	}
    }
      }


static void linear_programming(
    ObjectiveFunction   objective_function,
    int         num_constraints,
    Constraint      *constraints,
    Solution        solution)
{
    int         i,
		j,
		k;
    Constraint  *active_constraints[3],
		*new_constraints[3];
    Solution    apex,
		new_apex,
		max_apex;
    double      apex_height,
		new_height,
		max_height;
    int         inactive_constraint_index;

    for (i = 0; i < 2; i++)
	active_constraints[i] = constraints + i;

    if (solve_two_equations(active_constraints, apex) == func_failed)
	uFatalError("linear_programming", "Dirichlet_basepoint");

    apex_height = EVALUATE_EQN(objective_function, apex);

    for (i = 0; i < num_constraints; i++)

	if (EVALUATE_EQN(constraints[i], apex) > CONSTRAINT_EPSILON)
	{

	    max_height = -ONE;
	    for (j = 0; j < 3; j++)
		max_apex[j] = ZERO;

	    for (j = 0; j < 2; j++)
	    {

	      for (k = 0; k < 2; k++)
		new_constraints[k] =
		  (k == j ? &constraints[i] : active_constraints[k]);

	      if (solve_two_equations(new_constraints, new_apex) == func_failed)
		continue;

	      new_height = EVALUATE_EQN(objective_function, new_apex);

	      
	      if (apex_is_higher(new_height, apex_height, new_apex, apex) == TRUE)
		continue;

	    
	      if (apex_is_higher(new_height, max_height, new_apex, max_apex) == TRUE)
		{
		  inactive_constraint_index = j;
		  max_height = new_height;
		  for (k = 0; k < 3; k++)
		max_apex[k] = new_apex[k];
		}
	    }

	    
	    active_constraints[inactive_constraint_index] = &constraints[i];

	    
	    for (j = 0; j < 3; j++)
	      apex[j] = max_apex[j];
	    apex_height = max_height;

	    i = -1;
	      }


    for (i = 0; i < 3; i++)
      solution[i] = apex[i];

    return;
    

}


static Boolean apex_is_higher(
    double      height1,
    double      height2,
    Solution    apex1,
    Solution    apex2)
{
    int i;

    if (height1 > height2)
	return TRUE;
    if (height1 < height2)
	return FALSE;

    for (i = 0; i < 2; i++)
    {
	if (apex1[i] > apex2[i])
	    return TRUE;
	if (apex1[i] < apex2[i])
	    return FALSE;
    }

    return FALSE;
}


static FuncResult solve_two_equations(
    Constraint  *equations[2],
    Solution    solution)
{
    
    int     r,
	    c,
	    p;
    double  equation_storage[2][4],
	    *eqn[2],
	    *temp,
	    pivot_value;


    for (r = 0; r < 2; r++)
	eqn[r] = equation_storage[r];

    for (r = 0; r < 2; r++)
	for (c = 0; c < 4; c++)
	    eqn[r][c] = (*equations[r])[c];


    for (p = 0; p < 2; p++)
    {
	
	for (r = p + 1; r < 2; r++)
	    if (fabs(eqn[r][p]) > fabs(eqn[p][p]))
	    {
		temp    = eqn[p];
		eqn[p]  = eqn[r];
		eqn[r]  = temp;
	    }

    
    
	pivot_value = eqn[p][p];

	if (fabs(pivot_value) < MIN_PIVOT)
	    return func_failed;

    
	for (c = p + 1; c < 4; c++)
	    eqn[p][c] /= pivot_value;


	for (r = p + 1; r < 2; r++)
	    for (c = p + 1; c < 4; c++)
		eqn[r][c] -= eqn[r][p] * eqn[p][c];
    }

    
    for (c = 2; --c >= 0; )
	for (r = c; --r >= 0; )
	    eqn[r][3] -= eqn[r][c] * eqn[c][3];

    for (r = 0; r < 2; r++)
	solution[r] = - eqn[r][3];
    solution[2] = 0;

    return func_OK;
}


void conjugate_matrices(
    AlgMatrixPairList       *gen_list,
    double              displacement[3])
{
    

    O31_matrix      t2;
    AlgMatrixPair   *matrix_pair;


    initialize_t2(displacement, t2);

    t2.O31_GramSchmidt();

    for (matrix_pair = gen_list->begin.next;
	 matrix_pair != &gen_list->end;
	 matrix_pair = matrix_pair->next)
    {
      matrix_pair->m[0]->matrix.conjugate_by(t2); 
      matrix_pair->m[1]->matrix = inverse(matrix_pair->m[0]->matrix); 
      matrix_pair->height = matrix_pair->m[0]->matrix(0,0);
    }
}


static void initialize_t2(
    Solution    solution,
    O31_matrix& t2)
{


    int     i,
	    j;


    for (i = 0; i < 3; i++)
	t2(0,i+1) = t2(i+1,0) = solution[i];


    t2(0,0) = ONE;
    for (i = 0; i < 3; i++)
	t2(0,0) += HALF * solution[i] * solution[i];

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    t2(i+1,j+1) = (i == j ? ONE : ZERO)
			 + HALF * solution[i] * solution[j];
}


static void sort_gen_list(
    AlgMatrixPairList   *gen_list,
    int             num_matrix_pairs)
{
    AlgMatrixPair   **array,
		*matrix_pair;
    int         i;

    array = NEW_ARRAY(num_matrix_pairs, AlgMatrixPair *);

    for (matrix_pair = gen_list->begin.next,
	    i = 0;
	 matrix_pair != &gen_list->end;
	 matrix_pair = matrix_pair->next,
	    i++)

	array[i] = matrix_pair;

    if (i != num_matrix_pairs)
	uFatalError("sort_gen_list", "Dirichlet_basepoint");

    qsort(  array,
	    num_matrix_pairs,
	    sizeof(AlgMatrixPair *),
	    compare_image_height);


    gen_list->begin.next = array[0];
    array[0]->prev = &gen_list->begin;
    array[0]->next = array[1];

    for (i = 1; i < num_matrix_pairs - 1; i++)
    {
	array[i]->prev = array[i-1];
	array[i]->next = array[i+1];
    }

    array[num_matrix_pairs - 1]->prev = array[num_matrix_pairs - 2];
    array[num_matrix_pairs - 1]->next = &gen_list->end;
    gen_list->end.prev = array[num_matrix_pairs - 1];

    my_free_array(array);
}


static int CDECL compare_image_height(
    const void  *ptr1,
    const void  *ptr2)
{
    double  diff;

    diff = (*((AlgMatrixPair **)ptr1))->height
	 - (*((AlgMatrixPair **)ptr2))->height;

    if (diff < ZERO)
	return -1;
    if (diff > ZERO)
	return +1;
    return 0;
}


static double length3(
    double  v[3])
{
    double  length;
    int     i;

    length = ZERO;

    for (i = 0; i < 3; i++)
	length += v[i] * v[i];

    length = sqrt(length);  //  no need for safe_sqrt()

    return length;
}


static double inner3(
    double  u[3],
    double  v[3])
{
    double  sum;
    int     i;

    sum = ZERO;

    for (i = 0; i < 3; i++)
	sum += u[i] * v[i];

    return sum;
}


static void copy3(
	    Solution    dest,
    const   Solution    source)
{
    int     i;

    for (i = 0; i < 3; i++)
	dest[i] = source[i];
}







