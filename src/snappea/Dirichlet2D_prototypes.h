#ifndef _Dirichlet2D_prototypes_
#define _Dirichlet2D_prototypes_

#include "Dirichlet2D.h"
#include "kernel.h"

#define MATRIX_EPSILON  1e-5

struct AlgMatrixPair
{

    AlgMatrix           *m[2];

    double              height;

/*
 *  The left_ and right_child fields are used locally in
 *  compute_all_products() in Dirichlet_compute.c to build a binary tree
 *  of AlgMatrixPairs.  Normally AlgMatrixPairs are kept on a doubly linked
 *  list, using the prev and next fields.  The next_subtree field is
 *  used even more locally within tree-handling routines, to avoid doing
 *  recursions on the system stack (for fear of stack/ heap collisions).
 */
    struct AlgMatrixPair *left_child,
		*right_child,
		*next_subtree;

    struct AlgMatrixPair *prev,
		*next;
};


struct AlgMatrixPairList
{

    AlgMatrixPair   begin,
	    end;

};


WE2DPolygon *compute_Dirichlet2D_domain(AlgMatrixPairList *gen_list, double vertex_epsilon);

void maximize_injectivity_radius(AlgMatrixPairList *gen_list, Boolean *basepoint_moved, Boolean interactive);
void conjugate_matrices(AlgMatrixPairList *gen_list, double displacement[3]);
void free_matrix_pairs(AlgMatrixPairList *gen_list);
void precise_Alg_product(AlgMatrix *a, AlgMatrix *b, AlgMatrix **product);
void precise_generators(AlgMatrixPairList *gen_list);
void print_word(CyclicWord* word, FILE* out);
void alt_print_word(CyclicWord* word);
void print_vector(O31_vector const& x);

#endif
