#ifndef _length_spectrum_
#define _length_spectrum_

#include "SnapPea.h"

/*
 *  The tiling of hyperbolic space by translates gD of the Dirichlet domain
 *  is stored on a binary tree.  Each node in the tree is a Tile structure
 *  containing the group element g associated with the given translate.
 */

struct Tile
{
    /*
     *  A translate gD of the Dirichlet domain is determined
     *  by the group element g.
     */
    O31_matrix       g;

    FGWord            word;
    struct Tile*    parent; 

    /*
     *  Please see complex_length.c for details on how the
     *  complex length is defined and computed.
     */
    Complex         length;

    /*
     *  Is the group element g an orientation_preserving
     *  or orientation_reversing isometry?
     */
    MatrixParity    parity;
    /*
     *  Is the geodesic topologically a circle or a mirrored interval?
     */
    Orbifold1       topology;

    /*
     *  The to_be_eliminated flag is used locally in eliminate_powers()
     *  and eliminate_conjugates().
     */
    Boolean         to_be_eliminated;

    Boolean         peripheral_tile; 

    /*
     *  The tiles are organized in two different way.
     */

    /*
     *  Organization #1.
     *
     *  To make checking for duplicates easy, the Tiles are kept on
     *  a binary tree.  The sort and search key is a more or less
     *  arbitrary function defined in the code.  The next_subtree field
     *  is used locally within already_on_tree() and free_tiling()
     *  to avoid doing recursions on the system stack;  the latter run
     *  the risk of stack/heap collisions.
     */
    struct Tile     *left_child,
		    *right_child;
    double          key;
    struct Tile     *next_subtree;

    /*
     *  Organization #2.
     *
     *  The function tile() needs to keep track of which Tiles have
     *  not yet had their neighbors checked.  Its keeps pending Tiles
     *  on a doubly-linked list.
     */
    struct Tile     *prev,
		    *next;

};

// Usage:
// for (it = root; it; it++) {...}
// or
// TilingIterator it(root);
// while (it) { ... it++; }

class TilingIterator {
  Tile* subtree_stack;

  void init() { if (subtree_stack) subtree_stack->next_subtree = 0; }
public:
  TilingIterator() : subtree_stack(0) {}
  TilingIterator(Tile* root) : subtree_stack(root) { init(); }
  void operator = (Tile* root) { subtree_stack=root; init(); }
  void operator ++ ();
  void operator ++ (int) { ++(*this); }
  Tile& operator * () { return *subtree_stack; }
  Tile* operator ->() { return subtree_stack; }
  Tile* tile() { return subtree_stack; }
  operator bool () { return subtree_stack; }
};


void tile(const WEPolyhedron *polyhedron, double tiling_radius, Tile **tiling);
int count_translates(Tile* root); 
double tiling_radius(Tile* root);
void expand_tiling(const WEPolyhedron *polyhedron, double tiling_radius, Tile *root);
void free_tiling(Tile *root);
void list_geodesics(Tile *tiling, double spine_radius, double cutoff_length, 
		    GeodesicWord **geodesics, int *num_geodesics, Boolean get_words);
double distance_to_origin(O31_matrix const& g, Complex length, bool or_rev);
void eliminate_conjugates(Tile **geodesic_list, int *num_good_geodesics, Tile *tiling, int num_translates, double spine_radius);

#endif
