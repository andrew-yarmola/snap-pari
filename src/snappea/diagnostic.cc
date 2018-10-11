/* 
 *  Added by oag. functions which print stuff on standard out. 
 *  Called by snap/tube and useful for debugging. 
 */

#include "kernel.h"

using std::cout; 

void print_holonomies(Triangulation* manifold)
{
  int n = get_num_cusps(manifold);
  int cusp;
  Complex mh, lh; 

  // Print out the holonomies and corresponding eigenvalues. 
  cout << "H=[";
  for (cusp=0; cusp < n; ++cusp) {
    if (cusp > 0) cout << ", ";
    get_holonomy(manifold, cusp, &mh, &lh, 0, 0);
    cout << "[" << mh << ", " << lh << "]";
  }
  cout << "]\n"; 
  cout << "E=[";
  for (cusp=0; cusp < n; ++cusp) {
    if (cusp > 0) cout << ", ";
    get_holonomy(manifold, cusp, &mh, &lh, 0, 0);
    cout << "[" << complex_exp(mh/2.0) << ", " << complex_exp(lh/2.0) << "]";
  }
  cout << "]\n"; 
}

void print_edge_angle_sums(Triangulation* manifold)
{
    EdgeClass   *edge;

    cout << "[";
    for (edge = manifold->edge_list_begin.next;
	 edge != &manifold->edge_list_end;
	 edge = edge->next) {
      
      if (edge != manifold->edge_list_begin.next) cout << ", ";
      cout << edge->edge_angle_sum; 
    }
    cout << "]\n";
}

