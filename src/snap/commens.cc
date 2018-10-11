/*
** Copyright (C) 2003 Oliver A. Goodman <oag@ms.unimelb.edu.au>
**                    and Damian Heard. 
**  
** This file is part of Snap.
** 
** Snap is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software 
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "commens.hh"
#include "snappea/kernel.h"

#include <cstdlib>
#include <cstdio>

using std::cout;
using std::cerr;
using std::endl;

using std::sscanf;
using std::exit;

/*
 *  Max number of tetrahedron in a manifold.
 */
#define MAX_TET 250 

/****************************************************************************/
/*                          GraphNode                                       */
/****************************************************************************/

struct graphnode
{
  int index;
  /*
   *      What is the tetrahedron map?
   */
  Tetrahedron *preimage, *image;

  /*
   *      Under what permutation?
   */
  Permutation permutation;
  
  /*
   *      neighbors.
   */
  GraphNode *nbr[4],
    *next_map,
    *next_queue,
    *next_mapping_graph;
};


/****************************************************************************/
/*                             Tables                                       */
/****************************************************************************/

const int edge_num[] = {0, 0, 1, 2, 1, 0};

/* 40*/
const int arith5[] = {0,1,2,3,4,9,10,25,124,125,126,127,128,129,130,131,
           132,133,134,135,136,139,140,202,203,204,205,206,207,
           208,405,406,407,408,409,410,411,412,413,414};
/* 27*/
const int arith6[] = {118,119,594,595,596,772,773,774,775,776,777,778,779,
           780,781,782,784,786,787,859,955,956,957,958,959,960,961};
/* 8*/
const int arith7[] = {1858, 1859, 2787, 2788, 2789, 2873, 2874, 3551};

#if 0
/* 3*/
const int com_cusp5[]={379,381,382};

/* 4*/
const int com_cusp6[]={843,910,920,930};

/* 51*/
const int com_cusp7[]={2208,2444,2492,2531,2533,2601,2618,2646,2652,2732,2779,
                        2781,2809,2810,2811,2830,2853,2860,2890,2892,2924,2945,2970,2992,
                        3021,3035,3039,3053,3073,3099,3114,3115,3148,3163,3360,3419,3446,
                        3488,3507,3547,3548};
#endif

/*****************************************************************************/
/*                         function prototypes                               */
/*****************************************************************************/

GraphNode        *attempt_tiling_isometry
(Tetrahedron *,Tetrahedron *, Permutation );

GraphNode    *map_list
( int, int, Permutation, GraphNode *map_listing[MAX_TET][MAX_TET]);

Boolean     same_shape 
(Tetrahedron *, Tetrahedron *,Permutation);

Boolean        permutation_type    
(Permutation );

GraphNode *new_node
(Tetrahedron *, Tetrahedron *, Permutation , int);

void print_mapping_graph
(GraphNode *);

void print_permutation
(Permutation );

void print_tet_shapes
(Triangulation *);

Boolean same_mapping_graph
(GraphNode* , GraphNode* );

#if 0
Boolean is_arithmetic
(int , int );
#endif

static void store_transformations
(Triangulation *, GraphNode *, vector<MoebiusTransformation>&);

void name
(int , int , char *);

/*
 *    For more information see README.
 */


/*
 * Both manifolds assumed to be endowed with canonical cell decomposition. 
 * If there are finite vertices, the number of commensurabilities will 
 * be returned but no MoebiusTransformations will be computed. 
 */

int commensurabilities(Triangulation* manifold0, 
		       Triangulation* manifold1, 
		       vector<MoebiusTransformation>& TR, 
		       bool& or_pres, 
		       GraphNode** mapping_graphs)
{
  Tetrahedron *tet0, *tet1;
  int i, num_tiling_isometries=0, num_mapping_graphs=0;

  GraphNode *mapping_graph, *head_mapping_graph=NULL, *mg;

  or_pres = true; 

  number_the_tetrahedra(manifold0);
  if (manifold1!=manifold0) number_the_tetrahedra(manifold1);

  if (get_num_tetrahedra(manifold0) > MAX_TET || 
      get_num_tetrahedra(manifold1) > MAX_TET) {
    cout << "max num tetrahedra exceeded " << get_num_tetrahedra(manifold0) << ' ';
    cout << get_num_tetrahedra(manifold1) << endl; 
    return 0; 
  }

  int nc = get_num_cusps(manifold0); 
  vector<int> fixed(nc);
  vector<int> fix_count(nc,0); 

  /* 
   *    Try mapping an arbitrary but fixed tetrahedron of
   *    manifold0 ("tet0") to each tetrahedron of manifold1
   *    in turn, with each possible permutation.  See which
   *    of these lift to the tilings of H^3 and extended to
   *    isometries.  We're guaranteed to find all possible
   *    isometries this way ( and typically quite a lot of
   *    other garbage as well).
   */
  
  /*
   *    Let tet0 be the first Tetrahedron on manifold0's list.
   */
  tet0=manifold0->tet_list_begin.next;
  
  /*
   * Consider each tetrahedron in manifold1.
   */
  for(tet1=manifold1->tet_list_begin.next;
      tet1!=&manifold1->tet_list_end;
      tet1=tet1->next){
    
    /*
     * Consider all 24 possible mappings from tet0 to tet1.
     */
   
    for(i=0;i<24;i++){
      
      /*
       *    Does mapping tet0 to tet1 via
       *    permutation_by_index[i] define an
       *    isometry of the tilings of H^3?
       */
      mapping_graph=attempt_tiling_isometry(tet0,tet1,permutation_by_index[i]);
      
      if (mapping_graph== NULL) continue; 
      
      if (!permutation_type(permutation_by_index[i])) or_pres=false;

      num_tiling_isometries++;
      
      /*
       *    Check if we've already found this mapping graph. 
       */
      for (mg=head_mapping_graph; mg!=NULL; mg=mg->next_mapping_graph)
	if (same_mapping_graph(mapping_graph, mg)) {
	  free_mapping_graphs(mapping_graph); break; 
	}
      
      if (mg!=NULL) continue; 

      /*
       *    store mapping graph.
       */
      mapping_graph->next_mapping_graph=head_mapping_graph;
      head_mapping_graph=mapping_graph;
      num_mapping_graphs++;
      // if (print_flag) print_mapping_graph(mapping_graph); 
    }                
  }
  
  if (manifold0->solution_type[filled]==geometric_solution)
    store_transformations(manifold0, head_mapping_graph,TR);
  
#if 0
  if (same_manif) {
    for (j=0; j<nc; j++) {
      cout << fix_count[j];
      if (j<nc-1) cout << ':';
    }
    cout << " num mapping graphs " << num_mapping_graphs << endl; 
  }
#endif
  /*
   * free structures
   */

  if (mapping_graphs) *mapping_graphs = head_mapping_graph; 
  else free_mapping_graphs(head_mapping_graph);

  return num_tiling_isometries;
}

struct bc_tet {
  Tetrahedron* T; 
  Permutation P; 
  bc_tet *rep; 

  bc_tet() : rep(0) {}
  bc_tet* rrep();
  int cusp() const { return T->cusp[P & 0x3]->index; }
  friend void set_same_class(bc_tet& a, bc_tet& b); 
};

bc_tet* bc_tet::rrep()
{
  if (!rep) return this; 
  if (!rep->rep) return rep; 
  rep = rep->rrep(); 
  return rep;
}

void set_same_class(bc_tet& a, bc_tet& b)
{
  bc_tet* rra = a.rrep();
  bc_tet* rrb = b.rrep();
  if (rra != rrb) rrb->rep = rra; 
}

void make_bc_tet_array(Triangulation* manifold, vector<bc_tet>& A)
{
  A.resize(get_num_tetrahedra(manifold)*24);

  int i; 
  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next) {
    for (i=0; i<24; i++) {
      bc_tet& X(A[24 * t->index + i]);
      X.T = t; 
      X.P = permutation_by_index[i]; 
    }
  }
}

void set_same_class(int a, int b, vector<int>& classes)
{
  // range check. 
  if (a<0 || a>=classes.size() ||
      b<0 || b>=classes.size()) {
    // printf("index out of range in set_same_class\n"); 
    return; 
  }

  int ca = classes[a], cb = classes[b], i; 
  if (ca==cb) return; 

  // ensure ca < cb.
  if (ca > cb) { i=cb; cb=ca; ca=i; } 
  
  // change cb to ca wherever it appears and decrement all classses 
  // bigger than cb to fill the gap. 
  int n = classes.size();
  for (i=0; i<n; i++) {
    if (classes[i]==cb) classes[i]=ca; 
    else if (classes[i]>cb) classes[i]--; 
  }
}

void quotient_bc_tet_arrays(GraphNode* head_graph, vector<bc_tet>& pr, vector<bc_tet>& im)
{
  // Traverse the mapping graphs merging equivalent bc_tet's. 
  GraphNode* mg;
  int i, pi; 
  for(; head_graph!=NULL; head_graph=head_graph->next_mapping_graph) {
    for(mg=head_graph; mg!=NULL; mg=mg->next_queue) {

      for (i=0; i<24; i++) {
	bc_tet& x(pr[mg->preimage->index * 24 + i]);
	pi = index_by_permutation[compose_permutations(mg->permutation, x.P)];
	set_same_class(x, im[mg->image->index * 24 + pi]); 
      }
    }
  }
}

vector<int> cusp_covering_degrees(int num_cusps, vector<bc_tet>& A)
{
  vector<int> degree(num_cusps, 0); 

  int c, i, N=A.size(); 
  for (c=0; c<num_cusps; c++) {

    // look for one triangle in this cusp. 
    for (i=0; i<N; i++) {
      if (A[i].cusp()==c) break; 
    }

    // check we found one. 
    if (i==A.size()) {
      cerr << "cusp " << c << " not found in bc_tet array\n"; 
      for (i=0; i<num_cusps; i++) degree[i]=0; 
      return degree; 
    }

    // count how many triangles in this cusp have same quotient image. 
    bc_tet* rr = A[i].rrep(); 
    for (; i<N; i++) {
      if (A[i].cusp()!=c) continue; 
      if (A[i].rrep()==rr) degree[c]++; 
    }
  }
  return degree; 
}

vector<int> cusp_hidden_symmetry_classes(int num_cusps, 
					 GraphNode* head_graph)
{
  // First give each cusp its own class. 
  int i; 
  vector<int> classes(num_cusps); 
  for (i=0; i<num_cusps; i++) classes[i]=i; 

  // Traverse the mapping graphs finding cusp mappings. 
  GraphNode* mg;
  int pi; 
  for(; head_graph!=NULL; head_graph=head_graph->next_mapping_graph) {
    for(mg=head_graph; mg!=NULL; mg=mg->next_queue) {

      for (i=0; i<4; i++) {
	pi = EVALUATE(mg->permutation, i); 
	set_same_class(mg->preimage->cusp[i]->index,
		       mg->image->cusp[pi]->index,
		       classes);
      }

    }
  }

  return classes; 
}

vector<int> cusp_covering_degrees(Triangulation *m0, Triangulation *m1, GraphNode* head_graph)
{
  vector<bc_tet> bc0, bc1;

  make_bc_tet_array(m0, bc0); 
  make_bc_tet_array(m1, bc1); 

  quotient_bc_tet_arrays(head_graph, bc0, bc1); 

  return cusp_covering_degrees(get_num_cusps(m0), bc0); 
}

double perimeter(VertexCrossSections* vcs, int i)
{
  int j; 
  double p=0.; 
  for (j=0; j<4; j++) {
    if (j==i) continue; 
    p+= vcs->edge_length[i][j]; 
  }
  return p; 
}

void set_size_ratio(int a, int b, double sa, double sb, n_vector& sizes, vector<int>& classes)
{
  if (a<0 || a>=classes.size() ||
      b<0 || b>=classes.size()) {
    printf("index out of range in set_size_ratio\n"); 
    return; 
  }

  if (classes[a]==classes[b]) {
    if (fabs(sizes[a]/sizes[b] - sa/sb) > 1e-5) 
      cout << "inconsistent cusp cross section sizes found.\n";
    return;
  }

  // Adjust ratio between sizes in class a and sizes in class b. 
  int i, n=classes.size(), cb = classes[b];
  double rat = sb/sa / (sizes[b]/sizes[a]); 
  for (i=0; i<n; i++)
    if (classes[i]==cb) sizes[i] *= rat; 

  set_same_class(a,b,classes); 
}

vector<int> get_cusp_size_vector(Triangulation* manifold, GraphNode* head_graph, 
				 n_vector& sizes)
{
  int num_cusps = get_num_cusps(manifold); 
  sizes = n_vector(num_cusps, 1.); 

  allocate_cross_sections(manifold); 
  compute_cross_sections(manifold); 

  // First give each cusp its own class. 
  int i; 
  vector<int> classes(num_cusps); 
  for (i=0; i<num_cusps; i++) classes[i]=i; 

  // Traverse the mapping graphs finding cusp mappings. 
  GraphNode* mg;
  int pi; 
  double x, y; 
  for(; head_graph!=NULL; head_graph=head_graph->next_mapping_graph) {
    for(mg=head_graph; mg!=NULL; mg=mg->next_queue) {

      for (i=0; i<4; i++) {
	pi = EVALUATE(mg->permutation, i); 

	x = perimeter(mg->preimage->cross_section, i);
	y = perimeter(mg->image->cross_section, pi); 

	set_size_ratio(mg->preimage->cusp[i]->index,
		       mg->image->cusp[pi]->index,
		       y, x, sizes, classes);
      }

    }
  }

  free_cross_sections(manifold); 
  return classes; 
}

bool is_integer_vector(n_vector const& v)
{
  int i;
  for (i=0; i<v.dim; i++) 
    if (v[i]!=floor(v[i]+.5)) return false; 
  return true;
}

// returns the normalized area sum. 

int normalize_independent_cusps(n_vector& areas, vector<int> const& chsc)
{
  int i, c=chsc.size(), mcn=0; // get max class number
  for (i=0; i<c; i++) if (chsc[i] > mcn) mcn = chsc[i]; 

  ++mcn;

  n_vector area_sum(mcn,0.);
  for (i=0; i<c; i++) area_sum[chsc[i]] += areas[i]; 

  if (mcn==1) return int(area_sum[0]);

  if (is_integer_vector(area_sum)) { // it should be
    int prod=1, hcf; 

    // get product of area sums
    for (i=0; i<mcn; i++) {
      prod *= int(area_sum[i]); 
    }
    // normalize areas
    for (i=0; i<c; i++)
      areas[i] *= (prod/int(area_sum[chsc[i]]));

    // remove any common factor. 
    hcf = int(areas[0]); 
    for (i=1; i<c; i++) {
      hcf = gcd(hcf, int(areas[i])); 
    }
    areas /= double(hcf); 

    return prod/hcf; 

  } else {
    cerr << "Warning: cusp area vector contains non-integer values\n"; 

    for (i=0; i<c; i++)
      areas[i] /= area_sum[chsc[i]];
  }
  return 1; 
}



bool has_finite_vertices(
    Triangulation   *manifold)
{
    Cusp    *cusp;

    for (cusp = manifold->cusp_list_begin.next;
     cusp != &manifold->cusp_list_end;
     cusp = cusp->next)
      if (cusp->is_finite) return true;

    return false;
}

#if 0
/******************************************************************************/
void print_tet_shapes(Triangulation *manifold)
/******************************************************************************/
{
  char s[20];
  int n;

  n = get_num_tetrahedra(manifold);

  Tetrahedron *t; 
  for (t = manifold->tet_list_begin.next;
       t != &manifold->tet_list_end;
       t = t->next)
    {
      cout << t->shape[complete]->cwl[ultimate][0].rect << endl;
      cout << " " << t->shape[complete]->cwl[ultimate][1].rect << endl;
      cout << "  " << t->shape[complete]->cwl[ultimate][2].rect << endl;
    }
}
#endif 

/******************************************************************************/
GraphNode    *attempt_tiling_isometry
/******************************************************************************/
(    Tetrahedron *tet0,
     Tetrahedron *tet1,
     Permutation permutation
)
{
  GraphNode *root=NULL,
    *queue=NULL,
    *tail=NULL,
    *node=NULL,
    *map_listing[MAX_TET][MAX_TET];
  Boolean OK;
  Tetrahedron *nbr0=NULL, *nbr1=NULL; 
  FaceIndex face0, face1;
  int i,j,index=1;
  Permutation gluing0, gluing1, nbr0_permutation;

  /*
   *    Check the intial map works.
   */

  OK = same_shape(tet0,tet1,permutation);

  if (!OK) return NULL;
  
  /*
   *    Intialize map_listing.
   */
  for(i=0;i<MAX_TET;i++){
    for(j=0;j<MAX_TET;j++)
      map_listing[i][j]=NULL;
  }


  root=new_node(tet0, tet1, permutation,index++);
  map_listing[tet0->index][tet1->index]=root; 

  
  /*
   *    The queue will hold the nodes of the graph whose neighbors haven't
   *    been checked.
   */
  queue=root;
  tail=root;

  /*
   *    While there are incomplete nodes continue.
   */

  for(;queue!=NULL;queue=queue->next_queue) {
    
    /*
     *    Look at the neighbors of queue->preimage in manifold0
     *  and the map induced on them by
     *    by mapping queue->preimage to queue->image
     *  via queue->permutation.
     */    
    permutation=queue->permutation;
    
    for( face0=0; face0<4; face0++) {
      nbr0=queue->preimage->neighbor[face0];
      face1=EVALUATE(permutation, face0);
      nbr1=queue->image->neighbor[face1];    
      
      /*
       *    Let gluing0 be the gluing which identifies face0 of tet0
       *    to nbr0, and similarly for gluing1.
       */
      
      gluing0=queue->preimage->gluing[face0];
      gluing1=queue->image->gluing[face1];
      
      /*
       *            queue->preimage------------>nbr0
       *                        |   gluing0      |
       *    queue->permutation  |                |    nbr0->permutation
       *                        |                |
       *                        V                V
       *            queue->image------------->nbr1
       *                            gluing1
       *
       *    We want to ensure that image and nbr1 enjoy the same
       *    to each other in manifold1's triangulation that preimage 
       *    and nbr0 do in manifold0's triangulation. 
       */    
      
      nbr0_permutation=compose_permutations(
                compose_permutations(gluing1,permutation),
                        inverse_permutation[gluing0]);

      node=map_list(nbr0->index,nbr1->index,nbr0_permutation,map_listing);

      /*
       *    make sure this map hasn't been seen before, if it has then
       *    then we can link the graph.           
       */
      if (node==NULL) {

	/*
	 *    This hasn't been seen before so check the tetrahedron
	 *    have the same shape.
	 */
	OK = same_shape(nbr0,nbr1,nbr0_permutation); 
	if (!OK) {free_mapping_graphs(root); return NULL; }
	
	node=new_node(nbr0, nbr1, nbr0_permutation, index++);
	
	/*
	 *    Add this new node to the map listing.
	 */
      
	node->next_map=map_listing[nbr0->index][nbr1->index];
	map_listing[nbr0->index][nbr1->index]=node;
            
	/*
	 *    add this to the queue.
	 */
	tail->next_queue=node;
	tail=node;
	
      }

      /*
       *    If we have seen this node before we should link the graph 
       *    back to it.
       */    
      
      queue->nbr[face0]=node;
    }
  }


  /*
   *    If this is a tiling isometry then we should return the root.
   */     
  return root;
}


/******************************************************************************/
GraphNode *map_list
/******************************************************************************/
(int     index0,
 int      index1,
 Permutation     permutation,
 GraphNode *map_listing[MAX_TET][MAX_TET])
{
  GraphNode *cur;

  for (cur=map_listing[index0][index1]; cur!=NULL; cur=cur->next_map){
   if (cur->permutation==permutation)
      return cur;
  }

  return cur;
}


/*****************************************************************************/
Boolean same_shape 
/*****************************************************************************/
( Tetrahedron *tet0,
  Tetrahedron *tet1,
  Permutation permutation
)
{
    int i,j,ii,ji;

    for(i=0; i<4; i++) {
      ii = EVALUATE(permutation,i);
      if (tet0->cusp[i]->is_finite!=
	  tet1->cusp[ii]->is_finite)
	return FALSE;
      for(j=i+1;j<4;j++) {
	ji = EVALUATE(permutation,j);
	if (tet0->edge_class[edge_between_vertices[i][j]]->order !=
	    tet1->edge_class[edge_between_vertices[ii][ji]]->order)
	  return FALSE;
      }
    }
    return TRUE;    

}


/***************************************************************************/
Boolean permutation_type
/***************************************************************************/
( Permutation permutation )
{
  int i,j, counter=0;
  int image[4];

  for (i=0; i<4; i++)
    image[i]=EVALUATE(permutation,i);

  for(i=0;i<3;i++)
    for(j=i+1;j<4;j++)
      if (image[i]>image[j]) counter++;
  
  if ((counter % 2)==0) return TRUE;
  else return FALSE;
}


/******************************************************************************/
void free_mapping_graphs
/******************************************************************************/
( GraphNode *head_graph )
{
 GraphNode  *next_graph,
                      *next_node;

 for(; head_graph!=NULL; head_graph=next_graph){
    next_graph=head_graph->next_mapping_graph;

    for(; head_graph!=NULL; head_graph=next_node){
        next_node=head_graph->next_queue;

        my_free(head_graph); 
    }
 }    

    return;
}


/******************************************************************************/
GraphNode *new_node
/******************************************************************************/
(
Tetrahedron *tet0,
Tetrahedron *tet1,
Permutation permutation,
int index
)
{
    GraphNode *node=NULL;
    int i;

    node = NEW_STRUCT(GraphNode);

    if (node==NULL){
        fprintf(stderr,"ERR:    memory.\n");
        exit(EXIT_FAILURE);
    }

    node->preimage=tet0;
    node->image=tet1;
    node->permutation=permutation;
    node->index=index;
    for(i=0;i<4;i++)
        node->nbr[i]=NULL;
    node->next_queue=NULL;
    node->next_map=NULL;
    node->next_mapping_graph=NULL;

    return node;
}


/******************************************************************************/
void print_mapping_graph
/******************************************************************************/
(
GraphNode *root
)
{
  GraphNode *cur;

  char pm; 
  for (cur=root; cur!=NULL; cur=cur->next_queue) {
    printf("tet %d -> %d ", cur->preimage->index, cur->image->index); 
    print_permutation(cur->permutation); 
    pm = permutation_type(cur->permutation) ? '+' : '-';
    printf(" %c\n", pm);
  }

  printf("\n"); 

  return;
}


/****************************************************************************/
Boolean same_mapping_graph(GraphNode* r1, GraphNode* r2)
/****************************************************************************/
{
  GraphNode *cur;
  for (cur=r2; cur!=NULL; cur=cur->next_queue) {
    if (r1->preimage->index == cur->preimage->index &&
    r1->image->index == cur->image->index &&
    r1->permutation == cur->permutation) 
      break; 
  }

  return cur!=NULL; 
}

/******************************************************************************/
void print_permutation
/******************************************************************************/
(
Permutation permutation
)
{
    int i;

    for(i=0;i<4;i++)
        printf("%d", EVALUATE(permutation,i));

    return;
}


static void next_on_face(Tetrahedron*& tet, VertexIndex& cv, VertexIndex& fv)
{
  Tetrahedron *ntet = tet->neighbor[fv];
  VertexIndex ncv = EVALUATE(tet->gluing[fv], cv);
  VertexIndex nfv = EVALUATE(tet->gluing[fv], 5-fv);
  tet = ntet;
  cv = ncv;
  fv = nfv; 
}

static void edge_partner(Tetrahedron*& tet, VertexIndex& cv, VertexIndex& fv)
{
  Tetrahedron *ptet = tet->neighbor[1-cv];
  VertexIndex pcv = EVALUATE(tet->gluing[1-cv], cv);
  VertexIndex pfv = EVALUATE(tet->gluing[1-cv], 5-fv);
  tet = ptet;
  cv = pcv; 
  fv = pfv;
}

inline static Complex the_corner(Tetrahedron* tet, VertexIndex cv, VertexIndex fv)
{
  return tet->corner[(fv==cv+2) ? fv : cv];
}

static void get_icell_corner_quad(Tetrahedron* tet, VertexIndex cv, VertexIndex fv, Complex c[4])
{
  if (cv > 1) goto error;
  c[0] = the_corner(tet,cv,fv);

  next_on_face(tet,cv,fv);
  if (cv > 1) goto error;
  c[1] = the_corner(tet,cv,fv);

  edge_partner(tet,cv,fv);
  if (cv > 1) goto error;
  c[2] = the_corner(tet,cv,fv);

  next_on_face(tet,cv,fv);
  if (cv > 1) goto error;
  c[3] = the_corner(tet,cv,5-fv);

#if 0
  int i;
  cout << "[";
  for (i=0; i<4; i++) {
    print(c[i]);
    if (i<3) cout << ", ";
  }
  cout << "]\n";
#endif

  return; 
 error:
  fprintf(stderr, "Error in get_icell_corner_quad: central vertex not 0 or 1.\n"); 
}

/******************************************************************************/
static void store_transformations
/******************************************************************************/
(
 Triangulation *manifold,
 GraphNode *root,
 vector<MoebiusTransformation>& TR
)
{
  MoebiusTransformation trans;
  int i,j;
  Complex c[4], pc[4];
  GraphNode *cur;
  int i0,i2;

#if 0
  /* 
   * These should already have been called. 
   * Can't call them now since triangulation might have finite vertices. 
   */
  find_complete_hyperbolic_structure(manifold);
  choose_generators(manifold,TRUE,FALSE);
#endif 

  bool finite_vertices = has_finite_vertices(manifold); 

  for (cur=root; cur!=NULL; cur=cur->next_mapping_graph) {
    if (finite_vertices) {
      /* Since cur->preimage won't be changing, could really do this 
	 outside the loop. */ 
      get_icell_corner_quad(cur->preimage, 0, 2, pc);
      i0 = EVALUATE(cur->permutation,0);
      i2 = EVALUATE(cur->permutation,2);
      get_icell_corner_quad(cur->image,  i0, i2, c);
      trans = get_transform(pc, c); 

#if 0
      cout << "[";
      for (i=0; i<4; i++) {
	print(c[i]);
	if (i<3) cout << ", ";
      }
      cout << "]\n";
      cout << "[";
      for (i=0; i<4; i++) {
	print(trans*pc[i]);
	if (i<3) cout << ", ";
      }
      cout << "]\n";
#endif

    } else {
      for (i=0;i<4;i++) {
	j=EVALUATE(cur->permutation,i);
	c[i]=cur->image->corner[j];
      }
      trans=get_transform(cur->preimage->corner,c);
    }
    TR.push_back(trans); 
  }
}

#if 0
/******************************************************************************/
Boolean is_arithmetic
/******************************************************************************/
(int census, int n)
{
  int i; 
  if (census==5) {
    for (i=0; i<40; i++) 
      if (n==arith5[i]) return TRUE; 
  } else if (census==6) {
    for (i=0; i<27; i++) 
      if (n==arith6[i]) return TRUE; 
  } else if (census==7) {
    for (i=0; i<8; i++) 
      if (n==arith7[i]) return TRUE; 
  } else {
    fprintf(stderr, "There is no %d census!\n", census);
    return FALSE;
  }
  return FALSE;
}
#endif

// scans at most v.dim, return number scanned.  
int scan_ratio_vector(const char* str, n_vector& v)
{
  int i, n=100; 
  if (*str=='(') str++; //')'
  for (i=0; i<v.dim; i++) {
    if (!*str) return i; 
    if (sscanf(str, "%lf %n ", &v[i], &n)!=1) return i; 
    if (n!=100) str+=n; 
    else { while (*str && *str!=':') str++; }
    if (*str==':') str++; 
  }
  return i; 
}

void sprint_vector(char* buf, vector<int> const& v, int w)
{
  char* cp=buf; 
  int i, n=v.size();
  *cp++='[';
  for (i=0; i<n; i++) {
    *cp++=char('0'+v[i]); 
    if (i<n-1) *cp++ = ',';
  }
  *cp++=']';
  while (cp-buf < w) *cp++ = ' '; 
  *cp = '\0'; 
}

void sprint_ratio_vector(char* buf, n_vector const& v)
{
  char* cp=buf;
  int i, nc; 
  for (i=0; i<v.dim; i++) {
    nc = (i<v.dim-1) ? sprintf(cp, "%g:", v[i]) : sprintf(cp, "%g", v[i]); 
    cp += nc; 
  }
}

void set_to_integer_multiple(n_vector& v)
{
  int i;
  long lcd=1, n, d; 
  const double eps=1e-9; // use 9 of 16 digits. 
  for (i=0; i<v.dim; i++) {
    if (!appears_rational(v[i]-eps,v[i]+eps,1e-4,&n,&d)) return; 
    lcd *= (d/gcd(lcd,d)); 
  }
  for (i=0; i<v.dim; i++) {
    appears_rational(v[i]-eps,v[i]+eps,1e-4,&n,&d); 
    v[i] = n*lcd/d;
  }
}

void remove_common_factor(vector<int>& v)
{
  int n = v.size(); 
  if (!n) return; 
  int i, cd=v[0]; 
  for (i=1; cd>1 && i<n; i++) {
    cd = gcd(cd, v[i]); 
  }
  for (i=0; i<n; i++)
    v[i] /= cd; 
}

