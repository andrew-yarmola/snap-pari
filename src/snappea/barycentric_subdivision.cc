/*
** Copyright (C) 2003 Oliver A. Goodman <oag@ms.unimelb.edu.au>
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

#include "kernel.h"
#include <cstdlib>

using namespace std;

static Permutation swap_two(Permutation P, int k)
{
  unsigned char mask  = (0x03 << (2*k));
  unsigned char mask2 = (mask << 2);
  unsigned char bits  = P & mask;
  unsigned char bits2 = P & mask2;
  P -= (bits + bits2); 
  bits  <<= 2; 
  bits2 >>= 2;
  P += (bits + bits2);
  return P; 
}

static Permutation perm_product(Permutation B, Permutation A)
{
  unsigned char C=0x0;
  int i; 
  for (i=0; i<4; i++) {
    C >>= 2; 
    C |= ((B >> (2*(A & 0x3))) & 0x3) << 6;
    A >>= 2;
  }
  return C; 
}

// Low or high precision COMPLEX?
#ifndef NOPARI
typedef pari COMPLEX;
#else
typedef Complex COMPLEX;
#endif

struct vrep
{
  MoebiusTransformation M;
  vrep* rep;
  COMPLEX P; 
  
  vrep() : M(One), rep(0) { }
  
  COMPLEX compute_P(); 

  vrep* root();
};

vrep* vrep::root()
{
  if (!rep) return this;
  vrep* the_rep = rep;
  while (rep->rep) { // Compress the tree of reps.
    M = M * rep->M;
    rep = rep->rep; 
  }
  return rep;
}

// here set a = M * b. 
static void identify(vrep* a, vrep* b, MoebiusTransformation M = One)
{
  vrep* ar = a->root(); 
  vrep* br = b->root();

  if (ar==br) return; // Already in same class, hope M is right.

  ar->rep = br; 
  ar->M = inverse(a->M) * M * b->M; 
}

COMPLEX vrep::compute_P()
{
  if (!rep) {
    // printf("root rep\n"); 
    // P.print(); 
    return P; 
  }
  root(); // Sets rep to root(). 
  // M.acc_matrix.print();
  // rep->P.print(); 
  return M * rep->P; 
}

#ifdef NOPARI

static double unit_random()
{
  static double rmax = 1073741823.5;
  return double(random())/rmax - 1.0; 
}

static Complex random_on_sphere()
{
  double r, theta, x, y, z; 
  z = unit_random(); 
  r = sqrt(1.0 - z*z)/(1.0-z); 
  theta = PI * unit_random(); 
  x = r*cos(theta); 
  y = r*sin(theta); 
  return Complex(x,y); 
}

// Shape is bad if flat or degenerate
static bool bad_shape(Complex cr) {
  const double small = 1e-5;
  Complex z = Two/(cr+I)+I; // Transform to unit circle to test
  return fabs(complex_modulus(z)-1.0) < small; 
}

#else 

#include "pariwrap.hh"

static Complex to_complex(const pari& p)
{
  static pari big = pow(real(10),integer(10)); 
  if (gabs(p) > big) {
    return Infinity; 
  }

  Complex z;
  z.real = greal(p).double_value(); 
  z.imag = gimag(p).double_value(); 
  return z; 
}

static pari unit_random()
{
  static pari rmax = (gpow(pTWO,integer(31))-pONE)/pTWO;
  return real(random())/rmax - pONE;
}

static pari random_on_sphere()
{
  pari r, theta, x, y, z; 
  z = unit_random(); 
  r = sqrt(pONE - z*z)/(pONE-z); 
  theta = pPI * unit_random(); 
  x = r*cos(theta); 
  y = r*sin(theta); 
  return complex(x,y); 
}

pari cross_ratio(pari const& z)
{
  int i;
  static pari big = pow(real(10), integer(30));
  // z.print();
  // printf("CR:\n"); 
  for (i=0; i<4; i++) {
    if (gabs(z[i]) > big) break; // z[i]=Infinity.
  }
  // printf("CR:OK\n");

  pari cr; 
  switch (i) {
  case 4: // Nothing at infinity. 
    cr = (z[3] - z[1]) * (z[2] - z[0]) / ((z[3] - z[0]) * (z[2] - z[1]));
    break; 
  case 0:
    cr = (z[3] - z[1]) / (z[2] - z[1]);
    break;
  case 1:
    cr = (z[2] - z[0]) / (z[3] - z[0]);
    break;
  case 2:
    cr = (z[3] - z[1]) / (z[3] - z[0]);
    break;
  case 3:
    cr = (z[2] - z[0]) / (z[2] - z[1]);
    break;
  }
  return cr; 
}

// Shape is bad if flat or degenerate
static bool bad_shape(pari const& cr) {
  static pari small = pow(real(10),integer(-5));
  pari z = pTWO/(cr+pI)+pI; // Transform to unit circle to test
  // printf("BS:"); z.print(); 
  return gabs(gabs(z)-pONE) < small; 
}
#endif

// static void set_shape(Tetrahedron*, Complex z);
//
// should be defined here and used below -- except
// that we don't actually need the low precision
// shapes at present, and the shape structures take 
// up loads of extra memory. Omitting them means 
// we won't be able to call choose_generators on the
// subdivision. Instead we set up all generator 
// related structures from the input manifold. 

// Quick & dirty, no bound checking. 
class tet_queue {
  Tetrahedron** Q;
  int b, e; 
public:
  tet_queue(int sz) : b(0), e(0) 
  { Q = NEW_ARRAY(sz, Tetrahedron*); }
  ~tet_queue() { my_free_array(Q); }

  void push(Tetrahedron* p) { Q[e] = p; e++; }  
  Tetrahedron* pop() { return (e>b) ? Q[b++] : 0; }
  int size() const { return e-b; }
};

// The following is choose_generators 'lite'. 
// The manifold will already have generator statuses
// and vertex positions worked out. All we need
// to add are the generator paths. We do this rather
// that calling choose_generators because we want the
// vertex positions to stay where they are and 
// we don't want to have to set up the TetShape structures
// which choose_generators requires if it is to position
// vertices. 

static void setup_generator_paths(Triangulation* manifold)
{
  tet_queue Q(manifold->num_tetrahedra); 

  Tetrahedron* tet = manifold->tet_list_begin.next;
  Tetrahedron* nbr;
  int f; 

  // generator paths are all -1 for the base tet. 
  for (f=0; f<4; f++) {
    tet->generator_path = -1; 
  }

  Q.push(tet); 
  while (Q.size()) {
    tet = Q.pop(); 

    // check for neighbors with unassigned generator_paths. 
    for (f=0; f<4; f++) {
      if (f==tet->generator_path) continue; 
      if (tet->generator_status[f]!=not_a_generator) continue; 
      nbr = tet->neighbor[f]; 
      if (nbr->generator_path!=-2) continue; 

      // set path to face pointing back to tet. 
      nbr->generator_path = EVALUATE(tet->gluing[f],f); 
      Q.push(nbr); 
    }
  }
}

#define VR(T,v) ((vrep*)T->extra)[v]


// Normally vertices of each tetrahedron in a triangulation
// are ordered such that a shape parameter in the UHP corresponds
// to a positive volume contribution. For an oriented manifold
// this vertex ordering will be globally consistent with an 
// orientation. For an unoriented or nonorientable manifold 
// it won't be but the (sign of the) ordering is still forced 
// by requiring geometric tetrahedra to have UHP shapes. 
//
// In our barycentric subdivision, the goal is to obtain a
// very specific vertex ordering, so UHP shapes no longer need
// provide positive volume contributions. Since the volume
// code will now need to know which sign to take, we record 
// it in the has_correct_orientation flag. 


Triangulation* barycentric_subdivision(Triangulation* manifold)
{
  choose_generators(manifold, TRUE, FALSE); 

  int n = get_num_tetrahedra(manifold);
  int nc= get_num_cusps(manifold); 
  int ng= get_num_generators(manifold); 

  // Get the matrix generators which we'll need later. 
  MoebiusTransformation* gen = NEW_ARRAY(ng, MoebiusTransformation); 
  matrix_generators(manifold, gen, FALSE, TRUE); 
#ifndef NOPARI
  if (gen[0].acc_matrix.type()==t_INT) {
    printf("accurate info must be set for call to barycentric_subdivision");
  }
#endif

  int i, j, k, v, v0, v3, ni, nj, m, N = 24*n, gi;
  Permutation P, nP;
  MoebiusTransformation M; 
  GeneratorStatus gs;

  // Prepare the new triangulation. 
  Triangulation* new_m = NEW_STRUCT(Triangulation);
  initialize_triangulation(new_m); 
  for (i=0; i<2; i++) 
    new_m->solution_type[i] = manifold->solution_type[i];
  new_m->orientability  = manifold->orientability; 
  new_m->num_generators = ng; 

  // Create N = 24*n Tetrahedra. 
  Tetrahedron** T = NEW_ARRAY(N, Tetrahedron*); 
  for (m=0; m<N; m++) {
    T[m] = NEW_STRUCT(Tetrahedron); 
    initialize_tetrahedron(T[m]); 
    INSERT_BEFORE(T[m], &new_m->tet_list_end);

    // Initialize Extra field, used for computing vertex positions. 
    T[m]->extra = (extra*)new vrep[4];

#ifndef NOPARI 
    T[m]->acc_corners = rvector(4); 
#endif
  }
  new_m->num_tetrahedra = N; 

  // Make copies of the original cusps. 
  Cusp* cusp; 
  Cusp** C = NEW_ARRAY(nc, Cusp*);
  for (i=0; i<nc; i++)
    C[i] = NEW_STRUCT(Cusp);
  for (cusp = manifold->cusp_list_begin.next; 
       cusp != &manifold->cusp_list_end; 
       cusp = cusp->next) {
    i = cusp->index; 
    if (i>=nc || i<0) { // Check for bad cusp index. 
      printf("cusp->index out of range in barycentric_subdivision\n");
      free_triangulation(new_m); return 0; 
    }
    // Copy everything, then overwrite pointers.
    *C[i] = *cusp; 
    INSERT_BEFORE(C[i], &new_m->cusp_list_end); 
  }
  new_m->num_cusps = nc; 
  new_m->num_or_cusps = manifold->num_or_cusps; 
  new_m->num_nonor_cusps = manifold->num_nonor_cusps; 

  // printf("1. got to here\n"); 

  // Here begins the real work. 
  Tetrahedron* tet; 
  for (tet = manifold->tet_list_begin.next; // Each original tet. 
       tet != &manifold->tet_list_end;
       tet = tet->next) {
    i = tet->index; 

#if 0
    // Check all the face pairings do what we expect. 
    for (k=0; k<4; k++) {
      gs = tet->generator_status[k]; 
      if (gs!=not_a_generator) {
	gi = tet->generator_index[k]; 
	for (v=0; v<4; v++) {
	  if (v==k) continue; 
	  tet->acc_corners[v].print(); 
	  // fwprint(stdout, tet->corner[v], 8); printf(" "); 
	}
	printf("\n%d:\n", gi); 
	for (v=0; v<4; v++) {
	  if (v==k) continue;
	  M = (gs==inbound_generator) ? gen[gi] : inverse(gen[gi]); 
	  (M * tet->acc_corners[v]).print(); 
	  // fwprint(stdout, M * tet->corner[v], 8); printf(" "); 
	}
	printf("\n"); 
      }
    }
#endif

    for (j=0; j<24; j++) { // Each tet in subdivision. 
      m = 24*i+j;

      // Set all gluings to the identity permutation. 
      for (k=0; k<4; k++) { // Each face. 
	T[m]->gluing[k] = IDENTITY_PERMUTATION;
      }

      // Set the volume contribution flag. 
      P = permutation_by_index[j];
      T[m]->has_correct_orientation = (parity[P]==0); // 0==even.

      // Compute neighbors on internal faces. 
      for (k=0; k<3; k++) { // each internal face.
	nj = index_by_permutation[swap_two(P,k)];
	T[m]->neighbor[k] = T[24*i+nj];

	// Set internal vertex identifications.
	for (v=0; v<4; v++) {
	  if (v==k) continue; 
	  identify(&VR(T[m],v), &VR(T[24*i+nj],v)); 
	}
      }

      // Compute external neighbor.
      v3 = EVALUATE(P, 3); 
      ni = tet->neighbor[v3]->index; 
      nP = perm_product(tet->gluing[v3],P);
      nj = index_by_permutation[nP];
      T[m]->neighbor[3] = T[24*ni+nj];

      // Get the face pairing transformation. 

      // Setup the generator information on T[m]. 
      gs = tet->generator_status[v3];
      gi = tet->generator_index[v3];
      for (k=0; k<3; k++) {
	T[m]->generator_status[k] = not_a_generator; 
      }
      T[m]->generator_status[3] = gs; 
      T[m]->generator_index[3] = gi;
      T[m]->generator_parity[3] = tet->generator_parity[v3]; 

      // Set external vertex identifications.
      for (v=0; v<3; v++) {
	// Outbound generator maps nbr vertex to vertex. 
	if (gs==not_a_generator) {
	  identify(&VR(T[m],v), &VR(T[24*ni+nj],v));
	} else if (gs==outbound_generator) {
	  identify(&VR(T[m],v), &VR(T[24*ni+nj],v), gen[gi]); 
	} else {
	  identify(&VR(T[24*ni+nj],v), &VR(T[m],v), gen[gi]); 
	}
      }

      // Position the 0-vertices.
      v0 = EVALUATE(P,0); 
#ifndef NOPARI
      VR(T[m],0).P = tet->acc_corners[v0]; 
      // VR(T[m],0).P.print(); 
#else
      VR(T[m],0).P = tet->corner[v0];
#endif

      // Assign 0-vertices their appropriate cusps.
      T[m]->cusp[0] = C[tet->cusp[v0]->index]; 
    }
  }

  // printf("2. got to here\n"); 

  // Check we have the right number of vertex reps.
  int count[4]; // Number of reps by bc dimension. 
  for (v=0; v<4; v++) count[v]=0; 
  for (m=0; m<N; m++) {
    for (v=0; v<4; v++) {
      if (VR(T[m],v).rep==0) count[v]++; 
    }
  }
  bool err = false; 
  if (count[0]!=get_num_cusps(manifold)) {
    printf("Number of 0-vertices wrong in bc-subdivision.\n"); err=true; }
  if (count[1]!=n  ) { // num_edges == num_tetrahedra == n.
    printf("Number of 1-vertices wrong in bc-subdivision.\n"); err=true; }
  if (count[2]!=2*n) { // num_faces == 2*num_tetrahedra.
    printf("Number of 2-vertices wrong in bc-subdivision.\n"); err=true; }
  if (count[3]!=n  ) { // num_tetrahedra == n.
    printf("Number of 3-vertices wrong in bc-subdivision.\n"); err=true; }
  if (err) {
    for (v=0; v<4; v++) {
      printf("%d-vertices: %d\n", v, count[v]);
    }
  }
  if (err) {
    free_triangulation(new_m); 
    return 0;
  }

  create_edge_classes(new_m);  

  setup_generator_paths(new_m); 

  // Setup the peripheral curves (already 0 by initialize_tetrahedron). 
  int ml, rl, s, c, ca, cb, c1, c2, v1, v2, m1, m2;
  for (tet = manifold->tet_list_begin.next; // Each original tet. 
       tet != &manifold->tet_list_end;
       tet = tet->next) {
    i = tet->index; 

    for (ml=0; ml<2; ml++) {
      for (rl=0; rl<2; rl++) {
	for (v=0; v<4; v++) {
	  for (s=0; s<4; s++) {
	    if (s==v) continue; 

	    // Side s on vertex v is split into two parts.
	    // Find the two sub-tetrahedra on side s. 
	    v1 = remaining_face[v][s];
	    v2 = remaining_face[s][v];
	    P  = CREATE_PERMUTATION(0,v,1,v1,2,v2,3,s);
	    nP = CREATE_PERMUTATION(0,v,1,v2,2,v1,3,s);
	    m1 = 24*i + index_by_permutation[P];
	    m2 = 24*i + index_by_permutation[nP];

	    // We need the sum of the flows on these parts to equal c. 
	    c  =  tet->curve[ml][rl][v][s];
	    ca = T[m1]->curve[ml][rl][0][3];
	    cb = T[m2]->curve[ml][rl][0][3];
	    if (c != ca+cb) { // Maybe set up already by neighbor. 

	      // Split the flow roughly in halves. 
	      ca = c/2;
	      cb = (c+(c>0 ? 1:-1))/2;

	      // Assign flows to the sides and their neigbors. 
	      T[m1]->curve[ml][rl][0][3] = ca;
	      T[m1]->neighbor[3]->curve[ml][rl][0][3] = -ca;
	      T[m2]->curve[ml][rl][0][3] = cb; 
	      T[m2]->neighbor[3]->curve[ml][rl][0][3] = -cb;
	    }

	    // Get the other two flows on the original triangle. 
	    c1 = tet->curve[ml][rl][v][v1];
	    c2 = tet->curve[ml][rl][v][v2];

	    // Setup curves for internal sides. 
	    T[m1]->curve[ml][rl][0][2] = FLOW(c2,c);
	    T[m1]->curve[ml][rl][0][1] = -(ca + FLOW(c2,c)); 
	    T[m2]->curve[ml][rl][0][2] = FLOW(c1,c);
	    T[m2]->curve[ml][rl][0][1] = -(cb + FLOW(c1,c));
	  }
	}
      }
    }
  }

  // printf("3. got to here\n"); 

#if 0
  for (m=0; m<N; m++) {
    VR(T[m],0).P.print();
  }
#endif

  // Position the vertices and compute shapes.
#ifndef NOPARI
  pari cr; 
#else 
  Complex cr;
#endif
  int tries = 0; 
  bool bad = false; 
  while (tries < 3) {
    ++tries; 

    // Randomly assign positions to the 1,2 and 3-vertices. 
    for (m=0; m<N; m++) {
      for (v=1; v<4; v++) {
	if (VR(T[m],v).rep) continue; // Only position root reps. 
	VR(T[m],v).P = random_on_sphere(); 
	// VR(T[m],v).P.print(); 
      }
    }
    // Compute the vertex positions and tet shapes. 
    for (m=0; m<N; m++) {
#ifndef NOPARI
      for (v=0; v<4; v++) {
	T[m]->acc_corners[v] = VR(T[m],v).compute_P(); 
	T[m]->corner[v] = to_complex(T[m]->acc_corners[v]); 
      }
      // T[m]->acc_corners.print();
      cr = cross_ratio(T[m]->acc_corners);
      bad = bad_shape(cr);
      // printf("BS:OK\n");
      if (bad) break; // Avoid flat tetrahedra.
      T[m]->acc_shape = cr;
      // set_shape(T[m], to_complex(cr)); 

#else
      for (v=0; v<4; v++) {
	T[m]->corner[v] = VR(T[m],v).compute_P(); 
      }
      cr = cross_ratio(T[m]->corner);
      if (bad_shape(cr)) break; // Avoid flat tetrahedra.
      // set_shape(T[m], cr);
#endif
    }
    if (m==N) break; // No flat tetrahedra. 
  }

  // Tidy up. 
  for (m=0; m<N; m++) {
    my_free_array((vrep*)T[m]->extra); 
    T[m]->extra = 0;
  }

  my_free_array(C);
  my_free_array(T); 
  my_free_array(gen);
 
  if (tries==3) {
    printf("Couldn't position finite vertices in barycentric_subdivision\n");
    free_triangulation(new_m); return 0; 
  }

  return new_m; 
}



