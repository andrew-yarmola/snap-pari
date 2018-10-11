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
#include "find_commens.hh"
#include "commens.hh"
#include "helpers.hh"
#include "kernel_extras.hh"
#include "tilt_polytope.hh"

using std::cout;
using std::endl;

bool validate_class_spec(vector<int> const& cs)
{
  int i, n=cs.size(), j=-1;
  for (i=0; i<n; i++) {
    if (cs[i]>j) {
      j++; 
      if (cs[i]!=j) return false; 
    }
    if (cs[i]<0) return false; 
  }
  return true; 
}

static vector<int> def_cusp_ties(int nc)
{
  vector<int> ties(nc);
  int i; 
  for (i=0; i<nc; i++) ties[i]=i;
  return ties;
}


static void print_ratio_vector(vector<int> const& v)
{
  int i, nc=v.size(); 
  for (i=0; i<nc; i++) {
    if (i>0) cout << ':'; 
    cout << v[i]; 
  }
}

int bin_coeff(int n, int k)
{
  int l=n-k; 
  int bc=1; 
  for (;n>l;n--) bc*=n;
  for (;k>1;k--) bc/=k;
  return bc; 
}

static bool commensurable(Triangulation* m, n_vector const& csm, 
			  Triangulation* n, n_vector const& csn)
{
  bool m_copied, n_copied; 
  if (!my_canonize(&m, csm, m_copied) || !my_canonize(&n, csn, n_copied)) {
    cout << "Canonize failed\n"; 
    return false; 
  }

  vector<MoebiusTransformation> TR;
  bool or_pres;
  int nc = 0;

  // copied <=> non-tetrahedral cells so only possible if both agree. 
  if (m_copied==n_copied) 
    nc = commensurabilities(m, n, TR, or_pres); 

  if (m_copied) free_triangulation(m); 
  if (n_copied) free_triangulation(n); 

  return nc > 0;
}

void check_commensurability(i_triangulation& m, n_vector const& csm, i_triangulation& n, n_vector const& csn)
{
  if (num_filled_cusps(m.M()) > 0 || num_filled_cusps(n.M()) > 0) {
    cout << "Manifolds may not have any filled cusps\n"; 
    return; 
  }

  if (commensurable(m.M(), csm, n.M(), csn))
    cout << "Manifolds are commensurable\n";
  else {
    cout << "No common cover found\n";
  }
}

bool next_partition(vector<int>& p)
{
  // find a last part not equal to one. 
  int j, k=p.size(); 
  for (j=k-1;j>0;j--) {
    if (p[j]>1) break; 
  }
  if (j==0) return false; 

  p[j-1]++; 
  p[k-1] = p[j]-1; 
  for (;j<k-1;j++) p[j]=1; 
  return true; 
}

// prints stuff out if report is nonzero.
// 0x1 - print index, chirality, chsc and area vector.
// 0x2 - print generators
// 0x4 - refer to output as commensurator (rather than hidden symmetry gp)

static int compute_commensurator(Triangulation* m, n_vector const& cusp_sizes, int report, vector<int>* csc=0)
{
  int c = get_num_cusps(m); 

  bool copied; 
  if (!my_canonize(&m, cusp_sizes, copied)) {
    cout << "Canonize failed\n"; 
    return c; 
  }

  vector<MoebiusTransformation> TR;
  bool or_pres;
  GraphNode* mg; 
  int nc = commensurabilities(m, m, TR, or_pres, &mg); 

  vector<int> chsc = cusp_hidden_symmetry_classes(c, mg); 

  int i; 
  if (report & 0x2) {
    cout << "Generators: \n";
    for (i=0; i<TR.size(); i++)
      cout << TR[i] << endl; 
  }
  if (report & 0x1) {
    int chmax=0;
    for (i=0; i<chsc.size(); i++) 
      if (chsc[i] > chmax) chmax = chsc[i]; 

    if (report & 0x4 || chmax == 0)
      cout << "Commensurator: index " << nc; 
    else 
      cout << "Hidden symmetry group: index " << nc; 

    if (or_pres) 
      cout << ", chiral\n";
    else
      cout << ", amphicheiral\n";

    cout << "Cusp hidden symmetry classes: "; 
    print_vector(cout, chsc); 
    cout << endl; 

    n_vector sizes(c);
    vector<int> degrees = cusp_covering_degrees(m, m, mg); 
    remove_common_factor(degrees); 
    for (i=0; i<c; i++) sizes[i] = degrees[i]; 

    char buf[1000]; 
    sprint_ratio_vector(buf, sizes); 
    cout << "Cusp area ratios: " << buf << endl; 
  }
  free_mapping_graphs(mg); 

  int mcn=0; // get max class number
  for (i=0; i<c; i++) if (chsc[i] > mcn) mcn = chsc[i]; 

  if (copied) free_triangulation(m); 
  if (csc) *csc = chsc; 
  return mcn; 
}

typedef polytope<VS>::nface_list nface_list;

static bool seek_commensurator(i_triangulation const& m, n_vector& cusp_sizes, int report)
{
  int nc = get_num_cusps(m.M()); 
  cusp_sizes = n_vector(nc, 1.0); 

  if (num_filled_cusps(m.M()) > 0) {
    printf("Manifold may not have filled cusps\n"); 
    return false; 
  }

  if (nc==1) return true; 

  if (!m.tp()) {
    printf("seek_commensurator requires a valid tilt polytope\n"); 
    return false; 
  }

  if (m.tp()->FL.size() < 2) return true; // probably never happens

  // make a copy so original is not changed
  Triangulation* t;
  copy_triangulation(m.M(),&t);

  // get the lattice of faces of all dimensions
  int ntc = m.tp()->ntc; 
  vector<nface_list> L(ntc+1); // will have face lists for dim= 0,...,ntc
  m.tp()->P.get_lattice(L); 

  // seek triangulation which minimizes the number of cusp commensurability classes
  int d, mcc, minmcc=nc;
  vector<int> chsc; 
  n_vector barycenter(ntc); 
  nface_list::const_iterator i, best; 
  for (d=ntc-1; d>=0 && minmcc>0; d--) { // do faces of dimension d
    if (report) cout << "checking " << L[d].size() << " faces of dimension " 
		     << d << endl; 
    for (i=L[d].begin(); i!=L[d].end(); ++i) {
      m.tp()->get_cusp_sizes(cusp_sizes, *i);
      if (on_boundary(cusp_sizes)) continue; 
      mcc = compute_commensurator(t, cusp_sizes, 0, &chsc);
      if (mcc < minmcc) { 
	minmcc = mcc; 
	best=i; 
	if (report) { 
	  cout << "CHSC: "; print_vector(cout, chsc); 
	  cout << " at " << (*i) << endl; 
	}
	if (minmcc==0) { 
	  break; 
	}
      }
    }
  }

  // free the copy
  free_triangulation(t);

  // check we found something. 
  if (minmcc==nc) {
    printf("Failed to compute a single valid commensurability!\n");
    return false; 
  }

  // call compute commensurator again to print out the best result. 
  // if (report) cout << "max symmetry at face " << (*best) << endl; 
  m.tp()->get_cusp_sizes(cusp_sizes, *best);
  return true; 
}


void brute_force_hsymms(i_triangulation& m, int n, int N, vector<int> const& ties, vector<int> const& areas)
{
  int c = get_num_cusps(m.M()); 

  // set C to the number of tie classes. 
  int j; 
  int C=0; 
  for (j=0; j<c; j++) if (ties[j]>C) C=ties[j]; 
  C++; 

  if (n < C) n=C; 
  if (N < n) N=n; 

  vector<int> A(C), AOpt, chsc, Aexp(c); 

  n_vector S(c);
  int mcn, min_mcn = c-1;

  while (n <= N) {

    // set up the first area vector for this n. 
    for (j=0; j<C-1; j++) A[j]=1; 
    A[C-1] = n-(C-1); 

    cout << "Trying " << bin_coeff(n-1, C-1) << " sum " << n <<
      " area vectors\n"; 

    do {
      for (j=0;j<c;j++) {
	Aexp[j] = A[ties[j]]*areas[j];
	S[j] = sqrt(double(Aexp[j])); 
      }
      mcn = compute_commensurator(m.M(), S, 0, &chsc); 
      if (mcn < min_mcn) {
	cout << "CHSC: ";
	print_vector(cout, chsc); 
	cout << " at "; 
	print_ratio_vector(Aexp); 
	cout << endl; 
	
	min_mcn = mcn; 
	AOpt = A; 
	if (mcn == 0) break; 
      } 
    } while (next_partition(A)); 

    if (mcn == 0) break; 
    n++; 
  }

  if (min_mcn==c-1) cout << "No quotient found with fewer cusps\n";
  else {
    for (j=0;j<c;j++) S[j] = sqrt(double(AOpt[ties[j]]*areas[j])); 
    compute_commensurator(m.M(), S, 1);
  }
}

void print_hidden_symmetries(i_triangulation& m, n_vector const& csm)
{
  if (num_filled_cusps(m.M()) > 0) {
    cout << "Manifolds may not have any filled cusps\n"; 
    return; 
  }

  compute_commensurator(m.M(), csm, 3); 
}

void find_commensurator(i_triangulation& m)
{
  if (num_filled_cusps(m.M()) > 0) {
    cout << "Manifold may not have any filled cusps\n"; 
    return; 
  }

  Triangulation* t=m.M();
  int nc = get_num_cusps(t); 
  n_vector cusp_sizes(nc, 1.0); 

  // can we get out quickly? 
  if (nc==1) {
    compute_commensurator(t, cusp_sizes, 1); // cusp_sizes=[1], report=1. 
    return;
  }

  if (!m.tp()) { 
    cout << "Computing tilt polytope...\n"; 
    m.compute_tilt_polytope(def_cusp_ties(nc), cusp_sizes); 
    cout << "Tilt polytope has " << (m.tp()->FL.size()) << " faces\n";
  } else {
    cout << "Using existing tilt polytope..\n"; 
  }

  seek_commensurator(m, cusp_sizes, 1);

  compute_commensurator(t, cusp_sizes, 5);
}
