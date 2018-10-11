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
#include "tilt_polytope.hh"

using std::cout;
using std::cerr;
using std::endl;

#define EPSILON 1e-6

tp_face::tp_face(const tilt_polytope* tp, n_vector const& tcs)
  : T(tp), tied_cusp_sizes(tcs), cusp_size(T->ties.size()) 
{
  T->expand_with_ties(cusp_size, tied_cusp_sizes); 
}

void tp_face::insert_tilt(n_vector const& new_tilt)
{
  list<n_vector>::const_iterator it;
  for (it=tilt.begin(); it!=tilt.end(); ++it)
    if (*it==new_tilt) return; // already got it, don't need to add it
  tilt.push_back(new_tilt); 
  return;
}

void tilt_polytope::contract_with_ties(n_vector& tilt, const double tilt_vector[]) const
{
  int i, nc=ties.size(); 
  for (i=0; i<tilt.dim; i++) tilt[i]=0.;
  for (i=0; i<nc; i++) {
    tilt[ties[i]] += tilt_vector[i] * tie_ratios[i]; 
  }
}

// tilt_sum should be [0,...,0,1] when called. 

void tp_face::initialize(Triangulation* manifold, n_vector& tilt_sum)
{
  proto_canonize(manifold, cusp_size.V);
  double* tilt_mx = face_tilt_matrix(manifold); 

  int i, nc=T->ties.size(), nf=2*get_num_tetrahedra(manifold); 

  n_vector this_tilt(T->ntc);

  for (i=0; i<nf; i++) {
    T->contract_with_ties(this_tilt, &tilt_mx[i*nc]); 
    if (tied_cusp_sizes * this_tilt < -EPSILON) // check all tilts positive.
      cout << "proto_canonize seems to have failed!\n";
    insert_tilt(n_vector(this_tilt));
    tilt_sum -= this_tilt; 
  }
  my_free_array(tilt_mx); 
}

bool on_boundary(n_vector const& v)
{
  int i; 
  for (i=0; i<v.dim; i++) 
    if (v[i]<EPSILON) return true;
  return false; 
}

void tilt_polytope::expand_with_ties(n_vector& V, n_vector const& OV) const
{
  int i, nc = V.dim; 
  for (i=0; i<nc; i++) V[i] = OV[ties[i]] * tie_ratios[i]; 
}

// find a tied_cusp_size vector outside of the tilt constraints (V.dim == nc)
bool tp_face::locate_outside_vertex(n_vector& V) const
{
  list<n_vector>::const_iterator ti;
  polytope<VS>::IV_p vi;
  int i; 
  for (ti=tilt.begin(); ti!=tilt.end(); ++ti) {
    for (vi=F->I.begin(); vi!=F->I.end(); ++vi) {
      if ((*ti) * ((*vi)->V) < -EPSILON) {
	V = (*vi)->V; 

	// don't want any cusp size to be zero. 
	if (on_boundary(V)) {
	  double ipi = (*ti) * tied_cusp_sizes;
	  double ipo = (*ti) * V;
	  double t = ipo/(2*(ipo-ipi)); // go halfway between crossing point and V.
	  V += t*(tied_cusp_sizes - V); 
	}

	// cout << "\noutside vertex " << ((*vi)->V) << " tilt " << (*ti) << "\n\n";
	for (i=0; i<V.dim; i++) {
	  if (V[i] < 0.) {
	    cout << "Oops, bad cusp_size " << V << endl; 
	    V[i] = .02; 
	  }
	}

	return true; 
      }
    }
  }
  return false; 
}

static vector<int> count(int n)
{
  vector<int> res(n); 
  int i; 
  for (i=0; i<n; i++) res[i]=i; 
  return res; 
}

// assumes tie classes are numbered consecutively from 0. 
static int num_tie_classes(vector<int> const& t)
{
  int i, n=0, m=t.size();
  for (i=0; i<m; i++) if (t[i] > n) n=t[i]; 
  return n+1; 
}

void tilt_polytope::set_default_ties()
{
  ties = count(get_num_cusps(manifold)); 
  tie_ratios = n_vector(get_num_cusps(manifold),1.);
}

tilt_polytope::tilt_polytope(Triangulation* m, vector<int> const& t, n_vector const& r)
  : manifold(m), ties(t), tie_ratios(r)
{
  unsigned nc = get_num_cusps(manifold);

  // check on ties. 
  if (ties.size()==0 && tie_ratios.dim==0) {
    set_default_ties(); 
  } else if (ties.size()!=nc || tie_ratios.dim!=nc) {
    cout << "invalid tie specification in tilt_polytope constructor.\n";
    set_default_ties(); 
  } 
  ntc = num_tie_classes(ties);
  P.initialize(ntc+1); 
}

tilt_polytope::tilt_polytope(Triangulation* m)
  : manifold(m)
{
  ntc = get_num_cusps(manifold);
  set_default_ties(); 
  P.initialize(ntc+1); 
}

void tilt_polytope::add_face(n_vector const& tied_cusp_sizes)
{
  FL.push_back(tp_face(this, tied_cusp_sizes)); 

  n_vector tilt_sum(ntc+1,ntc); // tilt_sum = [0,...0,1];
  FL.back().initialize(manifold, tilt_sum); //gets the tilt_sum for this face
  P.cut(tilt_sum, &(FL.back().F), &dead_F); 
}

bool tilt_polytope::compute(int report, int limit)
{
  // for diagnostic printout from add_face.
  // int oldprec = cout.precision(15);

  // create initial face
  n_vector tied_cusp_sizes(ntc, 1.0);
  add_face(tied_cusp_sizes); 

  // add more faces to the tilt polytope until it is complete
  list<tp_face>::iterator fi;
  for (fi=FL.begin(); fi!=FL.end() && limit>0;) {
    if (!fi->F->keep) {
      FL.erase(fi++); // this face is already dead
      continue; 
    }
    if (!fi->locate_outside_vertex(tied_cusp_sizes)) {
      ++fi;
      continue; // this face is complete now
    }
    // add a new tp_face for the new tied_cusp_sizes vector
    // if (report) cout << "adding face for size_vector: " << tied_cusp_sizes << endl; 
    add_face(tied_cusp_sizes);
    --limit;
    // don't increment fi, check this face again. 
  }

  return limit > 0; 
  
  // cout.precision(oldprec); 
}

void tilt_polytope::print() const
{
  cout << "Boundary vertices:";
  n_vector org(ntc+1, ntc);
  polytope<VS>::vertex_cp vp; 
  for (vp = P.VL.begin(); vp!=P.VL.end(); ++vp) {
    if (!on_boundary(vp->V)) continue; 
    if (vp->V==org) continue; 
    cout << ' ' << vp->index;
  }
  cout << endl; 

  cout << "Faces:";
  list<tp_face>::const_iterator fi;
  for (fi=FL.begin(); fi!=FL.end(); ++fi) {
    cout << ' ';
    fi->F->print_incidences(cout); 
  }
  cout << endl; 
}

static const Complex Z[] = { Complex(0., 3.464101615137754587054892682),
			     Complex(-2.,0), 
			     Complex(2.,0) }; 

static Complex to_complex(n_vector const& v)
{
  double sum=0.;
  Complex z(0.,0.); 
  int i; 
  for (i=0; i<3; i++) {
    sum += v[i];
    z += v[i] * Z[i]; 
  }
  return z/sum;
}

class ArgSort {
  Complex o;
public:
  ArgSort(Complex org) : o(org) {}

  bool operator() (const Complex& a, const Complex& b) const; 
};

bool ArgSort::operator() (const Complex& a, const Complex& b) const
{
  return (a-o).arg() < (b-o).arg(); 
}

void tp_face::save_picture(picfile& pic) const
{
  n_vector v(3); 
  polytope<VS>::IV_p vp; 
  int count=0;
  Complex z, mid(Zero); 
  list<Complex> vl;

  for (vp = F->I.begin(); vp != F->I.end(); vp++) {
    v.dehomogenized_copy((*vp)->V);
    z = to_complex(v); 
    mid += z; 
    vl.push_back(z);
    count++;
  }
  if (count < 3) return; // face < 2-d. (error?)
  mid = (mid/(1.0*count));
  // cout << mid << endl; 
  ArgSort cmp(mid);
  vl.sort(cmp); // sort into ccw order around mid. 

  // print the vertices for debug first. 
  list<Complex>::const_iterator i; 
  cout << '{'; 
  for (i=vl.begin(); i!=vl.end(); i++) {
    if (i!=vl.begin()) cout << ", "; 
    cout << (*i);
  }
  cout << '}' << endl;

  vl.push_back(vl.front()); // close it up
  pic.print_line(vl ,Zero);
}

void tilt_polytope::save_picture(picfile& pic) const
{
  if (ntc != 3) {
    cerr << "I can only draw 3-cusp tilt polytopes!\n"; 
    return;
  }

  list<tp_face>::const_iterator f;
  for (f=FL.begin(); f!=FL.end(); f++) 
    f->save_picture(pic);
}

void tilt_polytope::get_cusp_sizes(n_vector& cusp_sizes, n_face const& f) const
{
  n_vector bc;
  f.get_dh_barycenter(bc); 
  expand_with_ties(cusp_sizes, bc); 
}

