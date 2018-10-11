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
#include "ortholine.hh"
#include <iomanip>

using std::cerr;
using std::endl;
using std::setw;
using std::swap;
using std::vector;

// static say xx("ortholine.cc"); 

// Ortholine coordinates express the position of a directed geodesic in H^3 relative
// to the geodesic [0,infinity]. The distance is the complex distance from [0,infinity]
// to l, and the position is the complex distance along [0,infinity] 
// from [-1,1] to the common orthogonal between [0,infinity] and l.

void get_ortholine_coordinates(const line& l, Complex& distance, Complex& position)
{
  double eps = 1e-7; 

  // Deal with all the degenerate cases first. 
  if (same_point(l.end[0],l.end[1],eps)) {
    position = complex_log(l.end[0], 0.0);
    distance = Infinity; 
    return; 
  }
  if (same_point(l.end[0],Zero,eps)) {
    distance = Zero; 
    if (same_point(l.end[1],Infinity,eps)) {
      position = Zero; 
    } else {
      position.real = -1e34;
      position.imag = complex_log(l.end[1], 0.0).imag;
    }
    return;
  } 
  if (same_point(l.end[0],Infinity,eps)) {
    distance = Complex(0.,PI); 
    if (same_point(l.end[1],Zero,eps)) {
      position = Complex(0.,PI); 
    } else {
      position.real = 1e34;
      position.imag = complex_log(l.end[1], 0.0).imag;
    }
    return;
  } 
  if (same_point(l.end[1],Zero,eps)) {
    distance = Complex(0.,PI);
    position.real = -1e34;
    position.imag = complex_log(l.end[0], 0.0).imag;
    return; 
  }
  if (same_point(l.end[1],Infinity,eps)) {
    distance = Zero; 
    position.real = 1e34;
    position.imag = complex_log(l.end[0], 0.0).imag;
    return; 
  }

  Complex f = complex_sqrt(l.end[0] * l.end[1]); 
  distance = complex_log((f + l.end[0])/(f - l.end[0]), 0.0);

  // Deal with special cases. 
  if (distance.real <= -eps) { // Real distance is negative. 
    f = -f; 
    distance = -distance; 
  } else if (distance.real < eps) { // Real distance is zero. 
    distance.real = 0.;
    if (distance.imag < 0.) { // Put imag part of distance in range [0,PI]. 
      f = -f;
      distance.imag = -distance.imag; 
    }
  }

  position = complex_log(f, 0.0);
  if (fabs(distance.imag + PI) < eps) distance.imag = PI; 
  if (fabs(position.imag + PI) < eps) position.imag = PI; 
}

#if 0
// This is discontinuous at +/-len.real/2, 
static void torus_normalize(Complex& z, Complex const& len)
{
  z      -= floor(z.real/len.real + 0.5) * len; 
  z.imag -= floor(z.imag/TWO_PI + 0.5) * TWO_PI;
}
#endif

// To get a well behaved comparison function on the torus
// we map continuously to a rectangle (rectangular
// torus) then apply folded_cmp below lexicographically
// on the components. 
// This function puts the real part of z into the range
// -len.real/2 +len.real/2, and the imaginary part into 
// the range -PI,PI.

static void torus_to_rect(Complex& z, Complex const& len)
{
  z      -= floor(z.real/len.real + 0.5) * len;
  z.imag -= len.imag/len.real * z.real; 
  z.imag -= floor(z.imag/TWO_PI + 0.5) * TWO_PI;
}

// Folded_cmp compares two real numbers, first
// by comparing their absolute values, and then, if
// those coincide, using the usual ordering. 
// This is useful for ordering points on a circle
// as represented on a symmetrically placed
// interval such as [-pi, pi]. Ordinarily the point
// at +/-pi could have either of these representatives
// and would therefore be indeterminately ordered
// with respect to any other point on the circle. 
// Using folded_cmp (and assuming a!=b on the circle)
// this will be the biggest point, 0 the smallest etc.

static int folded_cmp(double a, double b, double eps)
{
  double aa = fabs(a);
  double ab = fabs(b);
  if (fabs(aa-ab) > eps) 
    return (aa < ab) ? -1:1; 
  return (a < b) ? -1:1; 
}

// We need a comparison function for CGRays so that we
// can quickly test for duplicates. There is no particularly
// sensible way to order points on a torus so after testing
// for equality mod (g->length, 2*pi*i) we do something 
// fairly arbitrary. 

int cmp(CGRay const& a, CGRay const& b, double eps)
{
  if (a.g.index != b.g.index) 
    return (a.g.index < b.g.index) ? -1:1;

  Complex d = b.pos-a.pos;
  if (a.g.g) torus_to_rect(d, a.g.g->length); 
  // Allow for null geodesic pointer - don't bother to normalize.

  if (complex_small(d,eps)) return 0; 

  // Comparison has to be absolute, not relative; we compare
  // normalized versions of a.pos and b.pos, lexicographically. 

  Complex na = a.pos; if (a.g.g) torus_to_rect(na, a.g.g->length);
  Complex nb = b.pos; if (b.g.g) torus_to_rect(nb, b.g.g->length);

  if (fabs(d.real) > eps) 
    return folded_cmp(na.real, nb.real, eps);

  return folded_cmp(na.imag, nb.imag, eps);
}

// For output we just print the geodesic number and
// the position (which might not be normalized). A
// geodesic number of -3 is printed as an S, since this
// is used to represent an infinite spiraling geodesic. 

ostream& operator << (ostream& out, CGRay const& r)
{
  if (r.g.index==-3) out << " S"; else out << setw(2) << r.g.index;
  return out << ':' << r.pos;
}

// Constructors and modifiers all call set_low_end(). This means that
// although cmp(end[0],end[1]) might not always be -1, we always know
// whether it is this or 1 (if nonzero).

Ortholine::Ortholine(Complex const& d, CGRay const& a, CGRay const& b)
  : _distance(d)
{
  end[0] = a; 
  end[1] = b;
  set_low_end();
}

#if 0
Ortholine::Ortholine(Complex const& d, Complex const& a, Complex const& b, int na, int nb)
  : _distance(d) 
{ 
  end[0].pos = a; 
  end[1].pos = b; 
  end[0].gnum = na; 
  end[1].gnum = nb; 
  set_low_end();
}

Ortholine::Ortholine(const MoebiusTransformation& m, int na, int nb, FGWord const& w)
  : word(w)
{
  set_from_Moebius(m); 
  end[0].gnum = na; 
  end[1].gnum = nb;
  set_low_end(); 
}
#endif

Ortholine::Ortholine(const MoebiusTransformation& m, GSpec const& a, GSpec const& b)
{
  end[0].g = a; 
  end[1].g = b;
  set_from_Moebius(m); 
}

void Ortholine::set_to(Complex const& d, Complex const& a, Complex const& b)
{ 
  _distance=d; 
  end[0].pos = a; 
  end[1].pos = b; 
  set_low_end();
}

void Ortholine::set_from_Moebius(MoebiusTransformation const& m)
{
  Complex d0, d1; 
  MoebiusTransformation M = inverse(m); 
  get_ortholine_coordinates(line(m), d0, end[0].pos); 
  get_ortholine_coordinates(line(M), d1, end[1].pos);

  if (!complex_close(d0,d1,1e-7)) {
    cerr << "Problem with Ortholine initialization\n"; 
    cerr << "d0 = " << d0 << " d1 = " << d1 << endl;
    cerr << "l0 = " << line(m) << endl; 
    cerr << "l1 = " << line(M) << endl; 
    cerr << "m  = " << m << endl;
  }

  _distance = d0;
  set_low_end(); 
}

void Ortholine::sort_ends()
{
  if (low_end!=0) {
    swap(end[0],end[1]); 
    if (end[0].g.index==end[1].g.index) word.invert(); 
    low_end = 0; 
  }
}



#if 0
bool lxless(Complex const& a, Complex const& b, double eps)
{
  if (fabs(b.real-a.real) > eps) return a.real < b.real;
  if (fabs(b.imag-a.imag) > eps) return a.imag < b.imag;
  return false; // consider them equal.
}
#endif

// Ortholines are sorted first in order of orthodistance.
// Imaginary parts are treated symmetrically with respect
// to conjugation unless there is a tie, in which case
// one with negative imaginary part comes before its conjugate. 

bool operator < (const Ortholine& a, const Ortholine& b)
{
  if (fabs(a._distance.real - b._distance.real) > ort_end_eps)
    return a._distance.real < b._distance.real; 

  if (fabs(a._distance.imag - b._distance.imag) > ort_end_eps)
    return folded_cmp(a._distance.imag, b._distance.imag, ort_end_eps) < 0; 

  // Distances are equal so compare ends lexicographically.

  int sgn;
  sgn = cmp(a.end[a.low_end], b.end[b.low_end]);
  if (sgn) return sgn < 0; 

  sgn = cmp(a.end[1-a.low_end], b.end[1-b.low_end]);
  return sgn < 0; 
}

ostream& operator << (ostream& out, const Ortholine& o)
{
  out << o._distance << ' ';
  out << o.end[0] << ' ' << o.end[1] << ' ' << o.word;
  return out; 
}

typedef Tile* TilePtr;

line Ortholine::get_line(int i) const
{
  Complex a0,b0; 
  a0 = complex_exp(distance()); 
  b0 = -a0; 
  Complex a, b; 
  a = (a0 - One)/(a0 + One); 
  b = (b0 - One)/(b0 + One); 
  Complex zmul = complex_exp(position(i));
  return line(zmul * a, zmul * b); 
}

Ortholine::operator MoebiusTransformation() const
{
  line la = get_line(0); 
  line lb = get_line(1);

  Complex tr1[3], tr2[3]; 

  tr1[0] = la.end[0]; 
  tr1[1] = la.end[1];
  tr1[2] = Zero;
  tr2[0] = Zero;
  tr2[1] = Infinity;
  tr2[2] = lb.end[0];

  return MoebiusTransformation(tr1,tr2); 
}


#if 0
class interval_pair_accumulator {
public:
  virtual void operator() (interval const& a, interval const& b) =0; 
};
#endif

// In order to do sometimes lengthy ortholine computations incrementally
// we may search a tiling in order of radius for ortholines. This function
// looks for ortholines coming from the tiles between tile_rad_lo and 
// tile_rad_hi. It adds all found ortholines to ort_set. The process of 
// adding ortholines to an Ortholine_set removes duplicates. 

void get_ortholines(Tile* tiling, list<interval> const& ivls, double tile_rad_lo, double tile_rad_hi, Ortholine_set& ort_set, double cutoff, bool check_range)
{
  TilingIterator t; 
  Complex PiI(0.,PI);

  Ortholine ol; 
  list<interval>::const_iterator ivl, ivl2; 
  O31_matrix ivl_trans_inv; 
  MoebiusTransformation ort_tr; 

  for (ivl = ivls.begin(); ivl != ivls.end(); ivl++) {

    ivl_trans_inv = inverse(ivl->the_line().T()); 

    for (ivl2 = ivl; ivl2 != ivls.end(); ivl2++) {

      for (t = tiling; t; ++t) {

	// We compute the orthodistance from ivl to (tile_trans * ivl2),
	// = orthodist( ivl.trans * x_axis, tile_trans * ivl2.trans * x_axis ),
	// = orthodist( x_axis , inverse(ivl.trans) * tile_trans * ivl2.trans * x_axis )
	// = orthodist( inverse(ivl.trans) * tile_trans * ivl2.trans )

	ort_tr = ivl_trans_inv * t->g * ivl2->the_line().T();
	ol = Ortholine(ort_tr, ivl->geodesic(), ivl2->geodesic());

	// Ignore zero distance lines. 
	if (complex_small(ol.distance(),ort_end_eps)) continue; 
	if (complex_small(ol.distance()-PiI,ort_end_eps)) continue; 

	// Check if we`re inside the cutoff distance. 
	if (fabs(ol.distance().real) > cutoff) continue; 

	// Check that the ortholine meets the interval ivl (inside the
	// Dirichlet domain) so as to avoid getting too many duplicates.
	// We expand the interval slightly so as not to miss any hitting
	// an endpoint of the interval. 

	// Geodesics with numbers <= -3 are exempted from 
	// range checking: they are assumed to be infinite, not closed. 

	if (ol.geodesic_num(0) > -3 && ! ivl->in_range(ol.position(0).real)) continue; 
	if (ol.geodesic_num(1) > -3 && !ivl2->in_range(ol.position(1).real)) continue; 

	// Save the conjugacy but only if it means something. 
	if (ol.geodesic_num(0)==ol.geodesic_num(1))
	  ol.word = inverse(ivl->the_word()) * t->word * ivl2->the_word(); 

	ol.sort_ends();

	ort_set.insert(ol); 
      }
    }
  }
}


static void normalize_angle(Complex& z, Complex const& clen)
{
  z -= floor(z.real/clen.real + 0.5 - O31_line::epsilon) * clen; 
  z.imag -= floor(z.imag/TWO_PI + 0.5 - O31_line::epsilon) * TWO_PI;
}

void new_ortholines(const WEPolyhedron* poly, list<interval> const& li, 
		    vector<GeodesicWord> const& GL, 
		    double radius, Ortholine_set& OS)
{
  set<line_wd> lifts;

  get_lifts(poly, li, radius + outradius(li), lifts);

  Complex clen; 
  Ortholine O;
  O31_matrix T0;
  int inum; 
  list<interval>::const_iterator ivl; 
  set<line_wd>::const_iterator i; 
  for (ivl = li.begin(); ivl != li.end(); ++ivl) {
    T0 = inverse(ivl->the_line().T());
    inum = ivl->gnum(); 

    for (i = lifts.begin(); i != lifts.end(); ++i) {

      O = Ortholine(MoebiusTransformation(T0 * i->L().T()), inum, i->gnum);

      if (!ivl->in_range(O.position(0).real)) continue; 
      if (O.distance().real > radius) continue; 
      if (complex_small(O.distance(), O31_line::epsilon)) continue; 
#if 0
      normalize_angle(O.position(0), GL[inum].length);
      normalize_angle(O.position(1), GL[i->gnum()].length);
#endif
      O.sort_ends(); 

      OS.insert(O);
    }
  }
}

