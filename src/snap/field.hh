#ifndef _field_
#define _field_
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


#include <string>
#include <vector>
#include "pariwrap.hh"

using std::string;
using std::vector;

int is_integral(const pari& elt);
pari canonical_poly(const pari& poly, pari& newroot_on_old, int report);
pari reduced_poly(const pari& poly, pari& newroot_on_old, int report);
pari poly_denominator(const pari& poly);

class field {
  pari nf; 

  int root_num; 
  vector<int> root_index; 

  bool _is_canonical; 

  pari numeric_ibasis; 
  pari ibasis_lllmat;

public:
  static int maxdeg; 
  static double fudge; 

  field();
  field(const field& f); 
  field(const pari& min_poly, int root_num, bool is_can); 
  field(const pari& min_poly_and_rn); 

  ~field();

  field& operator = (const field& f); 

  int is_set() const { return root_num; } 
  int is_rational() const { return length(nf[0])==2; }
  int is_real() const;
  bool is_canonical() const { return _is_canonical; }
  pari min_poly() const { return nf[0]; }
  int root_number(const pari& r) const; 
  int root_number() const { return root_num; } 
  pari root(int n) const;
  pari root() const { return root(root_num); }  
  pari p_nf() const { return nf; }
  int degree() const { return lgeff(min_poly()) - 1; }
  bool amphicheiral() const;
  pari spec() const; 

  int contains(pari elt) const;
  int contains(pari elt, pari& exact_elt, int report = 0) const;

  pari numeric_value(pari element, int emb = 0) const; 
  pari exact_value(const pari& x) const;
  void set_exact_value(pari& x) const; 

  void print(int how) const; 

  /* Modifiers. */ 
  void clear();
  void update_precision(); 

  int set_field(pari const& poly, pari const& root, int report);
  int set_field(pari const& poly, pari const& root, pari& exact_root, int cp=0, int report=0);

  // Returns 2 if field changed, 1 if unchanged, 0 if failure.
  int extend(pari const& elt, pari& exact_elt, int cp = 1, int report = 0, int gnum = 0); 

  int generated_by(const pari& elts, pari& exact_elts, int cp = 1, int report = 0); 
  int generated_by(const pari& elts, pari& exact_elts, int cp, const field *f, int relation = 0, int report = 0); 

  int canonize(pari& newroot_on_old, int report = 0); 

  int find_canonical_root(pari& new_root_on_old, int report);

private:

  void set_root_indices(); 
};

extern field rationals; /* these are just names, they aren't set to anything */ 
extern field complex_infinity; 


class named_field : public field {
  string _name; 
public:

  named_field& operator = (const field& f)
    { field::operator = (f); return *this; } 

  void set_name(const string& name) { _name = name; }
  string const& name() const { return _name; }
  void print(int n) const; 
};

// Objects of this class are comparison function objects. 
// We assume a fixed list of roots, supplied in the constructor: 
// eg. root_index_less rless(root_list); 
// Then rless(i,j) is true iff root_less(root_list[i],root_list[j]). 
//
// Using this, rather than sorting a C++ vector of pari (complex type) roots using 
// root_less, we can sort a vector of integer indices into a pari vector of roots.
//
bool root_less(const pari& a, const pari& b);

class root_index_less {
  pari rl; 
public: 
  root_index_less(const pari& root_list) : rl(root_list) {}

  bool operator() (int a, int b) const { return root_less(rl[a],rl[b]); } 
};

#endif
