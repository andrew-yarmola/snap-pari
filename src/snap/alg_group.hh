#ifndef _alg_group_
#define _alg_group_
/*
** Copyright (C) 2003,2005 Oliver A. Goodman <oag@ms.unimelb.edu.au>
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

#include "field.hh"
#include "snappea/SnapPea.h"

class alg_group {
  GroupPresentation* G;

  pari exact_gens; 

  pari tracegens[2], itgens[2]; // 0=numeric, 1=exact

  field trace, it, groupf; 
  int _integer_traces; 

  pari hilbert[2]; 

public:

  alg_group() : G(0), _integer_traces(-1) {}

  GroupPresentation* group() const { return G; }
  void set_group(GroupPresentation* g); 
  void clear(); 
  void clear_gf() { groupf.clear(); }
  void set_tracefield_gens();
  void update_precision(); 
  void compute_fields(int cc_poly = 0, int interactive = 1, unsigned int which = 0x1F); 

  void print() const; 
  void print_fields(int how=0) const; 
  void print_group_generators() const;
  int integer_traces() const;

  // field& itfield() { return it; }
  const field& itfield() const { return it; }
  const field& tracefield() const { return trace; }
  const field& groupfield() const { return groupf; } 
  const pari& alg_traces() const { return tracegens[1]; }
  pari alg_group_element(const FGWord& word) const; 
  pari exact_trace(FGWord const& w) const;
  const pari& get_exact_gens() const { return exact_gens; }
  int verify() const; 

  bool compute_hilbert_symbol(); 
  const pari* hilbert_symbol() const { return hilbert; }
  void print_hilbert_symbol(int show) const; 

  pari real_ramification() const;
  pari finite_ramification(int max_time=0) const;
  pari ramification(int max_time=0) const;

  int arithmetic() const;

  bool imquad_integer_itgens() const;

private:

  void set_exact_gens(); 
  int has_integer_traces() const; 
};

pari ramification(const pari hilbert[2], pari itf, int int_traces, int report, int max_time); 

#endif
