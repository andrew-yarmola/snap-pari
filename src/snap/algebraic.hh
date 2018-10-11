#ifndef _algebraic_
#define _algebraic_
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


#include "snap.hh"
#include "field.hh"

class alg_snap : public snap {

  pari _alg_shapes; 

  field shapef; 
  vector<named_field> _cusp_fields; 

public:
  alg_snap() {}
  alg_snap(Triangulation *m) : snap(m) {}

  alg_snap* orientable_double_cover() const;

  /* algebraic.cc */ 
  virtual void update_precision(); 
  void compute_fields(int cc_poly = 0, int interactive = 1, unsigned int which = 0x1F); 

#if 0
  void save_fields(const char* filename) const; 
  int load_fields(const char *directory, const char* name, int which = 0xF); 
#endif

  field& tetfield() { return shapef; }
  vector<named_field>& cusp_fields() { return _cusp_fields; }
  const pari& alg_shapes() const { return _alg_shapes; }

  int verify() const;
  int verify_and_report() const;

  pari borel_regulator(field const& it, const field* e = 0) const; 
  pari cusp_info() const; 
  vector<int> cusp_commensurability_classes(); 

  virtual void print() const;
  void print_fields(int how=0) const; 

private:

  void clear_algebraic_info();
protected:

  virtual void clear();
  virtual void triangulation_changing();
  virtual void clear_shape_dependencies();
};

#endif
