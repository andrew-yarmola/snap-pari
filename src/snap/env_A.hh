#ifndef _env_A_
#define _env_A_

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

#include "env_U.hh"

class alg_snap;
class alg_group; 

class env_A : public env_U {
  static double acc_vertex_eps; 
  static int cc_poly; 
  static int print_fields_long; 
  static int show_field_comp; 
  static int hilb_print; 
  static bool cs_diagnostic; 
  static text_menu fm; // menu of fields
  static int sally_verbosity;
  static bool save_sally_results;

public:

  env_A() : AG(0) {}

  static void setup_menu(); 

protected:

  void validate_event(int& what);
  void process_event(int what);
  virtual void print_settings() const;
  virtual bool read_group(); 

  void normalize_group(vector<FGWord> const& wl, bool report=false); 

  alg_snap* get_AS() { return (alg_snap*)T; }
  const alg_snap* get_AS() const { return (const alg_snap*)T; }
  virtual env_A* get_A(int i) =0;
  virtual void set_G(GroupPresentation* g); 

  virtual i_triangulation* new_manifold(Triangulation* tri) const;

  void eigenvalue_fields(bool canonical, int report, bool find_unit, double max_len) const;
  void update_precision(); 
  void sync_AG();
  void sally(int gn, int report);  

  alg_group *AG;
};

#endif
