#ifndef _env_T_
#define _env_T_
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

#include "menu.hh"
#include "n_vector.hh"

using std::string; 

class Triangulation; 
class GroupPresentation;
class i_triangulation;

class env_T {
protected:

  static string path; 
  static string prog; 
  static int snl; // shape normalization
  static bool batchmode;
  static bool echomode; 
  static bool simplify;
  static int tp_max; // tilt polytope, max num faces
  static text_menu mm;
  static int cm; // current manifold. 

private:

  static void show_menu();
  static void show_help();

public:

  static void init(string prg);
  static void set_options(int argc, char* argv[]); 
  static void setup_menu();
  static bool in_batchmode() { return batchmode; }

protected:

  i_triangulation* T;
  GroupPresentation* G;
  string gp_name; 
  int closed_num; 

  env_T() : T(0), G(0), closed_num(0) {}
  virtual ~env_T(); 

  int choose_save_manifold();
  int choose_manifold(string const& prefix = "");
  void set_manifold_from_surgery_description(string const& s);


  bool get_group(bool required=true); 
  void clear_group(); 
  virtual bool read_group(); 

  int get_event(); 
  void validate_event(int& what); 
  void process_event(int what); 
  int num_cusps() const; 

  virtual void print_settings() const; 

  virtual env_T* get_T(int i) =0;

  i_triangulation* get_TT(int i) { return get_T(i)->T; }

  virtual i_triangulation* new_manifold(Triangulation* tri) const; 

  virtual void T_shape_changed();
  virtual void T_changed();
  virtual void set_G(GroupPresentation* g); 
public:

  void set_T(i_triangulation* m);
  bool has_manifold() const { return T!=0||G!=0; }
  string name(int style=0) const; 

  GroupPresentation* group() const { return G; }
};

enum {
  Nothing=0, 
  Help, 
  Quit,
};

extern const char* solution_types[7];

struct SymmetryGroup;
void print_symmetry_group(SymmetryGroup* symm_gp);
string choose_save_file(string def_dir, string def_name);
n_vector get_cusp_sizes(i_triangulation const& m, bool always=false);

#endif
