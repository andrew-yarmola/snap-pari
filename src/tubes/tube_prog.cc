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
#include "env_U.hh"
#include <iomanip>

using std::ios;

class env : public env_U {
protected:
  static env M[10]; 

public:
  static void event_loop(); 

protected:
  env() {}
  virtual ~env() {}

  virtual env_T* get_T(int i) { return &M[i]; }
  virtual env_D* get_D(int i) { return &M[i]; }
  virtual env_T* get_env(int i);
};

void env::event_loop()
{
  int what;

  while (1) {

    what = M[cm].get_event();

    if (what==Quit) break; 

    M[cm].validate_event(what);

    M[cm].process_event(what); 
  }
}

int main(int argc, char* argv[])
{
  cout.setf(ios::fixed); // print 1e-9 as 0.00000

  env::init("tube");

  env::set_options(argc, argv); 
  env::setup_menu(); 
  env::event_loop(); 

  return 0; 
}

env_T* env::get_env(int i)
{ 
  return &M[i]; 
}

env env::M[10];

