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

#include <iostream>
#include "env_A.hh"
#include "pariclass.hh"
#include <cstdio>

using std::sscanf;

#ifndef VERSION
#define VERSION "unknown"
#endif

class env : public env_A {
  static env M[10];

public:
  static void event_loop(); 

protected:
  virtual env_T* get_T(int i) { return &M[i]; }
  virtual env_D* get_D(int i) { return &M[i]; }
  virtual env_A* get_A(int i) { return &M[i]; }
};

env env::M[10];

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

#ifdef init
#undef init
#endif

int main(int argc, char *argv[])
{
  std::ios::sync_with_stdio(); 

  def_print.dec = 28; 

  env::init("snap"); 

  env::setup_menu(); 
  env::set_options(argc, argv); 

  // don't want to instantly break all tests!
  if (!env::in_batchmode()) 
    printf("snap %s\n", VERSION); 

  int j;
  for (j=1; j<argc; j++) {
    if (strcmp(argv[j], "-tg")==0) trace_garbage = true;
    if (strcmp(argv[j], "-s")==0) {
      j++; 
      long mb; 
      if (sscanf(argv[j],"%ld",&mb)==1) {
	if (mb<2) mb=2; 
	if (mb>10000) mb=10000; // this is probably too big anyway
	pari::resize_stack(1000000*mb);
      }
    }
  }

  env::setup_menu(); 
  env::event_loop(); 

  return 0; 
}
