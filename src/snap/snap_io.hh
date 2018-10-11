#ifndef _snap_io_hh_
#define _snap_io_hh_

/*
** Copyright (C) 2004 Oliver A. Goodman <oag@ms.unimelb.edu.au>
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
#include <cstdio>

using std::string;

class Triangulation;

Triangulation* read_manifold_file(string const& path, string const& name);
Triangulation* read_closed(string const& path, int index);
Triangulation* read_surgery_description(string const& path, string const& sd, int mod);

string lookup_link_dtcode(string const& path, int ncross, char alt, int index);
string find_file_name(const string& path, const string& name, const string& opt = "");
FILE* open_data_file(string const& path, string const& name, string const& prefix="");
char* get_line(FILE* fp, int& index, int skip=0);

string find_unused_name(string dir, string base, string ext);
bool file_exists(string const& file); 

#endif
