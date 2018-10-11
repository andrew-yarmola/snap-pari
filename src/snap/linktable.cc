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

#include "file_array.hh"
#include "snappea/SnapPea.h"
#include "snap_io.hh"

using std::string;
using std::cout;
using std::endl;

struct link_entry {
  double volume;
  char dowker[26];

  static double eps; 

  link_entry(double vol=0.) : volume(vol) { dowker[0] = '\0'; }

  void print(FILE* fp = stdout) const;

  friend bool operator < (link_entry const& a, link_entry const& b)
  { return a.volume < (b.volume-eps); }
};

double link_entry::eps = 1e-9; 

string find_in_linktable(string const& path, Triangulation* m)
{
  string name = find_file_name(path, "links_by_volume"); 
  if (!name.length()) {
    cout << "couldn't find file links_by_volume\n";
  }

  file_array<link_entry> table(name.c_str()); 

  int i, n, first; 
  n = table.find(volume(m,0), first);

  Boolean isometric; 
  Triangulation* m1; 

  string code; 
  for (i=first; i<first+n; i++) {
    code = table[i].dowker; 
    m1 = DT2Triangulation(code);
    if (compute_isometries(m,m1,&isometric,0,0)!=func_OK) {
      cout << "trouble checking isometry for code " << code << endl;
      continue; 
    }
    if (isometric) break; 
  }

  if (i<first+n) {
    return code; 
  }
  return string(); 
}

void link_entry::print(FILE* fp) const
{ 
  fprintf(fp, "%s \t%lf\n", dowker, volume); 
}
