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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "snap_io.hh"
#include "snappea/unix_io.h"
#include "helpers.hh"
#include <cstdio>

using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::sscanf;

#define BSIZE 512

Triangulation* read_surgery_description(string const& path, string const& s, int mod)
{
  char census_char; 
  double me=0., lo=0.;
  int N, mn;
  int o=0;
  while (s[o] && s[o]==' ') o++; 
  int narg = sscanf(s.c_str()+o, "%c%d( %lf , %lf )", 
		    &census_char, &mn, &me, &lo);
  if (narg!=2 && narg!=4) return 0; 
  
  /* get which census */ 
  if (census_char=='m') N = 5; 
  else if (census_char=='s') N = 6;
  else if (census_char=='v') N = 7;
  else return 0; 
  
  /* Read the manifold. */ 
  Triangulation* tri = read_census_manifold(path.c_str(), N, mn);
  if (!tri) return 0;

  if (mod==1) randomize_triangulation(tri); 
  else if (mod==2) proto_canonize(tri); 

  set_cusp_info(tri, 0, (me==0. && lo==0.), me, lo);
  do_Dehn_filling(tri);

  return tri;
}

// if line number is beyond length of file
// sets index to the number of the last line. 

char* get_line(FILE* fp, int& index, int skip)
{
  static char buf[BSIZE];
  int i=0; 
  while (i < skip && fgets(buf, BSIZE, fp)) i++; 
  i = 1; 
  while (i < index && fgets(buf, BSIZE, fp)) i++; 
  if (!fgets(buf, BSIZE, fp) || i!=index) {
    index = i-1; 
    fclose(fp); 
    return 0;
  }
  fclose(fp); 
  return buf;
}

static int can_rqd[] = { 
  386,
  1857,
  1929,
  1960,
  2077,
  2174,
  3223,
  3358,
  4066,
  4405,
  5409,
  5794,
  6537,
  6636,
  7090,
  7313,
  7350,
  7664,
  7962,
  8744
};

Triangulation* read_closed(string const& path, int index)
{
  FILE* fp = open_data_file(path, "ClosedManifolds"); 
  if (!fp) return 0; 

  char* line = get_line(fp, index, 31); // data starts after line 31
  if (!line) {
    cout << "index must be in range 1-11031\n";
    return 0;
  }

  int rand = 0; 
  int j; 
  for (j=0; j<20; j++) 
    if (can_rqd[j]==index) rand = 1; 

  // surgery description is in cols 72-83.
  line[84] = '\0';
  return read_surgery_description(path, line+72, rand);
}

string find_file_name(const string& path, const string& name, const string& opt)
{
  string path_name; 
  std::vector<string> dirs;

  int i = split(path, dirs);
  int j;
  FILE* fp; 
  for (j=0; j<i; j++) {
    path_name = dirs[j] + "/" + name;
    fp = fopen(path_name.c_str(), "r");
    if (fp) {
      fclose(fp); // We are not actually interested in the contents here. 
      return path_name; 
    }
    if (!opt.size()) continue; 
    path_name = dirs[j] + "/" + opt + "/" + name; 
    fp = fopen(path_name.c_str(), "r");
    if (fp) {
      fclose(fp); // We are not actually interested in the contents here. 
      return path_name; 
    }
  }
  return ""; 
}

Triangulation* read_manifold_file(string const& path, string const& name)
{
  string file = find_file_name(path, name, "manifolds"); 
  if (!file.size()) {
    cout << "File " << name << " not found!\n"; 
    return 0; 
  }

  FILE *fp = fopen(file.c_str(), "r"); 
  Triangulation* tri = read_manifold_file(fp); 
  fclose(fp); 

  if (!tri) {
    cout << "Unable to read: " << file << endl;
  }
  return tri; 
}

FILE* open_data_file(string const& path, string const& name, string const& prefix)
{
  string pathname = find_file_name(path, prefix+name); 

  if (!pathname.length()) { 
    cerr << "couldn't locate " << prefix << name << endl; 
    return 0; 
  }

  // use shorter name in error messages now. 

  FILE* fp = fopen(pathname.c_str(), "r"); 
  if (!fp) {
    cerr << "problem opening " << name << endl; 
    return 0; 
  }
  return fp; 
}

string lookup_link_dtcode(string const& path, int ncross, char alt, int index)
{
  string code; 

  if (alt!='a' && alt!='n') {
    cout << "invalid link table specifier\n"; 
    return code; 
  }

  if (ncross < 4 || (alt == 'n' && ncross < 6)) {
    cout << "invalid number of crossings specified\n"; 
    return code; 
  }

  if (index < 1) {
    cout << "index must be greater than 0\n"; 
    return code; 
  }

  char lfile[BSIZE]; 
  sprintf(lfile, "hyperbolic_data_%02d%c", ncross, alt);

  FILE* fp = open_data_file(path, lfile, "link_data/");
  if (!fp) return code; 

  char* line = get_line(fp, index);
  if (!line) {
    cerr << "index must be in range 1-" << index << endl;  
    return code;
  }

  double vol;
  if (sscanf(line, "%*s %lf", &vol)!=1) 
    cout << "Warning: manifold is non-hyperbolic!\n";

  // write a terminating '\0' after the code.
  int i = 0; 
  while (line[i] && line[i]!=' ') i++; 
  if (line[i]==' ') line[i]='\0';
  code = string(line); 

  return code; 
}

bool file_exists(string const& file)
{
  FILE* fp = fopen(file.c_str(), "r");
  if (fp) {
    fclose(fp); 
    return true; 
  }
  return false; 
}

static string stringof(int i)
{
  char buf[20];
  sprintf(buf, "%d", i); 
  return string(buf); 
}

string find_unused_name(string dir, string base, string ext)
{
  int i=1, lim=1000000; 
  string name; 
  dir += "/"; 
  while (i < lim) {
    name = base + stringof(i) + ext; 
    if (!file_exists(dir+name)) break;
    i *= 2; 
  }
  if (i>lim) return ""; 
  int j; 
  for (j=i/2+1; j<= i; j++) {
    name = base + stringof(i) + ext; 
    if (!file_exists(dir+name)) 
      return name; 
  }
  return ""; // should never get here. 
}

