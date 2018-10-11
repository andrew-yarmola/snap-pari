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
#include "mark.hh"
#include <cstdio>
#include <cctype>
#include <cstring>

using std::ostream;
using std::sscanf;

ostream& operator << (ostream& out, Mark m)
{
  int ml[2];
  m.get_ml(ml);

  int me = ml[0];
  int lo = ml[1];

  int face = m.num / 0x100;
  
  if (m.num==0) out << "-:-";
  else {
    out << face << ':';
    if (me==0&&lo==0) out << "0";
    else {
      if (me) {
	out << ((me>0) ? '+':'-');
	if (me<0) me = -me;
	if (me>1) out << me;
	out << 'M';
      }
      if (lo) {
	out << ((lo>0) ? '+':'-');
	if (lo<0) lo = -lo;
	if (lo>1) out << lo;
	out << 'L';
      }
    }
  }
  return out; 
}
  
bool Mark::set_to(const char* s)
{
  if (strcmp("-:-", s)==0) {
    num = 0; 
    return true;
  }

  int n; 
  char buf[50];
  if (sscanf(s, "%d:%s", &n, buf)!=2) return false; 
  int me=0, lo=0, t;
  int off=0;

  if (strcmp("0",buf)!=0) {

    while (buf[off]) {

      // get a sign
      if (buf[off]=='+') t=1;
      else if (buf[off]=='-') t=-1;
      else return false;

      // read optional digit
      if (isdigit(buf[off+1])) {
	t*= int(buf[off+1]-'0');
	++off;
      }

      // get M or L
      if (buf[off+1]=='M') me=t;
      else if (buf[off+1]=='L') lo=t;
      else return false; 
      off += 2;
    }
  }

  num = 0x100*n + 0x10*(me+8) + lo+8;
  return true; 
}

void Mark::get_ml(int H[]) const
{
  H[0] = ((num/0x10) & 0xf) - 8;
  H[1] =  (num & 0xf) - 8; 
}

int Mark::lensq() const
{
  int H[2];
  get_ml(H);
  return H[0]*H[0] + H[1]*H[1];
}

Complex Mark::holo(const Complex H[2]) const
{
  int ml[2];
  get_ml(ml);
  return double(ml[0])*H[0] + double(ml[1])*H[1];
}

ostream& operator << (ostream& out, da_spec const& da)
{
  return out << da.d << ' ' << da.a << ' ' << da.m; 
}

bool da_spec::read(FILE* fp)
{
  char sd[60], sa[60], sm[60];
  if (fscanf(fp, "%s %s %s", sd, sa, sm)!=3)
    return false; 

  if (!complex_from_string(sd, d)) return false; 
  if (!complex_from_string(sa, a)) return false; 
  return m.set_to(sm);
}

