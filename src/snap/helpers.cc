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
#include "helpers.hh"
#include <cstdio>

using std::sscanf;
using std::swap;


vector<int> read_numbers(const char* c)
{
  int n; 
  vector<int> numbers;

  while (*c && isspace(*c)) c++; // eat any leading space.
  while (*c) {
    if (!sscanf(c, "%d", &n)) // read a number. 
      return numbers; // something unreadable as a number.. return as much as we got. 
    numbers.push_back(n); // save the number. 
    while (*c && !isspace(*c)) c++; // eat to end of the number. 
    while (*c && isspace(*c)) c++; // eat space after it. 
  }
  return numbers;
}

string reverse(string s)
{
  int i, n = s.length(); 
  for (i=0; i<n/2; i++)
    swap(s[i],s[n-i-1]); 
  return s; 
}

void print_vector(ostream& out, vector<int> const& v)
{
  int i, n = v.size(); 
  out << '[';
  for (i=0; i<n; i++) {
    out << v[i]; 
    if (i<n-1) out << ','; 
  }
  out << ']'; 
}

bool scan_vector(const char* str, vector<int>& v)
{
  int i=0, n, x;
  if (str[0]!='[') return false; 
  str++;
  while (true) {
    if (!*str) return false; // make sure we didn't get to the end prematurely
    if (sscanf(str, "%d %n", &x, &n)!=1) return false;  // read a number
    str+=n; 
    if (i>=v.size()) v.push_back(x); else v[i]=x; 
    i++; 
    if (*str==',') str++; 
    else if (*str==']') break; 
  }
  if (i<v.size()) v.resize(i); 
  return true; 
}

int split(const string& in, string subs[], int nmax, char delim)
{
  int i=0,j,n=0;

  /* Split in into separate words. */ 
  while (n<nmax) {
    while (i < in.length() && in[i]==delim) ++i; /* Skip spaces. */
    if (i==in.length()) break; 
    j = i; 
    while (i < in.length() && in[i]!=delim) ++i; /* Go to the end of the word. */
    subs[n] = in.substr(j,i-j); 
    ++n; 
  }
  return n; 
}

int split(const string& in, vector<string>& words, char delim)
{
  int i=0,j,n=0;
  words.resize(0);

  /* Split in into separate words. */ 
  while (true) {
    while (i < in.length() && in[i]==delim) ++i; /* Skip spaces. */
    if (i==in.length()) break; 
    j = i; 
    while (i < in.length() && in[i]!=delim) ++i; /* Go to the end of the word. */
    words.push_back(in.substr(j,i-j));
    ++n; 
  }
  return n; 
}


