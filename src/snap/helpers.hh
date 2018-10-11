#ifndef _helpers_
#define _helpers_
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

#include <vector>
#include <string>
#include <iostream>

using std::ios;
using std::string;
using std::vector;
using std::ostream;

void print_vector(ostream& out, vector<int> const& v);
bool scan_vector(const char* str, vector<int>& v);
vector<int> read_numbers(const char* c);
string reverse(string s);
int split(const string& in, string subs[], int nmax, char delim = ' ');
int split(const string& in, vector<string>& words, char delim = ' '); 

#ifdef __DECCXX // DEC c++ doesn't define a type for this enumeration. 
typedef long ios_fmtflags;
#else
typedef ios::fmtflags ios_fmtflags;
#endif

#endif
