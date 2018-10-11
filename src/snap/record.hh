#ifndef RECORD_H
#define RECORD_H
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

#include <cstdio>
#include <map>
#include <string>

using std::map; 
using std::string;

#if 0
class string_map 
: public map<string,string>
{
public:
  string_map() {}
};

class sstring
: public string
{
public:
  // explicit sstring (): string() {}
  sstring (): string() {}
//  sstring (const sstring& str): string(str) {}
  sstring (const string& str): string(str) {}
  sstring (const char* s): string(s) {}
};
#endif

typedef map<string,string> string_map; 

class record {
  string _tag;
  record *next; 
  record *sub; 
public:
  ~record() { delete next; delete sub; }
  record(const string& tag) : _tag(tag), next(0), sub(0) {}

  string tag() const { return _tag; }
  int has_subrecs() const { return sub!=0; }
  record* find_tag(const string& tag) const; 

  void add_tag(const string& tag);
};

FILE* find_record(const record& type, const string& value, 
		  const string& path, const string& file);

int load_record(const string& path, const string& file, 
		const record& type, const string& value, string_map& m);

bool load_next_record(const record& type, string_map& m, FILE* file_ptr, 
		   string& line);
#endif
