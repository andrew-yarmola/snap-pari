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
#include "record.hh"
#include "helpers.hh"

using std::vector; 

// #include "say.hh"
// static say xx("record.cc"); 

void record::add_tag(const string& tag)
{
  record* n = sub; 
  if (!n) {
    sub = new record(tag); 
    return; 
  }

  while (n->next) n = n->next; /* go to the end of the list of records */ 
  n->next = new record(tag); 
}

record* record::find_tag(const string& t) const
{
  record* n = sub; 

  while (n) {
    if (n->tag()==t) return n; 
    n = n->next; 
  }
  return 0; 
}


static bool locate(FILE *f, const string& search)
{
  static char input[1001];
  int n = search.length(); 
  while (fgets(input, 1000, f)) {
	if (input[n]=='\n') input[n] = '\0'; 
#ifdef OLD_GNU_COMPARE
    if (search.compare(input,0,n) == 0) return true; 
#else
    if (search.compare(0,n,input) == 0) return true; 
#endif
  }
  fclose(f); 
  return false; 
}

int find_record(const record& type, const string& value, FILE* file)
{

/* look for a record with type "type" and value "value" in file "file". 
   on entry, file points at the start of the file. on exit it points at 
   the start of the line following the type:value pair found. return
   non-zero if the record is found, zero if not. If the record is not 
   found the file is closed and zero is returned. */ 

  string search = type.tag(); 
  search += ":" + value; 

  return locate(file, search); 
}

FILE* find_record(const record& type, const string& value, 
		  const string& path, const string& file)
{
  
/* path is assumed to be a blank delimited list of files in which to 
   search for record. */ 

  vector<string> dirs; 
  int n = split(path, dirs, ' ');
  string full_name; 

  int found = 0, i = 0;
  FILE *f; 
  for (;!found && i < n; i++) {
    full_name = dirs[i] + "/" + file; 
    f = fopen(full_name.c_str(), "r"); 
    if (!f) continue;
    found = find_record(type, value, f);
  }

  if (found) return f;
  return 0; 
}

/* when load_record encounters a tag which does not belong to the 
   record "type" it will have read a line intended for the calling function.
   Hence it returns the line so that the calling function can use it. */ 

string load_record(const record& type, string_map& m, FILE* file_ptr, 
		   string prefix = "")
{

  /* for each tag occurring in record type and present in the file,
     set m["tag"] to the value from the file. If a record has a
     subrecord of its own then we set m["tag:subtag"] */

  string line, tag; 
  record* rec; 

  static char buf[20001];
  while (line.length() || fgets(buf, 20000, file_ptr)) {
    line = buf; 
    tag = line.substr(0,line.find(":")); 
    rec = type.find_tag(tag); 
    if (!rec) return line; 
    if (rec->has_subrecs()) {
      line = load_record(*rec, m, file_ptr, prefix + tag + ":");
    } else {
      m[prefix + tag] = line.substr(line.find(":")+1); 
      line = ""; 
    }
  }

  return ""; /* if we get here we reached end of file so there is 
		no next line to return */ 
}

int load_record(const string& path, const string& file, 
		const record& type, const string& value, string_map& m)
{
  FILE *fp = find_record(type, value, path, file); 

  if (!fp) return 0; 

  load_record(type, m, fp);
  fclose(fp); 

  return 1; 
}

bool load_next_record(const record& type, string_map& m, FILE* file_ptr, 
		   string& line)
{
  string tag; 

  static char buf[20001];
  if (!line.length()) {
    if (!fgets(buf, 20000, file_ptr)) return false; 
    line = buf; 
  }
  tag = line.substr(0,line.find(":")); 
  if (type.tag()!=tag) return false; 
  m[tag] = line.substr(line.find(":")+1); 
  line = load_record(type, m, file_ptr);
  return true; 
}


