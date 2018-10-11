#ifndef _menu_
#define _menu_
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


#include <string>
#include <vector>
#include <iostream>

#include "warn.hh"

using std::string;

struct text_menu_item {
  string word; 
  text_menu_item* next; 
  text_menu_item* submenu; 
  int token; 

  text_menu_item() : next(0), submenu(0), token(0) {}
  text_menu_item(const string& w);
  ~text_menu_item();
};

class text_menu {
  text_menu_item root; 
public:
  text_menu() {}

  void add_item(string in, int token);

  int read_input(const string& prompt, char** line=0) const; 
  string expand_input(const string& prompt, int& tok) const; 
  void get_menu(std::vector<string>& items) const;
  void print_menu(std::ostream& out = std::cout) const; 
  string list_options(const text_menu_item* menu_ptr) const; 
};

std::ostream& operator << (std::ostream& out, text_menu const& m);
int get_input(string& s, const string& pr, int words = 1, char** line=0);
int match(const string& in, const string& tok, int pre_len);
int ask(string pr, int def);
void print_history(std::ostream& out = std::cout);
void die(const string& msg); 
void print_help(const string& help_topic, const string& name, const string& path);
int count_words(const string& s);
void get_name_path(string const& file, string& name, string& path);
void get_base_ext(string const& file, string& base, string& ext);

void get_documented_topics(const string& name, const string& path, std::vector<string>& items);

#endif
