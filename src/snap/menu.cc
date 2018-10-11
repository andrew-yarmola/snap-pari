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
#include "menu.hh"
#include <stdio.h>
#include <iostream>
#include <ctype.h>
#include <list>
#include <algorithm>
#include "helpers.hh"

extern "C" {
#include <readline/readline.h>
#include <readline/history.h>
}

using namespace std;

int get_input(string& s, const string& pr, int words, char** line)
{
  static char *buf =0;
  static int n = 0;
  s = ""; 

  if (!buf) { buf=(char*)malloc(1); *buf='\0'; }

  if (line) *line=buf; 
  while (1) {
    if (!buf[n]) { /* nothing waiting, give prompt */
      free(buf); 
      buf = readline(pr.c_str()); /* prompt and await input */ 
      if (line) *line=buf; 
      if (!buf) exit(0); /* ctrl-d ends program */ 
      if (*buf) add_history(buf); /* save non-empty lines */ 
      n = 0;
      if (!buf[n]) return 0; /* empty line */ 
    }
    while (buf[n] && isspace(buf[n])) n++; /* skip leading space */
    if (!buf[n]) continue; /* only space left, get more */ 

    if (words==-1) { /* get a whole line */ 
      while (buf[n] && buf[n]!='\n') 
	s += buf[n++];
      break;
    }

    if (!words) break; /* read all we wanted ? */
    while (buf[n] && !isspace(buf[n])) 
      s += buf[n++];
    words--; 
    if (words > 0) s += ' ';
    if (!words) break; 
  }
  return 1;
}

int match(const string& in, const string& tok, int pre_len)
{
  int n = in.length(); 
  if (n < pre_len) return 0; 
#ifdef OLD_GNU_COMPARE
  return tok.compare(in, 0, in.length())==0;
#else
  return tok.compare(0, in.length(), in)==0;
#endif
}

static void downcase(string& s)
{
  int i, n = s.length(); 
  for (i=0; i<n; ++i) s[i] = tolower(s[i]); 
}

int ask(string pr, int def)
{
  string s;
  pr += (def) ? " (yes) " : " (no) ";

  if (!get_input(s, pr)) return def; 
  downcase(s);
  return match(s, "yes", 1);
}

void die(const string& msg)
{
  cerr << msg;
  exit(0);
}

static int split(const string& in, string subs[], int nmax, string& rest, char delim = ' ')
{
  int i=0,j,n=0;

  /* Split in into separate words. */ 
  while (true) {
    while (i < in.length() && in[i]==delim) ++i; /* Skip spaces. */
    if (i==in.length() || n==nmax) break; 
    j = i; 
    while (i < in.length() && in[i]!=delim) ++i; /* Go to the end of the word. */
    subs[n] = in.substr(j,i-j); 
    ++n; 
  }
  rest = in.substr(i); 
  return n; 
}

static void remove_leading_whitespace(string& s)
{
  int i=0;
  while (i<s.length() && isspace(s[i])) ++i; 
  s = s.substr(i); 
}

static FILE* find_file(const string& path, const string& name, char *mode)
{
  string path_name; 
  vector<string> dirs;

  int i = split(path, dirs);
  int j;
  FILE* fp; 
  for (j=0; j<i; j++) {
    path_name = dirs[j] + "/" + name;
    fp = fopen(path_name.c_str(), mode);
    if (fp) return fp; 
  }
  return NULL; 
}

void print_history(ostream& out)
{
  HIST_ENTRY *the_entry;
  int i;
       
  for (i = 0; the_entry = history_get(i); i++)
    out << the_entry->line << endl; 
}

#define BUFSIZE 1024

void print_help(const string& help_topic, const string& name, const string& path)
{
  char lbuf[BUFSIZE];

  FILE* help_file = find_file(path, name, "r");
  if (!help_file) {
    cout << "Couldn't locate the help file " << name << '\n';
    cout << "in: " << path << '\n'; 
    return; 
  }
  bool topic_found = false; 
  int nchars = help_topic.length(); 
  list<string> headlines, empty; 
  
  while (!topic_found) {
    // Find a line starting with a '.'
    while (fgets(lbuf, BUFSIZE, help_file))
      if (lbuf[0] == '.') break;
    if (lbuf[0] != '.') break; // No such line found. 

    // Save first line starting with a '.'
    headlines.push_back(string(lbuf)); 
#ifdef OLD_GNU_COMPARE
    topic_found = (help_topic.compare(lbuf+2, 0, nchars) == 0);
#else
    lbuf[nchars+2] = '\0'; // Truncate input line
    topic_found = (help_topic.compare(0, nchars, lbuf+2) == 0);
#endif

    // Save all lines starting with a '.' and one more. 
    while (fgets(lbuf, BUFSIZE, help_file)) {
      headlines.push_back(string(lbuf)); 
#ifdef OLD_GNU_COMPARE
      if (help_topic.compare(lbuf+2, 0, nchars) == 0) topic_found = true;
#else
      lbuf[nchars+2] = '\0'; // Truncate input line
      if (help_topic.compare(0, nchars, lbuf+2) == 0) topic_found = true;
#endif
      if (lbuf[0] != '.') break;
    }

    if (topic_found) break; 
    headlines = empty; 
  }

  if (!topic_found) {
    cout << "Sorry, \"" << help_topic << "\" is undocumented.\n"; 
    fclose(help_file); 
    return; 
  }

  // Print headlines (which include the topic name), and all 
  // lines up to but not including the next line starting with a '.'
  cout << '\n';
  list<string>::const_iterator it; 
  for (it = headlines.begin(); it!=headlines.end(); it++)
    cout << *it;
  while (fgets(lbuf, BUFSIZE, help_file)) {
    if (lbuf[0]=='.') break;
    cout << lbuf;
  }

  fclose(help_file); 
}

void get_documented_topics(const string& name, const string& path, vector<string>& doc)
{
  char lbuf[BUFSIZE];

  FILE* help_file = find_file(path, name, "r");
  if (!help_file) {
    cout << "Couldn't locate the help file " << name << '\n';
    cout << "in: " << path << '\n'; 
    return; 
  }
  int l; 
  string s; 
  while (true) {

    // Find line starting with . 
    while (fgets(lbuf, BUFSIZE, help_file))
      if (lbuf[0] == '.') break;
    if (lbuf[0] != '.') break; // No such line found. 
    s = string(lbuf+2);
    l = s.length(); 
    doc.push_back(s.substr(0,l-1)); 
  }

  sort(doc.begin(), doc.end()); 
}

int count_words(const string& s)
{
  int i = 0, n = 0; 
  while (i<s.length()) {
    while (s[i] && isspace(s[i])) i++; 
    if (s[i]) n++; else return n;
    while (s[i] && !isspace(s[i])) i++;
  }
  return n;
}

static void print_submenu(ostream& out, const text_menu_item* menu_ptr, int indent)
{
  if (!menu_ptr) return; 
  int i; 
  while (1) {
    for (i = 0; i < indent; i++) 
      out << " "; 
    out << menu_ptr->word << endl; 
    print_submenu(out, menu_ptr->submenu, indent+2);
    if (!menu_ptr->next) return; 
    menu_ptr=menu_ptr->next; 
  }
}

static void get_submenu(vector<string>& items, const text_menu_item* menu_ptr, const string& prefix)
{
  string this_item; 
  for (;menu_ptr; menu_ptr=menu_ptr->next) {
    this_item = prefix + menu_ptr->word; 

    if (menu_ptr->submenu) 
      get_submenu(items, menu_ptr->submenu, this_item + ' ');
    else 
      items.push_back(this_item); 
  }
}

void text_menu::print_menu(ostream& out) const
{
  print_submenu(out, root.submenu, 0); 
}

void text_menu::get_menu(vector<string>& items) const
{
  get_submenu(items, root.submenu, ""); 
  sort(items.begin(), items.end());
}

ostream& operator << (ostream& out, text_menu const& m)
{
  m.print_menu(out); return out; 
}

int get_menu_item(const string& w, const text_menu_item*& menu_ptr)
{
  const text_menu_item* saved = menu_ptr; 

  if (!menu_ptr->submenu) return -1; 
  menu_ptr = menu_ptr->submenu; 

  /* look for the word on this level */ 
#ifdef OLD_GNU_COMPARE
  while ((menu_ptr->word.compare(w,0,w.length())!=0) && menu_ptr->next) 
#else
  while ((menu_ptr->word.compare(0,w.length(),w)!=0) && menu_ptr->next) 
#endif
    menu_ptr=menu_ptr->next;

#ifdef OLD_GNU_COMPARE
  if (menu_ptr->word.compare(w,0,w.length())!=0) {
#else
  if (menu_ptr->word.compare(0,w.length(),w)!=0) {
#endif
    menu_ptr = saved; 
    return -1; 
  }

  return menu_ptr->token; 
}

string first_word(const string& s)
{
  string word[1]; 
  if (split(s,word,1)) return word[0]; 
  return string(""); 
}

// Return values for tok are as follows:
// tok > 0 : a complete menu item was read. 
// tok== 0 : an incomplete menu item was read. 
// tok==-1 : an invalid keyword was read. 

int text_menu::read_input(const string& pr, char** line) const
{
  string w, sofar, pr2; 
  const text_menu_item *menu_ptr=&root; 
  int tok; 
  bool use_pr = (pr != "?");

  while (1) {
    if (use_pr) {
      if (!get_input(w, pr, 1, line)) return 0; 
      use_pr = false; 
    } else {
      pr2 = sofar + "(" + list_options(menu_ptr) + ") ? ";
      if (!get_input(w, pr2, 1, line)) return 0; 
    }
    tok = get_menu_item(first_word(w),menu_ptr); 
    if (tok) return tok; 
    sofar += menu_ptr->word + ' '; 
  }
}

// Return values for tok as in read_input. 

string text_menu::expand_input(const string& pr, int& tok) const
{
  string w, sofar, pr2; 
  const text_menu_item* menu_ptr = &root; 
  bool use_pr = (pr.length() != 0);
  tok = 0; 

  while (1) {
    if (use_pr) {
      if (!get_input(w, pr)) return sofar;
      use_pr = false; 
    } else {
      pr2 = sofar + "(" + list_options(menu_ptr) + ") ? ";
      if (!get_input(w, pr2)) return sofar; 
    }
    tok = get_menu_item(first_word(w), menu_ptr); 
    if (tok >= 0) sofar += menu_ptr->word; 
    else sofar += w; 
    if (tok) return sofar; 
    sofar += ' '; 
  }
}

string text_menu::list_options(const text_menu_item* menu_ptr) const
{
  string buf;

  menu_ptr = menu_ptr->submenu; 
  if (!menu_ptr) return buf; 

  while (1) {
    buf += menu_ptr->word;
    if (!(menu_ptr = menu_ptr->next)) break;
    buf += ' '; 
  }
  return buf;
}

static int pop_word(string& in, string& word)
{
  return split(in, &word, 1, in); 
}

void text_menu::add_item(string in, int token)
{
  text_menu_item *menu_ptr=&root; 

  string word; 
  while (in.length()) {
    if (!menu_ptr->submenu) break; 
    menu_ptr = menu_ptr->submenu; 

    /* look for the word on this level */ 
    pop_word(in, word);
    while (word != menu_ptr->word && menu_ptr->next) 
      menu_ptr=menu_ptr->next;

    if (word != menu_ptr->word) {
      menu_ptr = menu_ptr->next = new text_menu_item(word); 
      break; 
    }
  }

  while (pop_word(in,word)) {
    menu_ptr = menu_ptr->submenu = new text_menu_item(word); 
  }
  
  menu_ptr->token = token; 
}

text_menu_item::text_menu_item(const string& w) 
: word(w), next(0), submenu(0), token(0)
{}

text_menu_item::~text_menu_item() {
  delete next; 
  delete submenu; 
}

void get_name_path(string const& file, string& name, string& path)
{
  int i = file.rfind('/');
  name = file.substr(i+1);
  path = (i>0) ? file.substr(0,i) : string();
}

void get_base_ext(string const& file, string& base, string& ext)
{
  int i = file.rfind('.');
  ext = (i>0) ? file.substr(i+1) : string();
  base = file.substr(0,i);
}



/* $Id: menu.cc,v 1.2 2009/12/04 21:58:54 matthiasgoerner Exp $ */ 
