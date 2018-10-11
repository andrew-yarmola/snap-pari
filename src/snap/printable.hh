#ifndef _printable_hh_
#define _printable_hh_

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

#include <iostream>

class printable {
public:
  virtual void print(std::ostream& out) const =0;
};

inline std::ostream& operator << (std::ostream& out, printable const& p)
{
  p.print(out); return out;
}

template <class POb>
class printable_object : public printable {
  POb const& x;
public:
  printable_object(POb const& _x) : x(_x) {}
  virtual void print(std::ostream& out) const
  { x.print(out); }
};

template <class POb> 
inline printable_object<POb> PF(POb const& x) 
{
  return printable_object<POb>(x);
}

class seq_delims {
public:
  char open, *sep, close;
  seq_delims(const char* delims);
};

template <class Seq>
class printable_sequence : public printable, private seq_delims {
  Seq const& S;
public:
  printable_sequence(Seq const& S0, const char* delims)
    : seq_delims(delims), S(S0) {}

  virtual void print(std::ostream& out) const; 
};

template <class Seq>
void printable_sequence<Seq>::print(std::ostream& out) const
{
  typename Seq::const_iterator it;
  if (open) out << open;
  for (it=S.begin(); it!=S.end(); ++it) {
    if (it!=S.begin()) out << sep;
    out << (*it);
  }
  if (close) out << close;
}

template <class Seq> 
inline printable_sequence<Seq> PSeq(Seq const& S, const char* delim="[, ]")
{
  return printable_sequence<Seq>(S, delim);
}

template <class Seq>
class p_printable_sequence : public printable, private seq_delims {
  Seq const& S;
public:
  p_printable_sequence(Seq const& S0, const char* delims, bool esep=false)
    : seq_delims(delims, esep), S(S0) {}

  virtual void print(std::ostream& out) const; 
};

template <class Seq>
void p_printable_sequence<Seq>::print(std::ostream& out) const
{
  typename Seq::const_iterator it;
  if (open) out << open;
  for (it=S.begin(); it!=S.end(); ++it) {
    if (it!=S.begin()) out << sep;
    out << PF(*it);
  }
  if (close) out << close;
}

template <class Seq> 
inline p_printable_sequence<Seq> PPSeq(Seq const& S, const char* delim="[, ]")
{
  return p_printable_sequence<Seq>(S, delim);
}

template <class Seq>
class line_printable_seq : public printable {
  Seq const& S;
public:
  line_printable_seq(Seq const& S0)
    : S(S0) {}

  virtual void print(std::ostream& out) const 
  {
    typename Seq::const_iterator it;
    for (it=S.begin(); it!=S.end(); ++it)
      out << *it << std::endl; 
  }
};

template <class Seq>
inline line_printable_seq<Seq> LPSeq(Seq const& S)
{
  return line_printable_seq<Seq>(S);
}

#endif
