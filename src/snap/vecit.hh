#ifndef _vecit_hh_
#define _vecit_hh_
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

#include <iterator>

template <class VEC>
class const_vecit : 
  public std::iterator<std::random_access_iterator_tag, 
		       typename VEC::value_type, int> {
  const VEC* V;
  int pos; 
public:
  const_vecit(VEC const& v, int p) : V(&v), pos(p) {}

  // Forward iterator. 
  const_vecit() : V(0), pos(0) {}
  
  typedef typename VEC::value_type value_type; 

  value_type  operator * () { return  (*V)[pos]; }
  // value_type* operator ->() { return &(*V)[pos]; }

  const_vecit& operator ++ () 
  { ++pos; return *this; }
  const_vecit  operator ++ (int) 
  { const_vecit tmp; ++pos; return tmp; }

  friend bool operator == (const_vecit const& a, const_vecit const& b)
  { return a.V==b.V && a.pos == b.pos; }
  friend bool operator != (const_vecit const& a, const_vecit const& b)
  { return !(a==b); }

  // Bidirectional iterator.
  const_vecit& operator -- () 
  { --pos; return *this; }
  const_vecit  operator -- (int) 
  { const_vecit tmp(*this); --pos; return tmp; }

  // Random access iterator.
  const_vecit& operator += (int n) 
  { pos += n; return *this; }
  const_vecit& operator -= (int n) 
  { pos -= n; return *this; }
  
  friend const_vecit operator + (const_vecit const& a, int n)
  { const_vecit tmp(a); return tmp += n; }
  friend const_vecit operator + (int n, const_vecit const& a)
  { const_vecit tmp(a); return tmp += n; }
  friend const_vecit operator - (const_vecit const& a, int n)
  { const_vecit tmp(a); return tmp -= n; }
  friend int operator - (const_vecit const& a, const_vecit const& b)
  { return a.pos - b.pos; }

  value_type operator [] (int n) 
  { return (*V)[n+pos]; }

  friend bool operator < (const_vecit const& a, const_vecit const& b)
  { return a.pos < b.pos; }
  friend bool operator > (const_vecit const& a, const_vecit const& b)
  { return a.pos > b.pos; }
  friend bool operator <=(const_vecit const& a, const_vecit const& b)
  { return a.pos <= b.pos; }
  friend bool operator >=(const_vecit const& a, const_vecit const& b)
  { return a.pos >= b.pos; }
};

template <class VEC>
class vecit : 
  public std::iterator<std::random_access_iterator_tag, 
		       typename VEC::value_type, int> {
  VEC* V;
  int pos; 
public:
  vecit(VEC& v, int p) : V(&v), pos(p) {}

  operator const_vecit<VEC> () const 
  { return const_vecit<VEC>(*V,pos); }

  // Forward iterator. 
  vecit() : V(0), pos(0) {}
  
  typedef typename VEC::value_type value_type; 

  value_type& operator * () { return  (*V)[pos]; }
  value_type* operator ->() { return &(*V)[pos]; }

  vecit& operator ++ () 
  { ++pos; return *this; }
  vecit  operator ++ (int) 
  { vecit tmp; ++pos; return tmp; }

  friend bool operator == (vecit const& a, vecit const& b)
  { return a.V==b.V && a.pos == b.pos; }
  friend bool operator != (vecit const& a, vecit const& b)
  { return !(a==b); }

  // Bidirectional iterator.
  vecit& operator -- () 
  { --pos; return *this; }
  vecit  operator -- (int) 
  { vecit tmp(*this); --pos; return tmp; }

  // Random access iterator.
  vecit& operator += (int n) 
  { pos += n; return *this; }
  vecit& operator -= (int n) 
  { pos -= n; return *this; }
  
  friend vecit operator + (vecit const& a, int n)
  { vecit tmp(a); return tmp += n; }
  friend vecit operator + (int n, vecit const& a)
  { vecit tmp(a); return tmp += n; }
  friend vecit operator - (vecit const& a, int n)
  { vecit tmp(a); return tmp -= n; }
  friend int operator - (vecit const& a, vecit const& b)
  { return a.pos - b.pos; }

  value_type& operator [] (int n) 
  { return (*V)[n+pos]; }

  friend bool operator < (vecit const& a, vecit const& b)
  { return a.pos < b.pos; }
  friend bool operator > (vecit const& a, vecit const& b)
  { return a.pos > b.pos; }
  friend bool operator <=(vecit const& a, vecit const& b)
  { return a.pos <= b.pos; }
  friend bool operator >=(vecit const& a, vecit const& b)
  { return a.pos >= b.pos; }
};


#endif
