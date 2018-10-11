#ifndef _file_array_
#define _file_array_
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

#include <cstdio>
#include <algorithm>
#include "vecit.hh"

template <class T>
class file_array {
  FILE* fp;
  int _size; 

public:
  file_array(const char* name);
  ~file_array() { if (fp) fclose(fp); }

  T operator [] (int n) const; 
  int size() const { return _size; }
  int find(T const& key, int& first) const;  

  typedef T value_type;
  typedef const_vecit<file_array> const_iterator;

  const_iterator begin() const 
  { return const_iterator(*this,0); }
  const_iterator end() const 
  { return const_iterator(*this,_size); }
};

template <class T>
file_array<T>::file_array (const char* name)
 : _size(0)
{
  fp = fopen(name, "rb");
  if (!fp) {
    return; 
  }

  if (fseek(fp, 0, SEEK_END)==0) {
    _size = ftell(fp)/sizeof(T);
  } else {
    printf("problem determining size of file_array\n");
  }
}

template <class T>
T file_array<T>::operator [] (int n) const
{
  T L; 
  if (n < 0 || n >= _size) return L; 

  if (fseek(fp, n*sizeof(T), SEEK_SET)) {
    printf("problem finding position %d in the file_array\n", n); 
    return L; 
  }
  if (!fread((void*)&L,sizeof(T),1,fp)) {
    printf("problem reading an entry from the file_array\n");
  }
  return L;
}

template <class T>
int file_array<T>::find(T const& key, int& first) const
{
  std::pair<const_iterator, const_iterator> R = 
    std::equal_range(begin(), end(), key);

  first = R.first - begin(); 
  return R.second - R.first; 
}

#endif
