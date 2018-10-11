#ifndef _pfunc_
#define _pfunc_
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


#include <functional>

template <class Arg, class Result>
class unary_vfunction : public std::unary_function<Arg, Result> {
public:
  virtual Result operator() (const Arg& a) const = 0; 
};

template <class Arg1, class Arg2, class Result>
class binary_vfunction : public std::binary_function<Arg1, Arg2, Result> {
public:
  virtual Result operator() (const Arg1& x, const Arg2& y) const = 0; 
};

template <class Arg1, class Arg2, class Arg3, class Result>
class ternary_vfunction {
public:
  virtual Result operator() (const Arg1&, const Arg2&, const Arg3&) const = 0; 
};

#endif
