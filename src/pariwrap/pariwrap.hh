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

#ifndef _pariwrap_
#define _pariwrap_

#include "pariclass.hh"

// Wrapper functions: for each Pari library function taking and/or
// returning the pointer type GEN we define a function taking and/or
// returning paris, as defined in pariclass.hh. These are converted by
// a perl script from paridecl.h.

#include "inlines.hh"


// We also want to be able to do anything with a pari that we can do
// with a C double. 

// Pari versions of the usual comparison and equality tests

inline int operator < (pari const& a, pari const& b)
{ return gcmp(pari::tmp(a).g, pari::tmp(b).g) < 0; } 
inline int operator > (pari const& a, pari const& b)
{ return gcmp(pari::tmp(a).g, pari::tmp(b).g) > 0; } 
inline int operator <= (pari const& a, pari const& b)
{ return gcmp(pari::tmp(a).g, pari::tmp(b).g) <= 0; } 
inline int operator >= (pari const& a, pari const& b)
{ return gcmp(pari::tmp(a).g, pari::tmp(b).g) >= 0; } 
inline int operator == (pari const& a, pari const& b)
{ return gegal(pari::tmp(a).g, pari::tmp(b).g); } 
inline int operator != (pari const& a, pari const& b)
{ return !gegal(pari::tmp(a).g, pari::tmp(b).g); } 

// Pari equivalents of everything in math.h

inline pari sin(pari const& x)
{ return gsin(pari::tmp(x).g,prec); }
inline pari cos(pari const& x)
{ return gcos(pari::tmp(x).g,prec); }
inline pari tan(pari const& x)
{ return gtan(pari::tmp(x).g,prec); }
inline pari asin(pari const& x)
{ return gasin(pari::tmp(x).g,prec); }
inline pari acos(pari const& x)
{ return gacos(pari::tmp(x).g,prec); }
inline pari atan(pari const& x)
{ return gatan(pari::tmp(x).g,prec); }
inline pari atan2(pari const& x, pari const& y)
{ return garg(complex(x,y), prec); }

inline pari sinh(pari const& x)
{ return gsh(pari::tmp(x).g,prec); }
inline pari cosh(pari const& x)
{ return gch(pari::tmp(x).g,prec); }
inline pari tanh(pari const& x)
{ return gth(pari::tmp(x).g,prec); }
inline pari asinh(pari const& x)
{ return gash(pari::tmp(x).g,prec); }
inline pari acosh(pari const& x)
{ return gach(pari::tmp(x).g,prec); }
inline pari atanh(pari const& x)
{ return gath(pari::tmp(x).g,prec); }

inline pari exp(pari const& x)
{ return gexp(pari::tmp(x).g,prec); }
inline pari log(pari const& x)
{ return glog(pari::tmp(x).g,prec); }
inline pari log10(pari const& x)
{ return log(x)/log(integer(10)); }
inline pari pow(pari const& x, pari const& y)
{ return gpui(pari::tmp(x).g,pari::tmp(y).g,prec); }

inline pari sqrt(pari const& x)
{ return gsqrt(pari::tmp(x).g, prec); }
inline pari ceil(pari const& x)
{ return gceil(pari::tmp(x).g); }
inline pari floor(pari const& x)
{ return gfloor(pari::tmp(x).g); }
inline pari fabs(pari const& x)
{ return gabs(pari::tmp(x).g, prec); }
inline pari abs(pari const& x)
{ return gabs(pari::tmp(x).g, prec); }

inline pari ldexp(pari const& x, long n)
{ return gshift(pari::tmp(x).g, n); }
inline pari frexp(pari const& x, long& n)
{ GEN h = pari::tmp(x).g; n = expo(h); return gshift(pari::tmp(x).g, -n); }
inline pari trunc(pari const& x)
{ return gtrunc(pari::tmp(x).g); }
inline pari modf(pari const& x, pari& n)
{ n = trunc(x); return (x-n); }
inline pari fmod(pari const& x, pari const& y)
{ return (x - (y * floor(x/y))); }

#undef max
#undef min

inline pari max(pari const& a, pari const& b)
{ return gmax(a, b); }
inline pari min(pari const& a, pari const& b)
{ return gmin(a, b); }

// Other useful functions 
inline pari pi()
{ return mppi(prec); }

inline unsigned long stack_size()
{ return top - bot; }
inline unsigned long stack_free()
{ return avma - bot; }

int startindex(GEN p, long& len);

inline long length(pari const& x)
{ GEN g = pari::tmp(x).g; long len; return (lg(g) - startindex(g,len)); }
inline long lgeff(pari const& x)
{ GEN g = pari::tmp(x).g; long len; return (lgef(g) - startindex(g,len)); }

#endif

/* $Id: pariwrap.hh,v 1.1.1.1 2006/06/24 09:07:50 oag Exp $ */ 

