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

#ifndef _pariclass_
#define _pariclass_

// C++ headers for pariwrap. Wrappers for pari to avoid working
// with the Pari stack directly.

#ifdef __DECCXX
#include <stdcomp>
#endif

#include <stdio.h> 
#include "warn.hh"

extern "C" {
#ifdef __DECCXX

#define exception math_exeption
#include <pari/pari.h>
#undef exception

#else 

#include <pari/pari.h>

#endif
}

#undef init 

#define DEBUG 1


// If compiling with an old version of Pari (e.g. 2.0.10.beta) 
// comment the following #define and uncomment out the one after. 

#define outfile pari_outfile
// #define sizedigit gsize

#undef min
#undef max

extern long prec; // Default precision for Pari to work at.

// While it is generally a bad idea to have global variables 
// we don't worry too much in this case: Pari already has plenty
// of global variables as well as a global stack for all its
// memory requirements. 

void set_precision(long); 

// Interrupt handling. 
extern int pari_interrupt;
void catch_interrupts();

// The following are of no concern to the user but are
// needed in the definition of the class pari. 
extern long digits; 
extern long gc_off; 
extern long garbage_count; 
extern bool trace_garbage;
extern int garbage_limit; 

#ifdef DEBUG
#include <assert.h>
#define CH_S if (shallow()) warn("shallow copy of pari_link encountered\n")
#define CH_G if (gc_off < 0) warn("gc_off is now negative\n")
#else
#define CH_S 
#define CH_G
#endif

// Pari has to be initialized before the program tries to create any
// variables of type pari. The only way to ensure that this occurs is
// to call Pari's initialization function during the global static
// initialization segment. This can be done by creating a global
// object of type initpari and putting the required call into its
// constructor. Ugly but invisible to the end user. 

class initpari {
public:
  initpari(); 
  ~initpari();
};

// We keep track of all the pari objects we create on a big linked
// list so that come garbage collection time we can review what we
// still want to keep.

// Encode print formatting information in a convenient structure. 

class pari_printform {
public:
  int raw; 
  int newline; 
  char format;
  long dec;
  long field; 

  pari_printform() 
	{ raw = 0; newline = 1; format = 'g'; dec = -1; field = 0; }
  pari_printform(int r, int n, char fo, long d, long fi)
	{ raw = r; newline = n; format = fo; dec = d; field = fi; }
};

extern pari_printform def_print; 

// Forward declaration of pari_slot. 

class pari_slot; 

// The pari class itself, encapsulating the pari pointer type GEN.
// Constructors and operator overloading enable it to be used pretty
// much the way built-in arithmetic types are used. 

class pari {
  static pari *first, *last; 

  GEN p; 
  pari *prev, *next;

  bool linked() const { return prev || this==first; }
  bool owner() const { return !linked() && next; } 
  void set_owner();

  void link();
  void unlink(); 
protected:
  virtual void set_to(GEN g);

  pari(GEN g, int) : p(g), prev(0), next(0) {}
public:
  pari() : p(gzero), prev(0), next(0) {}

  pari(GEN g);
  pari(const pari& i);

  pari(const int& n) : p(stoi(n)), prev(0) { link(); }
  pari(const long& n) : p(stoi(n)), prev(0) { link(); }
  pari(const unsigned int& n) : p(stoi(n)), prev(0) { link(); }
  pari(const unsigned char& n) : p(stoi(n)), prev(0) { link(); }
  pari(const double& x) : p(dbltor(x)), prev(0) { link(); }

  virtual ~pari();

  int int_value() const { return gtolong(p); }
  double double_value() const { return gtodouble(p); }
  char* string_value() const { return GENtostr(p); } // caller free()'s. 

  // The following is used to turn Pari functions of GEN 
  // into functions of pari. See inlines.hh. 

  class tmp {
  public:
    GEN g;
    
    tmp(GEN h) : g(h) { gc_off++; } 
    tmp(pari const& p) : g(p.p) { gc_off++; }
    tmp(tmp const& t) : g(t.g) { gc_off++; }
    ~tmp() { gc_off--; CH_G; }
  private:
    tmp(); 
  };
  
  GEN P() const { return p; } 
  
  void print() const; 
  void print(const pari_printform& pf) const;

  unsigned long type() const { return typ(p); }
  long length() const { return lg(p); }

  void operator = (const pari& x);

  pari operator [] (int n) const;
  pari_slot operator [] (int n);

  // Standard arithmetic operators 
  friend pari operator - (pari const& x)
	{ return gneg(tmp(x).g); }
  friend pari operator + (pari const& a, pari const& b)
	{ return gadd(tmp(a).g,tmp(b).g); }
  friend pari operator - (pari const& a, pari const& b)
	{ return gsub(tmp(a).g,tmp(b).g); }
  friend pari operator * (pari const& a, pari const& b)
	{ return gmul(tmp(a).g,tmp(b).g); }
  friend pari operator / (pari const& a, pari const& b)
	{ return gdiv(tmp(a).g,tmp(b).g); }
  friend pari operator % (pari const& a, pari const& b)
	{ return gmod(tmp(a).g,tmp(b).g); }

  void operator += (const pari& b)
    { set_to(gadd(P(),b.P())); }
  void operator -= (const pari& b)
    { set_to(gsub(P(),b.P())); }
  void operator *= (const pari& b)
    { set_to(gmul(P(),b.P())); }
  void operator /= (const pari& b)
    { set_to(gdiv(P(),b.P())); }
  void operator %= (const pari& b)
    { set_to(gmod(P(),b.P())); }

  static void resize_stack(long newsize);

  // Garbage collection, user may ignore. 
  static void collect_garbage(long *newtop=0);
  static GEN collect_garbage(long *to, pari keep);

  // static void restore_stack_state(); 

  void clone(); 
  void unclone(); 
  static void clone_all(); 

  // Debugging
  // GEN getp() { return p; }
  static void print_forward(pari *);
  static void print_back(pari *);
  static void count_links(); 
  static void print_stack();
  static void print_all(); 
  int shallow() const { return prev && (prev->next != this); } 

private:

  // More garbage collection. 
  void zero_if_garbage();
  void mark_to_keep();
};

int numerical_zero(const pari& r);
pari psign(const pari& r);

// Constructors: since Pari and hence pariwrap do not use separate C
// types for each mathematical type, not all required types can be
// constructed using bona-fide C++ constructors in the pari class. A
// range of convenient "external" constructors is therefore provided.

pari polymod(const pari& pol, const pari& mod);
pari cvector(int n);
pari rvector(int n);
pari matrix(int c, int r);
pari integer(long n);
pari real(int n);
pari real(long n);
pari real(double x);
pari complex(const pari& r, const pari& i);
pari polynomial(long degree); 
void set_poly(pari& p); 

// A pari_slot is a closure: 
// a function with its arguments ready to evaluate. If pari M
// is a vector then M[n] is a pari_slot which can be evaluated to
// return an entry, or used on the left hand side of an assignment, in
// which case the unevaluated information is required.

class pari_slot : public pari {
  GEN pp;
  long k;

  virtual void set_to(GEN g);

  pari_slot();
public:
  pari_slot(GEN g) : pari(g,0), pp(0), k(0) { gc_off++; } // null slot
  pari_slot(GEN g, long m) : pari((GEN)g[m],0), pp(g), k(m) { gc_off++; }
  pari_slot(pari_slot const& s) : pari(s.P(),0), pp(s.pp), k(s.k) { gc_off++; }

  virtual ~pari_slot() { gc_off--; CH_G; } 

  // default assignment not OK.
  void operator = (const pari_slot& x) { set_to(x.P()); }
  void operator = (const pari& x) { pari::operator = (x); }
};

// Constants 
extern const pari pZERO; 
extern const pari pONE; 
extern const pari pI;  
extern const pari pTWO;
extern const pari pFOUR;
extern const pari pHALF;
extern const pari pQUARTER;
extern const pari pHUGE; 

// May be recomputed when precision changes by calling init_real_globals().
extern pari pPI;
extern pari pTWO_PI;
extern pari pPI_OVER_2; 
extern pari pEPSILON; 

void init_real_globals();

class to_keep; 

// Diagnostics 
void info(GEN); 
void print_kbs(to_keep*);
void stack_info(GEN); 
void set_shallow_copy_handling(int value);

#endif

