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
/* C++ code for pariwrap. wrappers for pari to avoid working
   with the pari stack directly */ 

#include <stdio.h> 
#include "pariwrap.hh"

#ifdef PARI_2_2_OR_LATER
#include "pari_oldnames.hh"
#endif

/* global precision constant, also declared extern in pariclass.hh */

long prec;

/* gc_off is switched to 1 when we're in the middle of evaluating 
   some large expression to prevent garbage collecting at that time */
#ifdef WITH_GP
long gc_off = 1; 
#else
long gc_off = 0;
#endif

/* garbage_count is incremented each time a pari is discarded. when
   it is low garbage collection doesn't take place. once it exceeds 
   garbage_limit and the opportunity to collect garbage arises then 
   garbage is collected and garbage_count is set back to 0. */
long garbage_count = 0; 
int garbage_limit = 1000; 

/* for debugging only, set garbage_trap to a value at which bug first
   appears in garbage collection, then re-run program in debugger with
   a breakpoint on warn. garbage_ct is the number of times collect_garbage
   has been called since the start of the program. */
int garbage_trap = 0; 
int garbage_ct = 0; 
bool trace_garbage = false;

/* this will be set to a nonzero value if an interrupt is received.
   program may test this in the middle of a long computation and 
   abort it if it is nonzero. */ 
int pari_interrupt = 0; 

/* digits is the number of decimal digits of accuracy to which 
   calculations are performed */ 
long digits = 50; 

/* paritop is the top of the stack which pari calculations are done 
   on. the stack builds downwards from this so paritop should remain
   constant throughout the life of the program. it is set in initpari. */
static long *paritop = 0; 

/* pari::first and pari::last are static members of the pari class. 
   the pari class just serves to keep a global data structure with all
   pari objects on it. pari is derived from pari so every pari contains
   a pari. the main point of the pari class is so that whenever we want
   to iterate through all live paris.. eg for garbage collection, we
   can do this. */
pari *pari::first = 0;
pari *pari::last = 0;

#ifndef WITH_GP
static initpari ip; 
#endif

/* links_made keeps count of the number of links that have been made. 
   its only purpose seems to be for debugging */
static int links_made = 0; 

/* the definitions of some constants and almost constants that
   snappea uses. the constants that aren't constant aren't becuase they 
   have to be recomputed if the precision changes. 
*/ 

const pari pZERO           = gzero;
const pari pONE            = gun; 
const pari pI              = gi;
const pari pTWO            = gdeux; 
const pari pFOUR           = integer(4);
const pari pHALF           = pONE/pTWO;
const pari pQUARTER        = pONE/pFOUR; 
const pari pHUGE           = pow(real(10),integer(1000)); 

pari pPI          = pi();
pari pTWO_PI      = pTWO * pPI;
pari pPI_OVER_2   = pPI/pTWO; 
pari pEPSILON     = real(1)/pow(real(10),integer((3*digits)/4)); 

/* called when precision is changed */ 
void init_real_globals()
{
  pPI = pi();
  pTWO_PI = pTWO * pPI;
  pPI_OVER_2 = pPI/pTWO; 
  pEPSILON = real(1)/pow(real(10),integer((3*digits)/4)); 
}

/* in_stack tells us whether the pointer p points to something within
   the limits of the pari stack. it seems that several "constant" GEN's
   are created above paritop which shouldn't be touched so we will often
   have pointers to these but we shouldn't follow them and modify anything. */
inline int in_stack(GEN p)
{ 
  return (p - paritop < 0 && p - (GEN)avma >= 0); 
}

/* startindex tells us whether a GEN contains pointers to other GEN's
   which are part of its value: eg. a vector contains pointers to several
   other GEN's which are its components. when we're collecting garbage
   we need to be able to track down all the GEN's a given GEN depends on. */
int startindex(GEN p, long& len)
{
  switch (typ(p)) {
  case 3: case 4: case 5: case 6:
  case 8: case 9: case 13: case 14:
  case 16: case 17: case 18: case 19:
    len = lg(p);
    return 1; 
    break;
  case 7: case 11:
    len = lg(p); 
    return 2; 
    break;
  case 10: 
	len = lgef(p); 
	return 2;
	break;
  default:
    break;
  }
  return 0; 
}

static bool startindex(GEN p, long& n, int& s, int i)
{
  s = startindex(p, n); 
  if (!s) {
    warn("attempt to access member of a non-composite pari type\n");
    return false; 
  }
  if (s + i >= lg(p) || i < 0) {
    char message[50]; 
    sprintf(message, "attempt to access slot: %d of pari with %ld slot(s)\n", i, lg(p) - s);
    warn(message); 
    return false;
  }
  return true; 
}

/* for debugging */ 
static int shallow_copy_count = 0; 
static int show_shallow_copies = 0; 

/* interrupt handling */

static void (*old_handler)(int) = 0; 

void my_handler(int val)
{
  val = val; 
  pari_interrupt = 1; 
  signal(SIGINT,old_handler);
}

void catch_interrupts()
{
  void (*handler)(int) = signal(SIGINT,my_handler);  

  /* avoid clobbering old_handler if this is called twice */ 
  if (handler != my_handler) old_handler = handler; 
  pari_interrupt = 0; 
}


void pari::link() // assume unlinked so that prev=0; 
{
#ifdef DEBUG
  assert(!linked());
#endif
  links_made++; 
  next = 0; 
  if (!last) { /* the first pari */ 
    first = last = this; 
    return; 
  }
  last->next = this;
  prev = last; 
  last = this; 
}

/* remove from the linked list and adjust all pointers accordingly */ 
void pari::unlink() // assume linked. 
{
#ifdef DEBUG 
  assert(linked()); 
#endif
  if (prev) {

    if (prev->next != this) {
      shallow_copy_count++; 
      if (show_shallow_copies) {
	warn("discarding shallow copy of a pari: \n"); 
	(this)->print(); 
      }
      return; /* get out without corrupting the linked list */ 
    }

    prev->next = next; 
  } else { 
    first = next; 
  }

  if (next) next->prev = prev; 
  else last = prev; 
  
  prev = next = 0; 
  links_made--; 
}

void pari::set_owner() 
{
  if (linked()) {
    warn("illegal attempt to set_owner for pari on stack\n"); 
    return; 
  }
  next = (pari*)true; 
  CH_S;
}

/* pari defs */ 
/* set_precision sets a new precision (in 32-bit words) for pari/pariwrap. 
   one side effect of this is that all global variables that a program
   uses may have to be reinitialized. the user has to provide
   a function to do this and call it at the appropriate time. */ 

void set_precision(long pr)
{
  if (pr < 1) pr = 1; 

  if (sizeof(long) == 8)
	  prec = (pr+1)/2 + 2; 
  else
	  prec = pr + 2; 

  digits = (long)(pariK * (prec-2)); 
}

/* initpari gets pari ready for work by calling init which sets the 
   stack size and everything. it also sets the precision to which pari
   computations will be done */ 

initpari::initpari()
{ 
  static bool done = false; 
  if (done) return; 
  long size = 20; 
  pari_init(1000000 * size, 500000); 
  lisseq("x;y;z"); 
  paritop = (long*)avma; 
  printf("pari initialized\n");
  done = true; 
  set_precision(6);
}

initpari::~initpari()
{
  if (!shallow_copy_count) return; 
  char message[50]; 
  sprintf(message, "shallow copies found: %d\n",shallow_copy_count);
  warn(message); 
}

pari::pari(GEN g) : p(g), prev(0), next(0)
{ 
  if (!p) { 
    warn("pari::p set NULL by pari(GEN) constructor\n"); 
    p = gzero; 
  }
  if (in_stack(p)) link(); 
  CH_S;
}

pari::pari(const pari& i) : p(i.p), prev(0), next(0)
{ 
  // printf("pari copy\n");
  if (i.shallow()) warn("pari copying a shallow copy\n"); 
  if (i.owner()) p = forcecopy(i.p); 
  if (in_stack(p)) link(); 
}

pari::~pari() 
{ 
  if (linked()) {
    unlink(); 
    garbage_count++; 
    // printf("pari destructor unlinking\n");

    if (garbage_count > garbage_limit && !gc_off) collect_garbage();

  } else {
    if (owner()) {
      // printf("pari clone destructor\n");
      gunclone(p); 
    } 
  }
}

void pari::clone() 
{
  // printf("pari clone\n");  
  if (!linked()) return; 

  unlink(); 
  p = gclone(p); 
  set_owner(); 
  CH_S; 
}

void pari::unclone() 
{
  // printf("pari unclone\n");
  if (!owner()) return; 

  GEN t = p; 
  p = forcecopy(t); 
  gunclone(t); 
  link(); 
  CH_S; 
}

void pari::clone_all()
{
  pari *k = first;
  pari *t; 
  while (k != 0) {
    t = k->next;
    k->clone(); 
    k = t;
  }
  avma = (long)paritop; // the stack should be empty now! 
}

/* components of a composite pari object can be accessed using
   the syntax vector[n]. vectors are numbered starting from 0 like
   standard C arrays. */

pari pari::operator [] (int i) const
{
  long n; 
  int s; 
  if (!startindex(p, n, s, i)) return gzero; 

  if (typ(p) == 10 && s + i >= n) {
    // polynomials are a special case:
    return gzero; 
  }

  GEN g = (GEN)p[i+s];

  // We can't return a pari that points to a bit of something cloned,
  // and since we can't distinguish between bits of cloned ones and universals
  // we have to take the safe option and always copy onto the stack. 
  if (!in_stack(g)) return forcecopy(g); 

  return g; 
}

pari_slot pari::operator [] (int i)
{
  long n; 
  int s; 
  if (!startindex(p, n, s, i)) return gzero; 

  if (owner()) unclone(); 

  if (typ(p) == 10 && s + i >= n) { 
    // polynomials are a special case:
    // if we're increasing the effective length, fill with zeros. 
    while (n <= s + i) p[n++] = (long)gzero; 
    setlgef(p, n);
  }

  return pari_slot(p, i+s);
}

// g must be universal (shared) or on the stack. 

void pari::set_to(GEN g)
{
  if (owner()) gunclone(p); 
 
  p = g; 
  if (in_stack(g)) {
    if (!linked()) link();
  } else {
    if (linked()) unlink(); 
  }    

  CH_S;
}

void pari::operator = (const pari& x)
{
  if (x.owner()) 
    set_to(forcecopy(x.p)); 
  else 
    set_to(x.p); 
}
  
void pari_slot::set_to(GEN g) // g must be universal or on stack. 
{ 
  if (!pp) return; // Null slot, error message already given. 
  pp[k] = (long)g; 
} 

/* diagnostics */ 

void set_shallow_copy_handling(int value)
{
  show_shallow_copies = value; 
}

/* compare the number of pari actually found with the number supposed
   to have been created according to the counter links_made. */ 

void pari::count_links()
{
  printf("links_made = %d", links_made); 

  pari* l = first;
  int ct = 0; 
  int loop = 2000; 
  while(l) {
    ct++; 
    l = l->next;
    if (loop-- < 0) break; 
  }
  
  if (loop < 0) 
    printf("there is a loop in the paris\n");
  else { 
    if (links_made!=ct) 
      printf("; actually have %d paris\n", ct);
    else 
      printf("\n"); 
  }
}

// We use the clone bit to mark which objects on the stack
// are in use during garbage collection. 

inline bool in_use(GEN g)
{
  return isclone(g)!=0; 
}

inline void set_in_use(GEN g)
{
  setisclone(g);
}

inline void set_not_in_use(GEN g)
{
  unsetisclone(g); 
}


// Without a field in which to do reference counting we 
// can't do this. 

#if 0
void count_refs(GEN g)
{
  pari *l;
  int irefs, erefs, i; 
  long n; 
  GEN sp; 

  erefs = irefs = 0; 
  l = first;
  while (l)
    {
      if (l->getp() == g) erefs++; 
      l = l->next; 
    }

  sp = (GEN)avma; 
  while (sp < paritop) {
    if (pere(sp) && (i = startindex(sp, n))) {
      while (i < n) {
	if (g == (GEN)sp[i]) irefs++; 
	i++; 
      }
    }
    sp += lg(sp); 
  }

  if (pere(g) != (irefs+erefs))
    printf("pari:%p pere:%d actual internal:%d actual external:%d\n",
	   g, pere(g), irefs, erefs);
}
#endif 

    
void pari::print_all()
{
  pari *k; 

  k = first;
  while (k != 0)
    {
      k->print();
      k = k->next;
    }
}

void pari::print_forward(pari *k)
{
  while (k)
    {
      k->print();
      k = k->next;
    }
}

void pari::print_back(pari *k)
{
  while (k)
    {
      k->print();
      k = k->prev;
    }
}
  

void pari::print_stack()
{
  long l;
  GEN g = (GEN)avma;
  int n = 20; 

  while (paritop - g > 0 && --n > 0) {
    printf("in_use flag:%d length:%ld type:%ld\n", 
	   in_use(g), l = lg(g), typ(g));
    outbeaut(g);
    g += l;
  }
}

void safe_info(GEN g)
{
  printf("address: %p in_use flag:%d length:%ld type:%ld\n", 
	 g, in_use(g), lg(g), typ(g));
}

void info(GEN g)
{
  safe_info(g); 
  outbeaut(g); 
}

void stack_info(GEN g)
{
  long l;
  int n = 20; 

  while (paritop - g > 0 && --n > 0) {
	safe_info(g); 
    g += l;
  }
}




/* useful constructors */


pari integer(long n)
{
  GEN g; 

  gaffsg(n, g = cgeti(3)); 

  return g;
}

pari polymod(const pari& a, const pari& b)
{
  GEN g = cgetg(3,9); 
  g[1] = (long)b.P();
  g[2] = (long)a.P(); 
  return g;
}

pari rvector(int n)
{
  GEN g = cgetg(n + 1, 17); 
  for (int i = 1; i <= n; i++) g[i] = (long)gzero; 
  return g; 
}

pari cvector(int n)
{
  GEN g = cgetg(n + 1, 18); 
  for (int i = 1; i <= n; i++) g[i] = (long)gzero; 
  return g; 
}

pari matrix(int c, int r)
{
  GEN g = cgetg(c + 1, 19), h; 
  int i,j;

  for (i = 1; i <= c; i++) {
    h = cgetg(r + 1, 18);
    for (j = 1; j <= r; j++) h[j] = (long)gzero; 
    g[i] = (long)h;
  }

  return g; 
}

pari complex(const pari& r, const pari& i)
{
  GEN cx = cgetg(3,6);
  cx[1] = (long)r.P();
  cx[2] = (long)i.P();
  return cx;
}

pari polynomial(long degree)
{
  GEN g = cgetg(degree+3,10);
  setsigne(g, 0);
  setlgef(g, 2); 
  setvarn(g, 0);
  int i;
  for (i = 2; i < degree+3; i++)
    g[i] = (long)gzero; 

  return g; 
}

void set_poly(pari& p)
{
  GEN g = p.P(); 
  int i = lg(g) - 1; 

  while (i > 1 && isexactzero((GEN)g[i])) i--; 
  int len = i+1;
  while (i > 1 && gcmp0((GEN)g[i])) i--; 
  setsigne(g, i!=1); 
  setlgef(g, len); 
}

pari real(int n)
{
  GEN g; 
  gaffsg(long(n), g = cgetr(prec)); 
  return g;
}

pari real(long n)
{
  GEN g; 
  gaffsg(n, g = cgetr(prec)); 
  return g;
}

pari real(double x)
{
  return dbltor(x);
}

int numerical_zero(const pari& r)
{ 
  if (r.type()==t_REAL || r.type()== t_COMPLEX) 
    return abs(r) < pEPSILON; 
  return isexactzero(r); 
}

pari psign(const pari& r)
{
  int rt = r.type(); 
  if (rt != t_REAL && rt != t_INT) {
    printf("invalid type in psign\n");
    return pZERO; 
  }
  if (numerical_zero(r)) return pZERO; 
  return r > pZERO ? pONE : -pONE; 
}

/* garbage collection */ 

class to_keep {
  to_keep *prev;
public:
  long *from, *to, offset; 
  to_keep(to_keep* pr) { prev = pr; }
  ~to_keep() { if (prev) delete prev; }
  to_keep* next() { return prev; }
};

void print_kbs(to_keep *kg)
{
  while (kg) {
    printf("top: %p bottom: %p -> top: %p bottom %p\n", 
		   kg->to, kg->from, kg->to + kg->offset, kg->from + kg->offset); 
    kg = kg->next(); 
  }
}

static void adjust_pointer(GEN& p, to_keep *kglast)
{
  to_keep *kg; 
  
  if (!in_stack(p)) return; 
  kg = kglast; 
  while (kg) {
    if (p - kg->from >= 0) {
      if (p - kg->to >= 0) 
		warn("found pointer to garbage\n");
      p += kg->offset; 
      break;
    }
    kg = kg->next(); 
  }
}

/* mark a GEN to keep */ 
void mark_to_keep(GEN p)
{
  if (!in_stack(p) || in_use(p)) return; 

  set_in_use(p);
  long n; 
  int i = startindex(p, n); 
  if (i) {
    while (i < n) 
      {
		mark_to_keep((GEN)(p[i]));
		i++; 
      }
  }
}  

void pari::mark_to_keep()
{ 
  ::mark_to_keep(p); 
}

static long 
fix_size(long a)
{
  /* BYTES_IN_LONG*ceil(a/BYTES_IN_LONG) */
  ulong b = a+BYTES_IN_LONG - (((a-1) & (BYTES_IN_LONG-1)) + 1);
  if (b > VERYBIGINT) err(talker,"stack too large");
  return b;
}

void pari::resize_stack(long newsize)
{
  if (paritop != (long*)top) {
    warn("paritop != top, don't know how to resize the stack.\n");
    return; 
  }

  long sizeold = top - bot;

  if (!newsize)
  {
    newsize = sizeold << 1;
    err(warner,"doubling stack size; new stack = %.1f MBytes",newsize/1E6);
  }
  else if ((long)newsize < 1000L)
    err(talker,"required stack memory too small");

  long size = fix_size(newsize);

  ulong newbot = (ulong) gpmalloc(size);
  ulong newtop = memused = newbot+size;
  
  collect_garbage((long*)newtop);
  
  free((void*)bot); 
  bot = newbot;
  top = newtop;
  paritop = (long*)top; 
}

/* If newtop is nonzero, the whole stack is copied to a new location 
   during garbage collection. Sufficient memory should already have 
   been allocated below newtop. It is up to the calling function to 
   reset the variables bot, top, and memused and to free any 
   memory belonging to the old stack. This function modifies avma. 
*/ 

void pari::collect_garbage(long* newtop)
{
  long *gtop, current_offset, *rp, n, i; 
  GEN p, wp; 
  to_keep *kg, *kglast = 0; 
  int stuff_moved = 0; 
  pari *k;

  // int print_stack=0; // can set in debugger.

  if (!newtop) newtop = paritop; 

  if (++garbage_ct == garbage_trap) warn("garbage trap triggered\n"); 

  if (trace_garbage) printf("Collecting garbage.\n"); 

  // printf("."); fflush(stdout); 

  // if (print_stack) print_all(); 

#if 0
  /* look for bad internal refs */
  p = (GEN)avma;
  while (paritop - p > 0) {
    if (in_use(p) && (i = startindex(p, n))) {
      while (i < n) {
	if (!in_use((GEN)p[i])) {
	  char message[50];
	  sprintf(message, "found bad internal reference at:%x\n", p+i);
	  warn(message); 
	}
	i++; 
      }
    }
    p += lg(p);
  }

  /* look for bad external refs */ 
  k = first;
  while (k)
    {
      if (!in_use(k->p)) {
	char message[50]; 
	sprintf(message, "found pointer to garbage at:%x\n", k->p);
	warn(message); 
	exit(3);
      }
      k = k->next;
    }

  if (gc_off) {
    warn("collect_garbage called while garbage collection was disabled\n"); 
    return;
  }

  if (pari::last->next) {
    warn("linked list corrupt, last pari has a successor\n"); 
    return; 
  }
#endif

  /* go through stack, set in_use flag to false */
  p = (GEN)avma; 
  while (p < paritop) {
    set_not_in_use(p); 
    p += lg(p);
  }

  /* go through paris incrementing refcount for each pari */
  k = first;
  while (k != 0)
    {
      k->mark_to_keep();
      k = k->next;
    }

  /* make a record of what we want to keep */
  p = (GEN)avma; 
  while (p < paritop) {
    /* go to the end of the garbage */ 
    while (p < paritop && !in_use(p)) p += lg(p);
    if (p >= paritop) break; 
    kglast = new to_keep(kglast); 
    kglast->from = p;
    while (p < paritop && in_use(p)) p += lg(p);
    kglast->to = p; 
  }

  /* compute the offset for each block of stuff to keep */ 
  kg = kglast; 
  gtop = paritop; 
  current_offset = newtop - paritop; 
  while (kg) {
    current_offset += (gtop - kg->to);
    kg->offset = current_offset; 
    gtop = kg->from; 
    kg = kg->next(); 
  }

  /* copy everything we want to keep */
  kg = kglast; 
  wp = newtop; 
  while (kg) {
    if (kg->offset == 0) 
      { wp = kg->from; kg = kg->next(); continue; }
    rp = kg->to;
    while (rp - kg->from > 0) { *(--wp) = *(--rp); }
    stuff_moved = 1; 
    kg = kg->next(); 
  }

  /* do we need to check pointers? */ 
  if (!stuff_moved) { 
    delete kglast;   
    avma = (long)wp; 
    return; 
  }

  /* adjust internal pointers */
  p = wp;
  while (newtop - p > 0) {
    i = startindex(p, n);
    if (i) {
      while (i < n) {
	adjust_pointer(((GEN*)p)[i], kglast);
	i++; 
      }
    }
    p += lg(p);
  }

  /* adjust external pointers */
  k = first;
  while (k != 0)
    {
      adjust_pointer(k->p, kglast);
      k = k->next;
    }

  avma = (long)wp; 

  /* give back memory we used */ 
  delete kglast;

  // if (print_stack) print_all(); 

  garbage_count = 0; 
}

#if 0

struct stack_state {
  long *_paritop;
  long _avma;
  stack_state *prev; 
};

static stack_state *saved_state = 0; 

void push_stack_state()
{
  collect_garbage();
  stack_state *state_now = new stack_state; 
  state_now->prev = saved_state; 
  state_now->_paritop = paritop; 
  state_now->_avma = avma; 
  paritop = (long*)avma; 
  saved_state = state_now; 
}

void pop_stack_state()
{
  if (!saved_state) return; 
  paritop = saved_state->_paritop; 
  stack_state *t = saved_state; 
  saved_state = saved_state->prev; 
  delete t; 
}

void restore_stack_state()
{
  if (!saved_state) return; 
  paritop = saved_state->_paritop; 
  avma = saved_state->_avma; 
}
#endif

#if 0 
{
  /* set any pari's left dangling to zero */ 
  int dangling = 0; 
  pari *t, *k = first; 
  while (k != 0) {
	t = k; 
	k = k->next; 
	if (!in_stack(t->p)) {
	  t->p = gzero; 
	  dangling++; 
	}
  }
  /* give a warning if any dangling paris were found */ 
  if (dangling) {
	char s[100];
	sprintf(s,"%d dangling paris set to 0\n", dangling); 
	warn(s); 
  }
}
#endif

void pari::zero_if_garbage()
{ 
  if (!in_use(p)) p = gzero; 
}

/* the following is specifically for turning pariwrap code into
   code that can be installed into gp. so given pari func(pari a,...)
   we do: 

   GEN gpfunc(GEN a,...,long pr) { GEN av = (GEN)avma; setprecr(pr); 
     return collect_garbage(av,func(a,...)); }

   it uses the following global variables: 
     avma, pari::first, pari::last
   */

GEN pari::collect_garbage(long *to, pari keep) 
{
  long *gtop, current_offset, *rp, n, i, *old_ptop; 
  GEN p, wp; 
  to_keep *kg, *kglast = 0; 
  int stuff_moved = 0; 
  pari *k;

/*
  static int ct = 0; 

  if (++ct == garbage_trap) warn("garbage trap triggered\n"); 
*/
  /* 
  if (pari::last && pari::last->next) {
    warn("linked list corrupt, last pari has a successor\n"); 
    return 0; 
  }
  */

  old_ptop = paritop; 
  paritop = to; 

  /* go through stack, set refcounts to 0 */
  p = (GEN)avma; 
  while (p < to) {
    set_not_in_use(p);
    p += lg(p);
  }

  /* set refcount nonzero for keep and anything depending on it */ 
  keep.mark_to_keep(); 

  /* point all variables we're not going to keep at 0 */ 
  k = first;
  while (k != 0)
    {
      k->zero_if_garbage();
      k = k->next;
    }

  /* make a record of what we want to keep */
  p = (GEN)avma; 
  while (p < to) {
    /* go to the end of the garbage */ 
    while (p < to && !in_use(p)) p += lg(p);
    if (p >= to) break; 
    kglast = new to_keep(kglast); 
    kglast->from = p;
    while (p < to && in_use(p)) p += lg(p);
    kglast->to = p; 
  }

  /* compute the offset for each block of stuff to keep */ 
  kg = kglast; 
  gtop = to; 
  current_offset = 0; 
  while (kg) {
    current_offset += (gtop - kg->to);
    kg->offset = current_offset; 
    gtop = kg->from; 
    kg = kg->next(); 
  }

  /* copy everything we want to keep */
  kg = kglast; 
  wp = to; 
  while (kg) {
    if (kg->offset == 0) 
      { wp = kg->from; kg = kg->next(); continue; }
    rp = kg->to;
    while (rp - kg->from > 0) { *(--wp) = *(--rp); }
    stuff_moved = 1; 
    kg = kg->next(); 
  }

  /* do we need to check pointers? */ 
  if (!stuff_moved) { 
    delete kglast;   
    avma = (long)wp; 
	paritop = old_ptop; 
    return keep.P(); 
  }

  /* adjust internal pointers */
  p = wp;
  while (to - p > 0) {
    i = startindex(p, n);
    if (i) {
      while (i < n) {
	adjust_pointer(((GEN*)p)[i], kglast);
	i++; 
      }
    }
    p += lg(p);
  }

  /* adjust pointer in keep */
  adjust_pointer(keep.p, kglast); 

  avma = (long)wp; 

  /* give back memory we used */ 
  delete kglast;
  paritop = old_ptop; 
  return keep.P(); 
}

/* print styles are controlled by pari_printform objects */ 

pari_printform def_print = pari_printform(1,1,'g',-1,0); 

void pari::print() const 
{
  CH_S; 
  print(def_print); 
}

void pari::print(const pari_printform& pf) const
{
  CH_S; 
  if (pf.raw)
	brute(p,pf.format,pf.dec);
  else
	sor(p,pf.format,pf.dec,pf.field);

  if (pf.newline && (pf.raw || type() < 17 || type() > 19)) 
	fprintf(outfile, "\n"); 
}

#if 0 
pari lindep2(pari x, long bit)
{ 
  GEN g; 
  long av; 
  try { 
	av = avma; 
	g = lindep2(x.g(), bit); 
  } catch (my_int) {
	avma = av; 
	throw my_int(); 
  }
  return g; 
}

pari algdep2(pari x, long n, long bit)
{ 
  GEN g; 
  long av; 
  try { 
	av = avma; 
	g = algdep2(x.g(), n, bit); 
  } catch (my_int) {
	avma = av; 
	throw my_int(); 
  }
  return g; 
}
#endif 

/* $Id: pariwrap.cc,v 1.2 2009/12/03 21:52:43 matthiasgoerner Exp $ */ 
