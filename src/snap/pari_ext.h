#ifndef _pari_ext_
#define _pari_ext_
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


/* a version of polredabsall that also returns a generator for the field in terms 
   of each reduced polynomial. */ 

GEN integer_solve(GEN mx, GEN d);

GEN lindeppart1(GEN x, long bit);

GEN lindeppart2(GEN part, GEN x, long bit);


#endif

/* to add a function which gp can call you have to: 
   write source and, if it's part of snap and lives in snap dir, put a
   link to it in from the pari source dir. 
   add a declaration of the function you have defined in the above format
   to this file. 
   add a line to anal.c
   add a line to helpmessages.c
   update Makefile so it will compile and link the required object file. 
   type make; no arguments are required. 
*/ 
