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
#include <pari/pari.h>
#include "pari_ext.h"

GEN
lindeppart1(GEN x, long bit)
{
  long  tx=typ(x),lx=lg(x),ly,i,j,flag,av,tetpil,e;
  GEN   y,p1,p2,p3,p4,p5;
  
  if((tx<17)||(tx>18)) err(talker, "Wrong type for first arg in lindeppart1.");
  av=avma;p1=greal(x);p2=gimag(x);
  ly=(flag=!gcmp0(p2)) ? lx+3 : lx+2;
  p4=cgetg(lx,19);bit=(long)(bit/L2SL10);
  for(i=1;i<lx;i++)
  {
    p5=cgetg(ly,18);p4[i]=(long)p5;
    for(j=1;j<lx+1;j++) p5[j]=(i==j) ? un : zero;
    p5[lx+1]=lcvtoi(gshift((GEN)p1[i],bit),&e);
    if(flag) p5[lx+2]=lcvtoi(gshift((GEN)p2[i],bit),&e);
  }
  p5=(GEN)lllint(p4);
  tetpil=avma;
  y=gmul(p4,p5);
  y=gerepile(av,tetpil,y);
  return y;
}

GEN
lindeppart2(GEN part, GEN x, long bit)
{
  long  tx=typ(x),lx=lg(part)+1,ly,j,flag,av,tetpil,e;
  GEN   y,p2,p3,p4,p5;
  
  av=avma;p2=gimag(x);
  ly=lg((GEN)part[1]); 
  bit=(long)(bit/L2SL10);

  p5=cgetg(ly,18);
  for(j=1;j<lx;j++) p5[j]=(lx-1==j) ? un : zero;
  p5[lx]=lcvtoi(gshift(greal(x),bit),&e);
  if(ly==lx+2) p5[lx+1]=lcvtoi(gshift(p2,bit),&e);

  p4=gtomat(concat(part,p5));
  p5=gmul(p4,lllint(p4));p3=(GEN)p5[1];
  tetpil=avma;y=cgetg(lx,17);
  for(j=1;j<lx;j++) y[j]=lcopy((GEN)p3[j]);
  y=gerepile(av,tetpil,y);
  return y;
}

