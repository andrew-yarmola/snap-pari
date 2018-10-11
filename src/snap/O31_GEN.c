#ifndef NOPARI

#include "O31_GEN.h"

#ifndef stack_lim
#define stack_lim(av,n) (bot + (((av)-bot)>>(n)))
#endif

GEN
O1n_normsq(GEN x)
{
  GEN y;
  long i,tx=typ(x),lx,av,lim;

  if (! is_vec_t(tx)) err(talker, "not a vector in O1n_normsq");
  lx=lg(x); if (lx==1) return gzero;

  av=avma; lim = stack_lim(av,1); y = gneg(gsqr((GEN) x[1]));
  for (i=2; i<lx; i++)
  {
    y = gadd(y, gsqr((GEN) x[i]));
    if (low_stack(lim, stack_lim(av,1)))
    {
      if (DEBUGMEM>1) err(warnmem,"O1n_normsq");
      y = gerepileupto(av, y);
    }
  }
  return gerepileupto(av,y);
}

GEN
O1n_innerprod(GEN u, GEN v)
{
  GEN y;
  long i,lu,av,lim;

  if (!is_vec_t(typ(u)) || !is_vec_t(typ(v))) 
    err(talker, "not a vector in O1n_innerprod");
  lu=lg(u); 
  if (lu != lg(v)) 
    err(talker, "vectors not same length in O1n_innerprod"); 
  if (lu==1) return gzero;

  av=avma; lim = stack_lim(av,1); y = gneg(gmul((GEN) u[1], (GEN) v[1]));
  for (i=2; i<lu; i++)
  {
    y = gadd(y, gmul((GEN) u[i], (GEN) v[i]));
    if (low_stack(lim, stack_lim(av,1)))
    {
      if (DEBUGMEM>1) err(warnmem,"O1n_normsq");
      y = gerepileupto(av, y);
    }
  }
  return gerepileupto(av,y);
}

GEN 
O1n_inverse(GEN y)
{
  GEN m, col;
  long lm, nr, c, r; 

  /* Check the input */
  if (typ(y)!=t_MAT) err(talker, "not a matrix in O1n_inverse"); 

  lm = lg(y); 
  if (lm==1) return gcopy(y);
  nr = lg((GEN)y[1]); 

  /* Matrix is transposed and sign changed on first row and column */ 
  m = cgetg(nr, t_MAT);
  for (c=1; c<nr; c++) {
    col = cgetg(lm, t_COL);
    m[c] = (long)col; 
    for (r=1; r<lm; r++) {
      if ((c==1) == (r==1)) 
	col[r] = lcopy(gcoeff(y,c,r)); 
      else 
	col[r] = lneg(gcoeff(y,c,r));
    }
  }

  return m; 
}

GEN
O1n_GramSchmidt(GEN y, long prec)
{
  GEN m, svec, col, col0, mul, minusone;
  long lm, c, c0, r, nr, ltop; 

  /* Check the input */
  if (typ(y)!=t_MAT) err(talker, "not a matrix in O1n_GramSchmidt"); 

  lm = lg(y); 
  if (lm==1) return gcopy(y);
  nr = lg((GEN)y[1]); 

  /* First copy the matrix. Matrix must consist entirely of reals. */
  m = cgetg(lm, t_MAT);
  for (c=1; c<lm; c++) {
    col = cgetg(nr, t_COL);
    m[c] = (long)col; 
    for (r=1; r<nr; r++) {
      col[r] = (long)cgetg(prec+1, t_REAL); 
      gaffect(gcoeff(y,r,c), (GEN)col[r]); 
    }
  }

  ltop = avma; 

  /* Make vector to keep track of signs. */
  svec = cgetg(lm, t_VEC);
  minusone = gneg(gun); 

  for (c=1; c<lm; c++) {

    /* Orthogonalize a column wrt previous cols. */
    col = (GEN)m[c];
    for (c0=1; c0<c; c0++) {
      col0 = (GEN)m[c0];
      mul = gmul(O1n_innerprod(col0,col), (GEN)svec[c0]);

      /* col = col - mul*col0 */ 
      for (r=1; r<nr; r++)
	gaffect(gsub((GEN)col[r], gmul((GEN)col0[r],mul)), (GEN)col[r]); 
    }

    mul = O1n_normsq(col);

    /* Record the sign. */
    if (gsigne(mul) > 0) {
      svec[c] = un; 
      mul = gsqrt(mul, prec); 
    } else {
      svec[c] = (long)minusone; 
      mul = gsqrt(gneg(mul), prec); 
    }

    /* Normalize the column. */
    for (r=1; r<nr; r++)
      gaffect(gdiv((GEN)col[r], mul), (GEN)col[r]); 
  }

  avma = ltop; 
  return m; 
}

GEN
SL2C_to_O31(GEN m, long prec)
{
  static GEN M[4] = {0,0,0,0};
  GEN mct, A, img;
  long ltop, i, j, thedet;

  thedet = 1; 
  if (typ(m)==t_VEC) {
    thedet = gsigne((GEN)m[2]); 
    m = (GEN)m[1];
  }

  /* Make an uninitialized matrix for the result. */
  A = cgetg(5, t_MAT);
  for (i=1; i<5; i++) {
    A[i] = (long)cgetg(5, t_COL); 
    for (j=1; j<5; j++) coeff(A,j,i) = (long)cgetg(prec+1,t_REAL); 
  }

  ltop = avma; 
  if (!M[0]) {
    M[0] = gclone(lisexpr("[1,0;0,1]"));
    M[1] = gclone(lisexpr("[1,0;0,-1]"));
    M[2] = gclone(lisexpr("[0,1;1,0]"));
    M[3] = gclone(lisexpr("[0,I;-I,0]"));
    avma = ltop; 
  }

  mct = gtrans(gconj(m));

  for (i=1; i<5; i++) {
    img = gmul(mct, gmul(M[i-1], m));

    gaffect(gmul2n(gadd(greal(gcoeff(img,1,1)),greal(gcoeff(img,2,2))),-1), 
	    gcoeff(A,i,1));
    gaffect(gmul2n(gsub(greal(gcoeff(img,1,1)),greal(gcoeff(img,2,2))),-1), 
	    gcoeff(A,i,2));
    gaffect(gmul2n(gadd(greal(gcoeff(img,1,2)),greal(gcoeff(img,2,1))),-1), 
	    gcoeff(A,i,3));
    gaffect(gmul2n(gsub(gimag(gcoeff(img,1,2)),gimag(gcoeff(img,2,1))),-1), 
	    gcoeff(A,i,4));
  }
  if (thedet < 0) {
    gaffect(gneg((GEN)A[4]), (GEN)A[4]);
  }
  avma = ltop; 

  return A; 
}

GEN
O31_to_SL2C(GEN B, long prec)
{
  GEN A, AM0A_00, AM1A_00, aa, bb, R;
  long ltop, i, j, old_col4, thedet;

  R = cgetg(3, t_VEC);

  /* First make an uninitialized 2 by 2 matrix of complexes for the result. */ 
  A = cgetg(3, t_MAT); 
  for (i=1; i<3; i++) {
    A[i] = lgetg(3, t_COL);
    for (j=1; j<3; j++) {
      coeff(A,j,i) = lgetg(3, t_COMPLEX);
      mael3(A,i,j,1) = lgetg(prec+1, t_REAL);
      mael3(A,i,j,2) = lgetg(prec+1, t_REAL);
    }
  }

  R[1] = (long)A; 
  R[2] = (long)stoi(1); 
  /* Or pres until proven otherwise. This makes space to assign -1. */ 

  ltop = avma; 

  AM0A_00 = gadd(gcoeff(B,1,1),gcoeff(B,2,1));
  AM1A_00 = gadd(gcoeff(B,1,2),gcoeff(B,2,2));
  aa = gadd(AM0A_00, AM1A_00);
  bb = gsub(AM0A_00, AM1A_00);

  thedet = 1; 
  if (gsigne(det2(B))==-1) {
    old_col4 = B[4];
    B[4] = lneg((GEN)B[4]); 
    thedet = -1; 
  }

  if (gcmp(aa,bb) == 1) {

    // printf("aa bigger\n"); 
    gaffect(aa,    gmael3(A,1,1,1)); 
    gaffect(gzero, gmael3(A,1,1,2)); 

    gaffect(gadd(gcoeff(B,1,3),gcoeff(B,2,3)), gmael3(A,2,1,1)); 
    gaffect(gadd(gcoeff(B,1,4),gcoeff(B,2,4)), gmael3(A,2,1,2)); 

    gaffect(gadd(gcoeff(B,3,1),gcoeff(B,3,2)),       gmael3(A,1,2,1)); 
    gaffect(gneg(gadd(gcoeff(B,4,1),gcoeff(B,4,2))), gmael3(A,1,2,2)); 

    gaffect(gadd(gcoeff(B,3,3),gcoeff(B,4,4)), gmael3(A,2,2,1)); 
    gaffect(gsub(gcoeff(B,3,4),gcoeff(B,4,3)), gmael3(A,2,2,2)); 

  } else {

    // printf("bb bigger\n"); 
    gaffect(gadd(gcoeff(B,1,3),gcoeff(B,2,3)),       gmael3(A,1,1,1)); 
    gaffect(gneg(gadd(gcoeff(B,1,4),gcoeff(B,2,4))), gmael3(A,1,1,2)); 

    gaffect(bb,    gmael3(A,2,1,1)); 
    gaffect(gzero, gmael3(A,2,1,2)); 

    gaffect(gsub(gcoeff(B,3,3),gcoeff(B,4,4)),       gmael3(A,1,2,1)); 
    gaffect(gneg(gadd(gcoeff(B,3,4),gcoeff(B,4,3))), gmael3(A,1,2,2)); 

    gaffect(gsub(gcoeff(B,3,1),gcoeff(B,3,2)), gmael3(A,2,2,1)); 
    gaffect(gsub(gcoeff(B,4,2),gcoeff(B,4,1)), gmael3(A,2,2,2)); 

  }

  /* Restore B if necessary */ 
  if (thedet < 0) {
    B[4] = old_col4; 
  }

  avma = ltop; 

  A = gdiv(A, gsqrt(det2(A), prec));
  R[1] = (long)gerepileupto(ltop, A);
  R[2] = (long)stoi(thedet); 
  return R; 
}

/* GEN to_pari_O31() */
#if 0
int main()
{
  long d=30, prec = 3; 
  GEN x, y, a, b, m; 

  pari_init(1000000,2);
  // printf("precision of the computation in decimal digits: "); 
  // d = itos(lisGEN(stdin));
  if (d>0) prec = (long) (d*pariK1+3);

#if 0
  printf("first vector in GP format? "); 
  a = lisGEN(stdin); 
  printf("second vector in GP format? "); 
  b = lisGEN(stdin); 
#endif

  printf("matrix in GP format? "); 
  m = lisGEN(stdin);

  x = SL2C_to_O31(m, prec); 

  // x = O1n_GramSchmidt(m, prec);

  d = 10; 

  // x = gmul(O1n_inverse(x), x); 


  y = O31_to_SL2C(x, prec); 

  sor(m, 'g', d, 0); 

  sor(x, 'g', d, 0); 

  sor(y, 'g', d, 0); 

  printf("\n");
  return 0; 
}
#endif

#endif 

/* Keep compiler happy if NOPARI is #defined. */
int dummy_F100() { return 0; }
