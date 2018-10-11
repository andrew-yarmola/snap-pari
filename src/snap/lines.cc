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
#include "lines.hh"
#include <stdio.h>

#define  PRINT_LINE_MATRIX 0
#define  PRINT_ORTHOLINES 0

// static say xx("lines.cc"); 

static void sl2c_minus(const SL2CMatrix a, const SL2CMatrix b, SL2CMatrix result);

line::line(const MoebiusTransformation& mt) 
{
  end[0] = mt.matrix[0][1]/mt.matrix[1][1]; 
  end[1] = mt.matrix[0][0]/mt.matrix[1][0];
}

void line::print() const 
{
  printf("["); ::print(end[0]); printf(","); ::print(end[1]); printf("]");
}

void reflection_in_line(const line& l, SL2CMatrix refl)
{
  if (!complex_big(l.end[0]) && !complex_big(l.end[1])) {
    Complex sum, prod, inv_det;
    sum = complex_plus(l.end[0],l.end[1]);
    prod = complex_mult(l.end[0],l.end[1]);
    
    inv_det = complex_div(I,complex_minus(l.end[0],l.end[1]));

    refl[0][0] = complex_mult(complex_minus(Zero, sum),inv_det);
    refl[0][1] = complex_mult(complex_mult(Two, prod),inv_det);
    refl[1][0] = complex_mult(complex_minus(Zero, Two),inv_det);
    refl[1][1] = complex_mult(sum,inv_det);
    
    return; 
  }

  int i = (complex_big(l.end[0]) ? 1 : 0);

  refl[0][0] = complex_mult(I,MinusOne);
  refl[0][1] = complex_mult(complex_mult(Two,I), l.end[i]);
  refl[1][0] = Zero;
  refl[1][1] = I;
}

line orthogonal_through_z_axis(const line& l)
{
  Complex e = complex_sqrt(l.end[0] * l.end[1]); 

  return line(-e, e); 
}

line orthogonal_through_origin(const line& l)
{
  // If l itself goes through the origin, there is no unique
  // perpendicular through the origin. We choose the
  // common orthogonal between l and the z axis. 
  if (same_point(l.end[0], -One/l.end[1]))
    return orthogonal_through_z_axis(l); 

  Complex b = (l.end[0] * l.end[1] - One) / (l.end[0] + l.end[1]); 

  Complex d = complex_sqrt(b * b + One);

  return line(b - d, b + d);
}

line ortholine(const line& line1, const line& line2)
{
  SL2CMatrix temp1, temp2, trans;
  line ol; 

  reflection_in_line(line1, temp1); 
  reflection_in_line(line2, temp2); 

  /* product of two reflections is translation and rotation 
     through double the distance and angle between the two lines */ 
  sl2c_product(temp2, temp1, trans); 

  fixed_points(trans, ol); 
  return ol; 
}

/* finds normalized line matrix for an oriented line,
    given the endpoints, as  in Fenchel, p.64. */
    
void find_line_matrix(
    const line& l,
    SL2CMatrix line_matrix)
{
    Complex factor;
#if PRINT_LINE_MATRIX
    FILE *fp;
    Complex fix[2];
    int     n_fp;
#endif
    
    if(!complex_infinite(l.end[0]) && !complex_infinite(l.end[1]))
    {
	factor = complex_div(I,complex_minus(l.end[1],l.end[0]));
	line_matrix[0][0] = complex_mult(factor,complex_plus(l.end[0],l.end[1]));
	line_matrix[0][1] = complex_mult(factor, 
				complex_real_mult(-2.0,complex_mult(l.end[0],l.end[1])));
	line_matrix[1][0] = complex_real_mult(2.0,factor);
	line_matrix[1][1] = complex_negate(line_matrix[0][0]);
#if PRINT_LINE_MATRIX
	fp = stdout;
	
	fprintf(fp, "%18.15lf %+18.15lfi %18.15lf %+18.15lfi \n",
	    line_matrix[0][0].real,line_matrix[0][0].imag,
	    line_matrix[0][1].real,line_matrix[0][1].imag);
	fprintf(fp, "%18.15lf %+18.15lfi %18.15lf %+18.15lfi \n",
	    line_matrix[1][0].real,line_matrix[1][0].imag,
	    line_matrix[1][1].real,line_matrix[1][1].imag);
	    fprintf(fp,"det = %lf %+lfi\n fp: ",sl2c_determinant(line_matrix));
	    n_fp= fixed_points(line_matrix, fix,stdout);
	    fprintf(fp,"\n"); 
#endif

	return;
    }
    if(!complex_infinite(l.end[0]) && complex_infinite(l.end[1]) )
    {
	line_matrix[0][0] = I;
	line_matrix[0][1] = complex_mult(I,complex_real_mult(-2.0,l.end[0]));
	line_matrix[1][0] = Zero;
	line_matrix[1][1] = complex_negate(line_matrix[0][0]);
#if PRINT_LINE_MATRIX
	fp = stdout;
	
	fprintf(fp, "%18.15lf %+18.15lfi %18.15lf %+18.15lfi \n",
	    line_matrix[0][0].real,line_matrix[0][0].imag,
	    line_matrix[0][1].real,line_matrix[0][1].imag);
	fprintf(fp, "%18.15lf %+18.15lfi %18.15lf %+18.15lfi \n",
	    line_matrix[1][0].real,line_matrix[1][0].imag,
	    line_matrix[1][1].real,line_matrix[1][1].imag);
	    fprintf(fp,"det = %lf %+lfi\n fp: ",sl2c_determinant(line_matrix));
	    n_fp= fixed_points(line_matrix, fix,stdout);            
	    fprintf(fp,"\n"); 
#endif

	return; 
    }
    if(!complex_infinite(l.end[1]) && complex_infinite(l.end[0]) )
    {
	line_matrix[1][1] = I;
	line_matrix[0][1] = complex_mult(I,complex_real_mult(2.0,l.end[1]));
	line_matrix[1][0] = Zero;
	line_matrix[0][0] = complex_negate(line_matrix[1][1]);
#if PRINT_LINE_MATRIX
	fp = stdout;
	
	fprintf(fp, "%18.15lf %+18.15lfi %18.15lf %+18.15lfi \n",
	    line_matrix[0][0].real,line_matrix[0][0].imag,
	    line_matrix[0][1].real,line_matrix[0][1].imag);
	fprintf(fp, "%18.15lf %+18.15lfi %18.15lf %+18.15lfi \n",
	    line_matrix[1][0].real,line_matrix[1][0].imag,
	    line_matrix[1][1].real,line_matrix[1][1].imag);
	    fprintf(fp,"det = %lf %+lfi\n fp: ",sl2c_determinant(line_matrix));
	    n_fp= fixed_points(line_matrix, fix,stdout);
	    fprintf(fp,"\n"); 
#endif

	return; 
	
    }
    printf("Big trouble in find_line_matrix. Both endpoints are infinite!");
    return;
}


	/*
	 *  Finds the "ortholine" perpendicular to two lines given by
	 *      axes of matrices. (The matrices don't need to be line matrices here!)
	 */

void find_ortholine(const SL2CMatrix line1, const SL2CMatrix line2, SL2CMatrix ortholine)
{
    SL2CMatrix temp1, temp2;
#if PRINT_ORTHOLINES    
    FILE *fp;
    Complex fix[2];
    int     n_fp;
#endif
	/*
	 *  The ortholine determined by two translation matrices
	 *  f and g is gf - fg  (see Fenchel, p. 62).  The sign
	 *  of the matrix may or may not be correct.
	 */

	sl2c_product(line2, line1, temp1);
	sl2c_product(line1, line2, temp2);
	sl2c_minus(temp1, temp2, ortholine);
		/* if det = 0, should print error message  */
	sl2c_normalize(ortholine);
#if PRINT_ORTHOLINES
	fp =stdout;
	n_fp= fixed_points(ortholine, fix,fp);
#endif

}

/* gives complex distance from line 1 to line 2 */

Complex orthodist(const line& line1, const line& line2)
{
  SL2CMatrix temp1, temp2, trans;
  line ol; 

  reflection_in_line(line1, temp1); 
  reflection_in_line(line2, temp2); 

  /* product of two reflections is translation and rotation 
     through double the distance and angle between the two lines */ 
  sl2c_product(temp2, temp1, trans); 

  fixed_points(trans, ol); 

  return cross_ratio(ol.end[0], ol.end[1], line2.end[1], line1.end[1]); 
}

Complex orthodist(const SL2CMatrix line1, const SL2CMatrix line2)
{
    SL2CMatrix  ml;
    Complex     tr,exp_mu, orthodist;
	/*
	 *  We use the formula from page 68 of Fenchel:
	 *
	 *      cosh(mu) = (-1/2) tr(ml)
	 *
	 * where l = normalized line matrix for line1, 
	 *       m = normalized line matrix for line 2.
	 *
	 */

    /* The formula gives:
	exp(mu) + exp(-mu) = - tr(ml) or exp(mu)^2 + tr*exp(mu) + 1 = 0, 
	 where tr = tr(ml). Solving this quadratic gives:
	    exp(mu) = ( -tr +- sqrt(tr^2-4) )/2.
	The two solutions are reciprocals, so give opposite values for the
	complex distance mu. We can take either one, and then adjust mu to 
	have positive real part.
	
    */
    
	    sl2c_product(line2, line1, ml);
	    tr = complex_plus( ml[0][0],ml[1][1]);
	    exp_mu =complex_real_mult(0.5, 
			    complex_minus(
				complex_sqrt(complex_minus(complex_mult(tr,tr),Four)),
						tr));
	    orthodist = complex_log(exp_mu,0.0); /* assume approx argument is zero */
	    if(orthodist.real < 0.0)
		orthodist = complex_negate(orthodist);
	    return(orthodist);
}

#if 0
Complex orthodist_from_ends(Complex axis1[2],Complex axis2[2])
{
    SL2CMatrix   line1,line2;
    
    find_line_matrix(axis1,line1);
    find_line_matrix(axis2,line2);
    return  orthodist(line1,line2);
}
#endif

Complex cosh_orthodist(const line& axis1, const line& axis2)
{
    SL2CMatrix line1, line2;
    
    find_line_matrix(axis1,line1);
    find_line_matrix(axis2,line2);
    return  cosh_orthodist(line1,line2);
}

void find_ortholine(const line& axis1, const line& axis2, SL2CMatrix ortholine_matrix)
{
    SL2CMatrix line1, line2;
    
    find_line_matrix(axis1,line1);
    find_line_matrix(axis2,line2);
    find_ortholine(line1, line2, ortholine_matrix);
}




/* Find cosh(orthodistance) from line1 to line 2 */

Complex cosh_orthodist(const SL2CMatrix line1, const SL2CMatrix line2)
{
    SL2CMatrix  ml;
    Complex     tr;
	/*
	 *  We use the formula on page 68 of Fenchel:
	 *
	 *      cosh(mu) = (-1/2) tr(ml)
	 *
	 * where l = normalized line matrix for line1, 
	 *       m = normalized line matrix for line 2.
	 *
	 */

	    sl2c_product(line2, line1, ml);
	    tr = complex_plus( ml[0][0],ml[1][1]);
	    return complex_real_mult(-0.5, tr);
}


static void sl2c_minus(
    const SL2CMatrix a,
    const SL2CMatrix b,
    SL2CMatrix result)
{
    int i,
	j;

    for (i=0; i<2; i++)
	for (j=0; j<2; j++)
	    result[i][j] = complex_minus(a[i][j], b[i][j]);

    return;
}

/* returns the fixed point  (a-d)/2c 
for a  parabolic SL2C matrix m = {{a,b},{c,d}}  */

Complex parabolic_fix(const SL2CMatrix m)
{
  if(complex_nonzero(m[1][0]))
    return(complex_div(complex_minus(m[0][0],m[1][1]),
		       complex_real_mult(2.0,m[1][0])));
  else
    return(Infinity);
}

int fixed_points(const SL2CMatrix m, line& fix, FILE *fp)
{
  Complex a,b,c,d, diff, discr, tr;

  a = m[0][0]; b = m[0][1]; 
  c = m[1][0]; d = m[1][1]; 

  tr = a + d; 
  discr = tr * tr - Four; 
  diff = a - d; 

  if(complex_small(c))  /* if c = 0 */ 
    {
      if(complex_small(diff))    /* a - d = 0 */
	{
	  if(complex_small(b))  /* b = 0, implies identity*/
	    { 
	      if (fp) fprintf(fp,"All points are fixed! (Identity matrix) \n");
	      return -1;
	    }
	  else    /* b != 0 implies 1 fixed point at infinity */
	    {
	      fix.end[0] = fix.end[1] = Infinity;
	      if(fp) fprintf(fp,"  Infinity    Infinity\n");
	      return 1;
	    }
	}
      else   /* a - d != 0, two fixed points: b/(d-a), infinity */
	{
	  fix.end[0] = b/(d-a); 
	  fix.end[1] = Infinity;
	  if(fp)
	    fprintf(fp,"%9.6lf %+9.6lfi    Infinity\n",
		    fix.end[0].real, fix.end[0].imag);
	}
    }
  else /* c!= 0 */
    {
      if(complex_small(discr))  
	{ /*   discr = 0 implies there is one fixed point = (a-d)/2c   */
	  fix.end[0] = fix.end[1] = (a-d)/(2.0 * c); 
	  if(fp)
	    fprintf(fp,"%9.6lf %+9.6lfi  %9.6lf %+9.6lfi\n",
		    fix.end[0].real,fix.end[0].imag,fix.end[1].real,fix.end[1].imag);
	  return 1;
	}
      else /* discr != 0, implies 2 finite fixed points: ( a-d +/- sqrt(discr) )/2c  */
	{
	  discr =  complex_sqrt(discr);
	  fix.end[0] = (diff + discr)/(2.0 * c); 
	  fix.end[1] = (diff - discr)/(2.0 * c); 
	  
	  if(fp) fprintf(fp,"%9.6lf %+9.6lfi  %9.6lf %+9.6lfi\n", fix.end[0].real,fix.end[0].imag,fix.end[1].real,fix.end[1].imag);
	}
    }
    
  /* put the fixed points into canonical order, repelling, attracting.
     derivative is 1/(cz+d)^2. we consider (c fix.end[0] + d) */ 
  Complex inv_sqrt_deriv = c * fix.end[0] + d; 

  if (complex_modulus_squared(inv_sqrt_deriv) - 1.0 < -1e-10) 
    return 2; /* fix.end[0] repelling */
 
  if (complex_modulus_squared(inv_sqrt_deriv) - 1.0 < 1e-10) { 
    /* pure rotation */
    Complex log_inv_deriv = complex_log(inv_sqrt_deriv, 0.0);
    if (log_inv_deriv.imag < 0.0) 
      return 2; /* positive rotation at fix.end[0] */ 
  }

  /* if we got here fix.end[0] was attracting or a point of pure 
     negative rotation */ 
  Complex tmp = fix.end[0]; 
  fix.end[0] = fix.end[1];
  fix.end[1] = tmp;
  return 2; 
}   



int  real_fixed_points(const SL2CMatrix m, line& fix, FILE *fp)
{
    Complex diff, discr, tr;
    
    tr = complex_plus(m[0][0],m[1][1]); /* trace = a+d */
    discr = complex_minus(complex_mult(tr,tr),Four); /* discr = tr^2 - 4 */
    diff = complex_minus(m[0][0],m[1][1]); /* diff =  a-d */
    if(complex_modulus_squared(m[1][0]) < 1e-30) /* if c = 0 */                         
    {
	if(complex_modulus_squared(diff) < 1e-30)    /* a - d = 0 */
		{
		if( complex_modulus_squared(diff) < 1e-30)  /* b = 0, implies identity*/
		    {
			if(fp)
			    fprintf(fp,"All points are fixed! (Identity matrix) \n");
			return(-1);

		    }
		else    /* b != 0 implies 1 fixed point at infinity */
		    {
			fix.end[0] = fix.end[1] = Infinity;
			if(fp)
			    fprintf(fp,"  Infinity                                Infinity                              \n");
			return(0);
		    }
		    
		}
	    else   /* a - d != 0, two fixed points: b/(d-a), infinity */
		{
		    diff = complex_negate(diff);
		    fix.end[0] = complex_div(m[0][1],diff);
		    fix.end[1] = Infinity;
		    if(fp)
			fprintf(fp,"%9.6lf    Infinity                              \n",
					    fix.end[0].real);
		    return(1);
		}
    
    }
    
    else /* c!= 0 */
    {
	if(complex_modulus_squared(discr) < 1e-30)  
	    { /*   discr = 0 implies there is one fixed point = (a-d)/2c   */
		fix.end[0]=fix.end[1] = complex_div(diff,complex_real_mult(2.0,m[1][0]));
		if(fp)
		    fprintf(fp,"%9.6lf  %9.6lf\n",
				fix.end[0].real,fix.end[1].real);
		return(1);
	    }
	else   /* discr != 0, implies 2 finite fixed points: ( a-d +/- sqrt(discr) )/2c  */
	    {
		discr =  complex_sqrt(discr);
		fix.end[0] = complex_div( complex_plus(diff,discr),complex_mult(m[1][0],Two));
		fix.end[1] = complex_div( complex_minus(diff,discr),complex_mult(m[1][0],Two));
		if(fp)
		    fprintf(fp,"%9.6lf  %9.6lf\n",
			fix.end[0].real,fix.end[1].real);
		return(2);
	    
	    }
	    
    }
	
}


line operator * (const MoebiusTransformation& m, const line& l)
{
  line ml;

  ml.end[0] = m * l.end[0];
  ml.end[1] = m * l.end[1];

  return ml;
}

bool operator == (const line& a, const line& b)
{
  return same_point(a.end[0],b.end[0]) && same_point(a.end[1], b.end[1]); 
}

bool operator < (const line& a, const line& b)
{
  return a.end[0] < b.end[0] || (same_point(a.end[0],b.end[0]) && a.end[1] < b.end[1]);
}

int compare(const line& a, const line& b, double eps)
{
  if (same_point(a.end[0],b.end[0],eps))
    return same_point(a.end[1],b.end[1],eps) ? 1 : 0; 
  else
    if (same_point(a.end[0],b.end[1],eps) && same_point(a.end[1],b.end[0],eps)) return -1; 
  
  return 0; 
}
