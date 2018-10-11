/*
 *  Alg_matrices.c
 *
 *  This file provides the following functions for working with AlgMatrices:
 *
 *      void        Alg_copy(AlgMatrix *dest, AlgMatrix *source);
 *      void        Alg_invert(AlgMatrix *m, AlgMatrix *m_inverse);
 *      FuncResult  gl4R_invert(GL4RMatrix m, GL4RMatrix m_inverse);
 *      double      gl4R_determinant(GL4RMatrix m);
 *      void        Alg_product(AlgMatrix *a, AlgMatrix *b, AlgMatrix *product);
 *      Boolean     Alg_equal(AlgMatrix* a, AlgMatrix* b, double epsilon);
 *      double      Alg_trace(AlgMatrix* m);
 *      double      Alg_deviation(AlgMatrix* m);
 *      void        Alg_GramSchmidt(AlgMatrix *m);
 *      void        Alg_conjugate(AlgMatrix *m, AlgMatrix *t, AlgMatrix *Tmt);
 */

#include "Alg_matrices.h"
#include "kernel.h"
#include "fundamental_group.h"

/*
 *  gl4R_invert will consider a matrix to be singular iff one of the
 *  absolute value of one of the pivots is less than SINGULAR_MATRIX_EPSILON.
 */
#define SINGULAR_MATRIX_EPSILON     double(1e-2)

#define COLUMN_PRODUCT(m, i, j)     \
    (-m[0][i]*m[0][j] + m[1][i]*m[1][j] + m[2][i]*m[2][j] + m[3][i]*m[3][j])



void Alg_copy(
    AlgMatrix   **dest,
    AlgMatrix   *source)
{
	*dest = NEW_STRUCT(AlgMatrix);
	(*dest)->matrix = source->matrix; 
    if(source->identity == FALSE)
      { init_name (&((*dest)->name));
	append_word(source->name,(*dest)->name);
	simplify_word(source->name);
	simplify_word((*dest)->name);
	 (*dest)->identity = FALSE;
      }
    else
      (*dest)->identity = TRUE;
}


void Alg_invert(
    AlgMatrix*  m,
    AlgMatrix** m_inverse)
{
    *m_inverse =  NEW_STRUCT(AlgMatrix) ;
    (*m_inverse)->matrix = inverse(m->matrix);
    if(m->identity == FALSE)
      { init_name (&((*m_inverse)->name));
	simplify_word(m->name);
	 append_inverse(m->name,(*m_inverse)->name);
	(*m_inverse)->identity = FALSE;
      }
    else
      (*m_inverse)->identity =  TRUE;

}



void Alg_product(
    AlgMatrix*  a,
	 AlgMatrix* b,
	 AlgMatrix**    product)
{
  *product = NEW_STRUCT(AlgMatrix);
  (*product)->matrix = a->matrix * b->matrix; 
   if(a->identity ==  FALSE && b->identity == FALSE)
    { init_name(&((*product)->name));
      simplify_word(a->name);
      simplify_word(b->name);
      append_word(a->name,(*product)->name);
      append_word(b->name,(*product)->name);
      (*product)->identity = FALSE;
    }
  else
    {
      if (a->identity ==  FALSE)
    {init_name(&((*product)->name));
     simplify_word(a->name);
     append_word(a->name,(*product)->name);
     (*product)->identity = FALSE;
       }
      else if(b->identity == FALSE)
    {init_name(&((*product)->name));
     simplify_word(b->name);
     append_word(b->name,(*product)->name);
     (*product)->identity = FALSE;
       }
      else
    (*product)->identity = TRUE;
    }
}

Boolean Alg_equal(
    AlgMatrix*  a,
    AlgMatrix*  b,
    double      epsilon)
{
	return close(a->matrix,b->matrix,epsilon);
}


double Alg_trace(
    AlgMatrix*  m)
{
    return (trace(m->matrix));
}


double Alg_deviation(
    AlgMatrix*  m)
{
    return m->matrix.deviation();
}


void Alg_GramSchmidt(
    AlgMatrix*  m)
{
    m->matrix.O31_GramSchmidt();
}


void Alg_conjugate(
    AlgMatrix*  m,
    AlgMatrix*  t,
    AlgMatrix** Tmt)
{
    
    AlgMatrix   **t_inverse,
		**temp;

    Alg_invert(t, t_inverse);
    Alg_product(*t_inverse, m, temp);
    Alg_product(*temp, t, Tmt);
}

