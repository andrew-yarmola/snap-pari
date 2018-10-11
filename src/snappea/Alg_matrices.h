#ifndef _Alg_matrices_
#define _Alg_matrices_

#include "SnapPea.h"

struct CyclicWord;

struct AlgMatrix
{
  O31_matrix matrix;
  CyclicWord *name;
  Boolean    identity;
};


extern AlgMatrix    Alg_identity;

extern void         Alg_copy(AlgMatrix **dest, AlgMatrix *source);
extern void         Alg_invert(AlgMatrix *m, AlgMatrix **m_inverse);
extern void         Alg_product(AlgMatrix *a, AlgMatrix *b, AlgMatrix **product);
extern Boolean      Alg_equal(AlgMatrix *a, AlgMatrix *b, double epsilon);
extern double       Alg_deviation(AlgMatrix *m);
extern void         Alg_GramSchmidt(AlgMatrix *m);
extern void         Alg_conjugate(AlgMatrix *m, AlgMatrix *t, AlgMatrix **Tmt);

extern void init_name(CyclicWord **);
extern void simplify_word(CyclicWord *);

//  extern void                     append_inverse(CyclicWord  *,CyclicWord *);
//  extern void                     append_word (CyclicWord *,CyclicWord *);
// extern void                     cancel_inverses_word(CyclicWord *);
// extern void                     free_cyclic_word(CyclicWord *);
#endif


