// -*-C++-*-
#ifndef _to_pari_
#define _to_pari_

#include "snappea/Moebius_transformations.h"
#include "pariclass.hh"

inline pari to_pari(const Complex& z)
{ return complex(pari(z.real), pari(z.imag)); }
Complex to_complex(const pari& p);

pari to_pari(const MoebiusTransformation& m);
void to_Moebius(const pari& p, MoebiusTransformation& m);

#endif
