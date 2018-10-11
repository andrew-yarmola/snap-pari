// -*-C++-*-
#ifndef _O31_pari_
#define _O31_pari_

#include "pariwrap.hh"

extern "C" {
#include "O31_GEN.h"
}

inline pari O1n_normsq(pari const& x)
{ return O1n_normsq(pari::tmp(x).g); }
inline pari O1n_innerprod(pari const& u, pari const& v)
{ return O1n_innerprod(pari::tmp(u).g, pari::tmp(v).g); }
inline pari O1n_inverse(pari const& y)
{ return O1n_inverse(pari::tmp(y).g); }
inline pari O1n_GramSchmidt(pari const& y, long prc = prec)
{ return O1n_GramSchmidt(pari::tmp(y).g, prc); }
inline pari SL2C_to_O31(pari const& m, long prc = prec)
{ return SL2C_to_O31(pari::tmp(m).g, prc); }
inline pari O31_to_SL2C(pari const& B, long prc = prec)
{ return O31_to_SL2C(pari::tmp(B).g, prc); }

#endif
