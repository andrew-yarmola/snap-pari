#include <pari/pari.h>

GEN O1n_normsq(GEN x);
GEN O1n_innerprod(GEN u, GEN v);
GEN O1n_inverse(GEN y);
GEN O1n_GramSchmidt(GEN y, long prec);
GEN SL2C_to_O31(GEN m, long prec);
GEN O31_to_SL2C(GEN B, long prec);

