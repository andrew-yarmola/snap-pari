#ifndef NOPARI

#include "to_pari.hh"
#include "pariwrap.hh"

Complex to_complex(const pari& p)
{
  static pari big = pow(real(10),integer(10)); 
  if (gabs(p) > big) {
    return Infinity; 
  }

  Complex z;
  z.real = greal(p).double_value(); 
  z.imag = gimag(p).double_value(); 
  return z; 
}

void to_Moebius(const pari& p, MoebiusTransformation& m)
{
  m.matrix[0][0] = to_complex(p[0][0][0]); 
  m.matrix[1][0] = to_complex(p[0][0][1]); 
  m.matrix[0][1] = to_complex(p[0][1][0]); 
  m.matrix[1][1] = to_complex(p[0][1][1]); 
  if (p[1]==pONE) m.parity = orientation_preserving; 
  else m.parity = orientation_reversing; 
  m.acc_matrix = p[0];
}

pari to_pari(const MoebiusTransformation& m)
{
  pari r = rvector(2); 

  if (m.acc_matrix.type()!=t_INT) {
    r[0] = m.acc_matrix;
  } else {
    r[0] = matrix(2,2); 
    r[0][0][0] = to_pari(m.matrix[0][0]); 
    r[0][0][1] = to_pari(m.matrix[1][0]); 
    r[0][1][0] = to_pari(m.matrix[0][1]); 
    r[0][1][1] = to_pari(m.matrix[1][1]); 
  }
  if (m.parity == orientation_preserving) 
    r[1] = pONE; 
  else r[1] = -pONE; 

  return r;
}

#endif 

/* Keep compiler happy if NOPARI is #defined. */
int dummy_F102() { return 0; }
