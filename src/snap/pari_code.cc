#include "pari_code.hh"
#include "pariwrap.hh"

bool small(pari const& x)
{
  return pow(real(10),integer(digits/2))*abs(x) < pONE;
}

bool is_imag_quad_integer(pari const& p)
{
  pari re = pTWO*greal(p);
  pari reint = ground(re);
  if (!small(re-reint)) return false; 

  pari imsq; 
  if (gmod(reint,pTWO)==pONE) {
    imsq = (gsqr(pTWO*gimag(p))+pONE)/pFOUR;
  } else {
    imsq = gsqr(gimag(p)); 
  }
  return small(imsq-ground(imsq)); 
}

