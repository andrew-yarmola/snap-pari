#ifndef _pari_code_hh_
#define _pari_code_hh_

#include "pariclass.hh"

bool small(pari const& x);
bool is_imag_quad_integer(pari const& p);
pari ramification(const pari hilbert[2], pari itf, int int_traces, int report, int max_time); 

const pari_printform no_newline = pari_printform(1,0,'g',-1,0);

#endif
