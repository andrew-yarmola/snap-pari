#ifndef _pari_oldnames_
#define _pari_oldnames_

/* The functions below were renamed in pari 2.2 but are still used
   throughout snap.
   Include this file like this
   #ifdef PARI_2_2_OR_LATER
   #include "pari_oldnames.hh"
   #endif
   so that it compiles with pari < 2.2 and pari >= 2.2. */

#define ker_mod_p FpM_ker

#define p_idmat p_matid

#define p_lisexpr(x) gp_read_str((char *) (x))
#define p_flisexpr(x) gp_read_str((char *) (x))

#ifndef pariK
#define pariK  (9.632959862*(BYTES_IN_LONG/4))  /* SL*log(2)/log(10) */
#endif

#endif
