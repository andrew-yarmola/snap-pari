AC_INIT(snap,1.11.3) 
AC_CONFIG_SRCDIR(src/snap/snap_prog.cc)
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS(config.h)
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

if test -d /sw; then
  AC_MSG_NOTICE([found /sw])
  LDFLAGS="$LDFLAGS -L/sw/lib"
  CFLAGS="$CFLAGS -I/sw/include"
  CXXFLAGS="$CXXFLAGS -I/sw/include"
fi

AC_ARG_WITH(readline-dir, AC_HELP_STRING([--with-readline-dir=DIR],
  [the directory under which READLINE is installed]), [
    RL_DIR="$with_readline_dir"
  ], [
    RL_DIR=
])

if test -n "$RL_DIR" -a "$RL_DIR" != /sw -a -d "$RL_DIR"; then
  AC_MSG_NOTICE([found "$RL_DIR"])
  LDFLAGS="$LDFLAGS -L$RL_DIR/lib"
  CFLAGS="$CFLAGS -I$RL_DIR/include"
  CXXFLAGS="$CXXFLAGS -I$RL_DIR/include"
fi

AC_CHECK_SIZEOF(long)
AC_SUBST(SIZEOF_LONG)

AC_LANG(C++)
OLD_LIBS="$LIBS"
LIBS="$LIBS -lreadline"
AC_MSG_CHECKING([readline (standalone)])
AC_LINK_IFELSE(
  AC_LANG_PROGRAM([[
#include <stdio.h>
#include <readline/readline.h>]],
  	          [[char* b=readline("");]]), 
  [ 
    AC_MSG_RESULT([yes])
  ], [
    AC_MSG_RESULT([no])
    LIBS="$LIBS -ltermcap"
    AC_MSG_CHECKING([readline with termcap])
    AC_LINK_IFELSE(
      AC_LANG_PROGRAM([[
#include <stdio.h>
#include <readline/readline.h>]],
  	              [[char* b=readline("");]]), 
      [
        AC_MSG_RESULT([yes])
      ], [
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([snap requires readline])
      ]
    )
  ]
)

AC_MSG_CHECKING([compliant string::compare])
AC_COMPILE_IFELSE(
  AC_LANG_PROGRAM([[
#include <string>]],
  	          [[[std::string s; char cs[2]; s.compare(0,1,cs);]]]), 
  [ 
    AC_MSG_RESULT([yes])
  ], [
    AC_MSG_RESULT([no])
    AC_DEFINE(OLD_GNU_COMPARE,1,
     [Older GNU libs have a non-standard string::compare syntax])
  ]
)


AC_ARG_WITH(pari-dir, AC_HELP_STRING([--with-pari-dir=DIR],
  [the directory under which PARI is installed]), [
    PARI_DIR="$with_pari_dir"
  ], [
    PARI_DIR=
])

AC_MSG_NOTICE([checking for PARI (paridecl.h)])
for pdir in $PARI_DIR /usr/local /usr /sw "$prefix"; do
  AC_CHECK_FILE([$pdir/include/pari/paridecl.h], [
    PARI_DIR="$pdir"
    break
  ])
done
if test -z "$PARI_DIR"; then
  AC_MSG_ERROR([[snap requires PARI]])
else
  AC_SUBST(PARI_DIR)
  CPPFLAGS="-I$PARI_DIR/include"
fi

AC_LANG(C++)
AC_MSG_CHECKING([pariinl.h])
FIX_PARI_HEADER=
AC_COMPILE_IFELSE([
  AC_LANG_SOURCE([[extern "C" { 
#include <pari/pari.h> 
} unsigned long overflow;]])],
  [AC_MSG_RESULT([ok])],
  [AC_MSG_RESULT([broken])
    FIX_PARI_HEADER="fix_pariinl"; 
])
AC_SUBST(FIX_PARI_HEADER)

PWR_DIR=pariwrap
AC_SUBST(PWR_DIR)

LIBS="$LIBS -lpari"
AC_LANG(C)
AC_MSG_CHECKING([FpM_ker matid gp_read_str])
AC_LINK_IFELSE(
   [AC_LANG_SOURCE([[#include <pari/pari.h>
int main() {
    GEN a,b;
    FpM_ker(a,b);
    matid(1);
    gp_read_str("[1,0;0,1]");
}]])],
   [AC_MSG_RESULT([found. Assuming pari 2.2 or later]);
    AC_DEFINE(PARI_2_2_OR_LATER,1,"pari version >= 2.2")],
   [AC_MSG_RESULT([not found. Assuming pari 2.1 or earlier])])

AC_LANG(C++)

AC_OUTPUT(Makefile \
	src/Makefile \
	src/pariwrap/Makefile \
	src/snappea/Makefile \
	src/tubes/Makefile \
	src/snap/Makefile)

