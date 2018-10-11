INCLUDES = -I.. -I@PARI_DIR@/include
AM_LDFLAGS = -L@PARI_DIR@/lib

noinst_LIBRARIES = libpwr.a
noinst_PROGRAMS = example

inlines.hh: @PARI_DIR@/include/pari/paridecl.h @FIX_PARI_HEADER@
	$(top_srcdir)/src/pariwrap/make_inlines $< > $@

fix_pariinl:
	cp -r @PARI_DIR@/include/pari ..
	cat @PARI_DIR@/include/pari/pariinl.h |\
	sed -e 's/^ulong over/extern ulong over/' |\
	sed -e 's/^ulong hire/extern ulong hire/' > ../pari/pariinl.h
	echo "" > fix_pariinl

CLEANFILES = ../pari/* fix_pariinl inlines.hh unused_decls

BUILT_SOURCES = inlines.hh
EXTRA_DIST = make_inlines

libpwr_a_SOURCES = inlines.hh pariwrap.cc pariwrap.hh pariclass.hh \
	warn.cc warn.hh 

example_SOURCES = example.cc
example_LDADD = libpwr.a -lpari -lm
