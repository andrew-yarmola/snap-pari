bin_PROGRAMS = snap make_linkdata
# tpt refl closed_fields
#  voronoi

shared = \
	color.cc			\
	color.hh			\
	bbox.cc				\
	bbox.hh				\
	picfile.cc			\
	picfile.hh			\
	record.cc			\
	record.hh			\
	vecit.hh			\
	file_array.hh			\
	printable.hh			\
	printable.cc			\
	linktable.cc 			\
	linktable.hh			\
	menu.cc				\
	menu.hh				\
	n_vector.cc			\
	n_vector.hh			\
	rc_matrix.cc			\
	rc_matrix.hh			\
	VS.hh				\
	VS.cc				\
	polytope_T.hh			\
	R_vector.hh			\
	R_matrix.hh			\
	R_poly.cc			\
	R_poly.hh			\
	int_matrix.cc			\
	int_matrix.hh			\
	complex_lll.cc			\
	complex_lll.hh			\
	tilt_polytope.cc		\
	tilt_polytope.hh		\
	commens.cc			\
	commens.hh			\
	find_commens.cc			\
	find_commens.hh			\
	kernel_extras.cc		\
	kernel_extras.hh		\
	terse_io.cc			\
	terse_io.hh			\
	helpers.cc			\
	helpers.hh			\
	snap_io.cc			\
	snap_io.hh			\
	tile.cc				\
	tile.hh				\
	orthoangles.cc			\
	orthoangles.hh			\
	i_triangulation.cc		\
	i_triangulation.hh		\
	lines.cc			\
	lines.hh			\
	O31_line.cc			\
	O31_line.hh			\
	closed_geodesics.cc		\
	closed_geodesics.hh		\
	ortholine.cc			\
	ortholine.hh			\
	dd_tiling.cc			\
	dd_tiling.hh			\
	gv_print_Dirichlet.cc		\
	gv_print_Dirichlet.hh		\
	drill_Dirichlet.cc		\
	drill_Dirichlet.hh		\
	eq_surface.cc			\
	eq_surface.hh			\
	adaptive_curve.hh		\
	mark.cc				\
	mark.hh				\
	eqs_interval.cc			\
	eqs_interval.hh			\
	vfunction.hh			\
	point.cc			\
	point.hh			\
	triangulation_builder.cc	\
	triangulation_builder.hh	\
	tube_face.cc			\
	tube_face.hh			\
	uc_polylist.cc			\
	uc_polylist.hh			\
	tube.cc				\
	tube.hh				\
	connected_faces.cc		\
	connected_faces.hh

with_pari = \
	pari_ext.c \
	pari_ext.h \
	O31_GEN.c \
	O31_GEN.h \
	O31_pari.hh \
	to_pari.cc \
	to_pari.hh \
	p_VS.cc \
	p_VS.hh \
	pari_code.cc \
	pari_code.hh \
	get_equations.cc \
	get_equations.hh \
	field.cc \
	field.hh \
	ramification.cc \
	snap.cc \
	snap.hh \
	polish_rep.cc \
	polish_rep.hh \
	dirichlet.cc \
	dirichlet.hh \
	algebraic.cc \
	algebraic.hh \
	alg_group.cc \
	alg_group.hh

snap_interface = \
	env_T.cc	\
	env_T.hh	\
	env_D.cc	\
	env_D.hh	\
	env_U.cc	\
	env_U.hh	\
	env_A.cc 	\
	env_A.hh  	\
	snap_prog.cc 

snap_SOURCES = $(shared) $(with_pari) $(snap_interface)

make_linkdata_SOURCES = make_linkdata.cc

INCLUDES = -I.. -I@PARI_DIR@/include -I$(srcdir)/../@PWR_DIR@ \
	-I$(srcdir)/..
AM_LDFLAGS = -L@PARI_DIR@/lib
AM_CPPFLAGS = -DDATADIR='"$(datadir)/snap_data"'
snap_LDADD = ../snappea/libpsnap.a ../@PWR_DIR@/libpwr.a -lpari -lm

CLEANFILES = snap_menu

linkdata = snap_data/links_by_volume

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(datadir)/snap_data
	echo "save menu snap_menu" | ./snap -b
	$(INSTALL_DATA) snap_menu $(DESTDIR)$(datadir)/snap_data/snap_menu
	cat $(srcdir)/../../snap_data/link_data/hyp*[an] | grep '\.' \
	| sort -n +1 | ./make_linkdata
	$(INSTALL_DATA) links_by_volume $(DESTDIR)$(datadir)/$(linkdata)
	rm -f snap_menu links_by_volume

uninstall-local:
	rm -f $(DESTDIR)$(datadir)/snap_data/snap_menu
	rm -f $(DESTDIR)$(datadir)/$(linkdata)

# voronoi_SOURCES = n_vector.cc rc_matrix.cc polytope.cc voronoi.cc

# tpt_SOURCES = n_vector.cc rc_matrix.cc p_VS.cc test_polytope.cc
# tpt_LDADD = ../@PWR_DIR@/libpwr.a -lpari -lm

# refl_SOURCES = refl.cc p_VS.cc
# refl_LDADD = ../snappea/libpsnap.a ../@PWR_DIR@/libpwr.a -lpari -lm

# closed_fields_SOURCES = $(shared) $(with_pari) \
# 	db.cc db.hh \
# 	f_db.cc f_db.hh \
# 	a_db.cc a_db.hh \
# 	closed_fields.cc

# closed_fields_LDADD = ../snappea/libpsnap.a ../@PWR_DIR@/libpwr.a -lpari -lm

