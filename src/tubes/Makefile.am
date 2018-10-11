bin_PROGRAMS = tube 
# drill check_face isometric tubal
# try_eqs
# eqs_test pictest find_classes find_cusp_density

shared = \
	../snap/record.cc			\
	../snap/record.hh			\
	../snap/printable.hh			\
	../snap/printable.cc			\
	../snap/file_array.hh			\
	../snap/snap_io.cc			\
	../snap/snap_io.hh			\
	../snap/terse_io.cc			\
	../snap/terse_io.hh			\
	../snap/linktable.cc 			\
	../snap/linktable.hh			\
	../snap/menu.cc				\
	../snap/menu.hh				\
	../snap/n_vector.cc			\
	../snap/n_vector.hh			\
	../snap/rc_matrix.cc			\
	../snap/rc_matrix.hh			\
	../snap/VS.hh				\
	../snap/VS.cc				\
	../snap/polytope_T.hh			\
	../snap/tilt_polytope.cc		\
	../snap/tilt_polytope.hh		\
	../snap/commens.cc			\
	../snap/commens.hh			\
	../snap/find_commens.cc			\
	../snap/find_commens.hh			\
	../snap/helpers.cc			\
	../snap/helpers.hh			\
	../snap/i_triangulation.cc		\
	../snap/i_triangulation.hh		\
	../snap/drill_Dirichlet.cc		\
	../snap/drill_Dirichlet.hh		\
	../snap/env_T.cc			\
	../snap/env_T.hh			\
	../snap/env_D.cc			\
	../snap/env_D.hh			\
	../snap/env_U.cc			\
	../snap/env_U.hh

dd_tiling = \
	../snap/color.cc			\
	../snap/color.hh			\
	../snap/kernel_extras.cc		\
	../snap/kernel_extras.hh		\
	../snap/int_matrix.cc			\
	../snap/int_matrix.hh			\
	../snap/lines.cc			\
	../snap/lines.hh			\
	../snap/O31_line.cc			\
	../snap/O31_line.hh			\
	../snap/closed_geodesics.cc		\
	../snap/closed_geodesics.hh		\
	../snap/ortholine.cc			\
	../snap/ortholine.hh			\
	../snap/dd_tiling.cc			\
	../snap/dd_tiling.hh			\
	../snap/gv_print_Dirichlet.cc		\
	../snap/gv_print_Dirichlet.hh

tube_code = \
	../snap/mark.cc				\
	../snap/mark.hh				\
	../snap/R_vector.hh			\
	../snap/R_matrix.hh			\
	../snap/R_poly.cc			\
	../snap/R_poly.hh			\
	../snap/adaptive_curve.hh		\
	../snap/eqs_interval.cc			\
	../snap/eqs_interval.hh			\
	../snap/vfunction.hh			\
	../snap/point.cc			\
	../snap/point.hh			\
	../snap/complex_lll.cc			\
	../snap/complex_lll.hh			\
	../snap/eq_surface.cc			\
	../snap/eq_surface.hh			\
	../snap/bbox.cc				\
	../snap/bbox.hh				\
	../snap/picfile.cc			\
	../snap/picfile.hh			\
	../snap/polylist_base.hh		\
	../snap/uc_polylist.cc			\
	../snap/uc_polylist.hh			\
	../snap/tube_face.cc			\
	../snap/tube_face.hh			\
	../snap/connected_faces.cc		\
	../snap/connected_faces.hh		\
	../snap/triangulation_builder.cc	\
	../snap/triangulation_builder.hh	\
	../snap/tile.cc				\
	../snap/tile.hh				\
	../snap/orthoangles.cc			\
	../snap/orthoangles.hh			\
	../snap/tube.cc				\
	../snap/tube.hh	

tube_only = warn.cc warn.hh

tube_SOURCES = $(shared) $(dd_tiling) $(tube_code) $(tube_only) tube_prog.cc
tube_LDADD = ../snappea/libsnap.a -lreadline -lm

INCLUDES = -I$(srcdir)/.. -I$(srcdir)/../snap
AM_CPPFLAGS = -DNOPARI -DDATADIR='"$(datadir)/snap_data"'

# find_classes_SOURCES = $(shared) $(tube_only) find_classes.cc
# find_classes_LDADD = ../snappea/libsnap.a -lreadline -lm

# find_cusp_density_SOURCES = $(shared) $(tube_only) find_cusp_density.cc
# find_cusp_density_LDADD = ../snappea/libsnap.a -lreadline -lm

# FRAGS = headerfrag.hh menufrag.cc statefrag.cc switchfrag.cc

# BUILT_SOURCES = $(FRAGS)

# $(FRAGS) test_help: eqs_test.cmds
# 	./cmds2frag $<

# eqs_test_SOURCES = eqs_test.cmds eqs_test.cc \
# 	warn.cc $(shared)
# eqs_test_LDADD = ../snappea/libsnap.a -lm

EQSFILES = ../snap/R_poly.cc ../snap/mark.cc ../snap/point.cc \
	../snap/printable.cc ../snap/bbox.cc \
	../snap/eq_surface.cc ../snap/complex_lll.cc ../snap/eqs_interval.cc

# try_eqs_SOURCES = $(EQSFILES) try_eqs.cc
# try_eqs_LDADD = ../snappea/libsnap.a -lm

# pictest_SOURCES = pictest.cc ../snap/picfile.cc ../snap/color.cc
# pictest_LDADD = ../snappea/libsnap.a -lm

# check_face_SOURCES = $(EQSFILES) check_face.cc
# check_face_LDADD = ../snappea/libsnap.a -lm

# drill_SOURCES = $(dd_tiling) $(tube_code) $(tube_only) 	\
# 	../snap/snap_io.cc 			\
# 	../snap/terse_io.cc			\
# 	../snap/linktable.cc			\
# 	../snap/printable.cc			\
# 	../snap/drill.cc 

# drill_LDADD = ../snappea/libsnap.a -lm

# isometric_SOURCES = $(dd_tiling) warn.cc ../snap/printable.cc \
# 	../snap/snap_io.cc ../snap/isometric.cc

# isometric_LDADD = ../snappea/libsnap.a -lm

# tubal_SOURCES = $(dd_tiling) $(tube_code) $(tube_only) \
# 	../snap/snap_io.cc \
# 	../snap/printable.cc \
# 	../snap/tubal.cc

# tubal_LDADD = ../snappea/libsnap.a -lm

CLEANFILES = tube_menu

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(datadir)/snap_data
	echo "save menu tube_menu" | ./tube
	$(INSTALL_DATA) tube_menu $(DESTDIR)$(datadir)/snap_data/tube_menu

uninstall-local:
	rm -f $(DESTDIR)$(datadir)/snap_data/tube_menu

