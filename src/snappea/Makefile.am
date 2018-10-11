noinst_LIBRARIES = libsnap.a libpsnap.a

snappea = Dehn_coefficients.cc \
	Dirichlet.cc Dirichlet.h \
	Dirichlet_basepoint.cc \
	Dirichlet_construction.cc \
	Dirichlet_extras.cc \
	Dirichlet_precision.cc \
	FGWord.cc FGWord.h \
	Moebius_transformations.cc Moebius_transformations.h \
	O31.cc O31.h \
	SnapPea.h \
	abelian_group.cc \
	barycentric_subdivision.cc \
	canonize.cc canonize.h \
	canonize_part_1.cc \
	canonize_part_2.cc \
	canonize_result.cc \
	change_peripheral_curves.cc \
	chern_simons.cc \
	choose_generators.cc \
	close_cusps.cc \
	complex.cc complex.h \
	complex_length.cc \
	continued_fractions.cc \
	core_geodesics.cc \
	cover.cc covers.h \
	current_curve_basis.cc \
	cusp_cross_sections.cc \
	cusp_neighborhoods.cc \
	cusp_shapes.cc \
	cusps.cc \
	decode_new_DT.cc \
	diagnostic.cc \
	direct_product.cc \
	double_cover.cc \
	drilling.cc \
	dual_curves.cc \
	dual_one_skeleton_curve.h \
	edge_classes.cc \
	elements_generate_group.cc \
	filling.cc \
	find_cusp.cc \
	finite_vertices.cc \
	fundamental_group.cc fundamental_group.h \
	gcd.cc \
	gluing_equations.cc \
	holonomy.cc \
	homology.cc \
	hyperbolic_structure.cc \
	identify_solution_type.cc \
	index_to_hue.cc \
	interface.cc \
	intersection_numbers.cc \
	isometry.cc isometry.h \
	isometry_closed.cc \
	isometry_cusped.cc \
	kernel.h \
	kernel_prototypes.h \
	kernel_typedefs.h \
	length_spectrum.cc length_spectrum.h \
	link_complement.cc \
	link_projection.h \
	matrix_conversion.cc matrix_conversion.h \
	matrix_generators.cc \
	my_malloc.cc \
	normal_surface_construction.cc \
	normal_surface_recognition.cc \
	normal_surface_splitting.cc \
	normal_surfaces.h \
	o31_matrices.cc o31_matrices.h \
	orient.cc \
	peripheral_curves.cc peripheral_curves.h \
	polyhedral_group.cc \
	positioned_tet.cc positioned_tet.h \
	precision.cc \
	punctured_torus_bundles.cc \
	quotients.cc \
	rehydrate_census.cc \
	representations.cc \
	shingling.cc \
	shortest_cusp_basis.cc \
	simplify_triangulation.cc \
	sl2c_matrices.cc sl2c_matrices.h \
	solve_equations.cc \
	subdivide.cc \
	symmetric_group.cc \
	symmetry_group.cc symmetry_group.h \
	symmetry_group_closed.cc \
	symmetry_group_cusped.cc \
	symmetry_group_info.cc \
	tables.cc tables.h \
	terse_triangulation.cc terse_triangulation.h \
	tersest_triangulation.cc tersest_triangulation.h \
	tet_shapes.cc \
	tidy_peripheral_curves.cc \
	transcendentals.cc \
	triangulation.h \
	triangulation_io.h \
	triangulations.cc \
	two_bridge.cc \
	unix_io.cc unix_io.h \
	update_shapes.cc \
	volume.cc \
	winged_edge.h

d2d = Alg_matrices.cc Alg_matrices.h \
	Dirichlet2D.cc Dirichlet2D.h \
	Dirichlet2D_construction.cc \
	Dirichlet2D_extras.cc \
	Dirichlet2D_prototypes.h \
	Dirichlet_basepoint2D.cc \
	Dirichlet_precision2D.cc \
	SL2RMatrix.cc SL2RMatrix.h \
	polygon.h

dt = dt2snap.cc \
	crossing_signs.cc \
	link_diagram.cc link_diagram.h 

INCLUDES = -I.. -I@PARI_DIR@/include -I$(srcdir)/../@PWR_DIR@
AM_LDFLAGS = -L@PARI_DIR@/lib -L$(top_builddir)/@PWR_DIR@

libsnap_a_SOURCES = $(snappea) $(d2d) $(dt)
libsnap_a_CXXFLAGS = -DNOPARI

libpsnap_a_SOURCES = $(libsnap_a_SOURCES)
libpsnap_a_CXXFLAGS = 
