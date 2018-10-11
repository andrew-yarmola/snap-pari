SUBDIRS = src

nobase_dist_data_DATA = snap_data/closed.fields \
	snap_data/cusped.fields \
	snap_data/closed_census_algebras \
	snap_data/cusped_census_fields \
	snap_data/eta_fudges \
	snap_data/sample_session \
	snap_data/snap_help snap_data/tube_help \
	snap_data/trs5 snap_data/trs6 snap_data/trs7 \
	snap_data/trn.m.txt \
	snap_data/ClosedManifolds

EXTRA_DIST = paper snap_data/link_data snap_data/manifolds

install-data-local:
	for d in link_data manifolds; do \
	  $(mkinstalldirs) $(DESTDIR)$(datadir)/snap_data/$$d ;\
	  ( cd $(top_srcdir) ; $(AMTAR) cf - snap_data/$$d ) | \
	  ( cd $(DESTDIR)$(datadir) ; $(AMTAR) xf - ) ;\
	  chmod -R u+w $(DESTDIR)$(datadir)/snap_data/$$d ;\
	done

uninstall-local:
	for d in link_data manifolds; do \
	  chmod -R a+w $(DESTDIR)$(datadir)/snap_data/$$d ;\
	  rm -rf $(DESTDIR)$(datadir)/snap_data/$$d ;\
	done
