#ifndef _terse_io_
#define _terse_io_

#include "snappea/SnapPea.h"
#include <string>

bool save_terse(std::string const& dir, Triangulation* T, int& index);
Triangulation* read_terse(std::string const& dir, int n_tet, int index);
int num_terse(std::string const& dir, int n_tet);

class manifold_database {
  class manifold_db* rep;
public:
  manifold_database(std::string const& dir);
  ~manifold_database();

  bool locate(Triangulation* T, int& index, int& nt) const;
  bool insert(Triangulation* T, int& index, int& nt);

  void save_pending_keys();
  void print() const; 
};

#endif
