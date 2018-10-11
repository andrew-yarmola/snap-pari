#ifndef _unix_io_
#define _unix_io_

#include "SnapPea.h"
#include <stdio.h>
#include <string>

Triangulation* read_manifold_file(FILE* fp);
Triangulation* read_census_manifold(const char *path, int which_census, int which_manifold);

// FILE* locate_file(const std::string& path, const std::string& name, char *mode);
// FuncResult read_manifold_file(const char *path, const char *file, Triangulation *manifold);
// FuncResult read_manifold_file(FILE* fp, Triangulation *manifold);
// FuncResult read_census_manifold(const char *path, int which_census, int which_manifold, Triangulation *manifold);
// Triangulation* read_manifold_file(const char *path, const char *file);


void write_manifold_file(FILE *fp, Triangulation *manifold);
void write_old_manifold_file(FILE *fp, Triangulation *manifold);
int read_matrix_generators(FILE *fp, O31_matrix **mats); 

bool find_census_manifold(const char* path, Triangulation* T, int& census, int& n);

#endif
