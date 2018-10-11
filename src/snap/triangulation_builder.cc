/*
** Copyright (C) 2003 Oliver A. Goodman <oag@ms.unimelb.edu.au>
**  
** This file is part of Snap.
** 
** Snap is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software 
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "triangulation_builder.hh"

using std::cout;
using std::endl;

/* The purpose of a triangulation_builder is to construct the ideal
   triangulation dual to a given tube domain (or Ford domain). A tube
   domain is a polyhedral torus with face pairings. To describe a tube
   domain we first number the faces and number the vertices within
   each face. Each corner of each face is then given by a pair of
   numbers. (We`ll ignore ideal vertices for the moment.) We also
   orient all the faces consistently, so that each corner has an
   incoming and an outgoing edge.  For each corner we give the corner
   it is paired with under the face pairings, the corner adjacent to
   its incoming edge, and the next corner around the same
   face. Currently we assume that the tube domain is for an orientable
   manifold and therefore that every face pairing reverses face
   orientation.

   We can also give several tori at the same time. The method
   of description is unchanged, and face pairings may be between
   faces of one torus or between different tori. In this situation 
   the dual ideal triangulation will have several cusps. 

   The input describes combinatorially a polyhedral decomposition of
   the manifold. The dual will be a decomposition of a cusped manifold
   into ideal cells. These cells can then be further subdivided into
   tetrahedra. If every vertex on the tube domain torus has order 3
   and, after face pairings, every edge has 3 faces around it, the
   cells will all be tetrahedra. This is the generic situation. 

   We shall attempt to cope with higher order vertices, leaving higher
   order edges aside for the moment. This means our cells are
   'deltahedra,' ie. all the faces are triangles. To turn such a cell
   into a collection of tetrahedra, we choose one vertex, and cone all
   non-adjacent (triangular) faces down to this vertex. Between
   tetrahedra there are now two kinds of face pairing which we shall
   call "internal": between tetrahedra in the same cell, and
   "external": arising from face pairings between cells. 

   In order to get the right picture in the following, you should
   assume that, looking at the tube from the outside, the faces have
   been oriented clockwise. The dual cells will be constructed with
   the same orientation, clockwise from the outside looking in.

   Our first task, given the input, is to compute the collection of
   dual cells. We group together face_corners which share a common
   vertex, ie. corners which are either paired or adjacent. We can
   identify each corner with a (face,vertex) pair of the cell as
   follows: the incoming edge of a corner determines a face. The
   right-hand rule gives an orientation for the face (cw as we look at
   the face towards the corner along the incoming edge). The corner
   itself cuts through one edge of the face. The vertex in question is
   the first one around from this edge in the direction of the
   orientation. A whole cell is then a set of cell_faces where each
   cell_face is a set of (face,vertex) pairs. (If we are dealing only
   with deltahedra, each cell_face is a triple of (face,vertex) pairs. 

   Now, each (face,vertex) is represented by one of our original tube
   face_corners. The relations of pairing, adjacency and next-ness can
   readily be obtained from those of the face_corners. They are as
   follows:

   (face,vertex)        face_corner
   -------------        -----------
   paired               adjacent.next
   adjacent             rev_adjacent
   next                 adjacent.paired

   Let`s stop here and see if we can write code to compute and print
   out a full cell structure.

   To turn cells into tetrahedra, we first simply count how many
   tetrahedra will be needed. For each cell it will be the total
   number of faces, minus the number of faces around it`s base
   vertex. We`ll use cell[0][0] as the base vertex. 

   

*/

void triangulation_builder::print_faces() const
{
  fc_map::const_iterator it; 
  for (it = paired.begin(); it != paired.end(); it++) {
    cout << (*it).first << ' ' << (*it).second << ' ' << ((fc_map&)adjacent)[(*it).first] << ' ';
    cout << ((fc_map&)next)[(*it).first] << endl; 
  }
}

inline bool contains(face_corner const& fc, fc_map const& m)
{
  return m.find(fc)!=m.end();
}

face_corner triangulation_builder::dual_edge(face_vertex const& fv) const
{ 
  fc_map::const_iterator it = tfc_dual_edge.find(fv);
  if (it != tfc_dual_edge.end()) {
    face_corner fc = (*it).second;
    fc.face = fc.face % 256;
    return fc;
  }
  return face_corner(); 
}

void triangulation_builder::print_cell_complex() const
{
  int ncells = CX.size(), nvert, nfaces; 
  cout << "Complex has " << ncells << " cells.\n";
  int i, j, k; 

  face_vertex fv, fv0, ind, pfv; 
  bool found; 
  fc_map done; 
  fc_map vx_class; 
  int nf, nvx; 

  for (i=0; i<ncells; i++) {

    // Find vertex classes for this cell. 
    nvx = 0; 
    while (true) {
      // Find corner not yet done. 
      nfaces = CX[i].size();
      found = false; 
      for (j=0; !found && j<nfaces; j++) {
	nvert = CX[i][j].size(); 
	for (k=0; !found && k<nvert; k++) {
	  fv = CX[i][j][k]; 
	  if (!contains(fv,done)) found=true; 
	}
      }
      if (!found) break; 

      // Find adjacency cycle starting at fv.
      fv0 = fv; 
      nf = 0; 
      do {
	ind = ((fc_map&)index)[fv]; 
	ind.face = ind.face % 256; 
	vx_class[fv] = face_vertex(ind.face, nvx); 
	done[fv] = fv; 
	fv = fv_adjacent(fv); 
      } while (fv!=fv0 && ++nf < 10); 

      ++nvx; 
    }
  }

  for (i=0; i<ncells; i++) {
    cout << "Cell: " << i << endl; 

    // Print the vertices on each face
    nfaces = CX[i].size();
    for (j=0; j<nfaces; j++) {
      nvert = CX[i][j].size(); 
      cout << j << "["; 
      for (k=0; k<nvert; k++) {
	fv = CX[i][j][k]; 
	if (k) cout << ',';
	cout << vx_class[fv].vertex; 
      }

      pfv = fv_paired(fv);
      cout << "] -> " << ((fc_map&)index)[pfv].face/256 << ':' << 
	((fc_map&)vx_class)[pfv].face << "["; 

      for (k=0; k<nvert; k++) {
	pfv = fv_paired(CX[i][j][k]); 
	if (k) cout << ',';
	cout << vx_class[pfv].vertex; 
      }
      cout << "]" << endl; 
    }
  }

}

static bool find_fv_in_set_diff(fc_map const& a, fc_map const& b, face_vertex& fv)
{
  fc_map::const_iterator it; 
  for (it = a.begin(); it!=a.end(); it++) {
    fv = (*it).first; if (!contains(fv, b)) break; 
  }
  return it!=a.end(); 
}

bool triangulation_builder::compute_cell_complex()
{
  if (!ok) return false; 

  CX.clear();
  CX.reserve(10); // Make a little bit of space to be going on with.
  index.clear(); 

  fc_map todo, done; 
  face_corner fc, fc0, fc2; 
  fc_map::const_iterator it; 
  int cnum=-1, fnum, vnum;  
  vector<cell>::iterator the_cell; 
  cell::iterator the_face; 

  while (done.size() < paired.size() && ++cnum < 100) {

    // Start a new cell. 
    CX.push_back(cell()); 
    (the_cell = CX.end())--; 
    the_cell->reserve(10); 

    // Find a face_vertex to start with. 
    if (!find_fv_in_set_diff(paired,done,fc)) { 
      cout << "Ran out of face_vertices prematurely.\n"; 
      return false; 
    }
    todo[fc] = fc; 

    // Add the faces. 
    fnum = -1; 
    while(todo.size() && ++fnum < 200) {

      fc = fc0 = (*todo.begin()).first; 

      // Start a new face. 
      the_cell->push_back(cell_face()); 
      (the_face = the_cell->end())--; 
      the_face->reserve(3); 

      // Add vertices to the_face. 
      vnum = 0; 
      do {
	the_face->push_back(fc); 

	index[fc] = face_vertex(256*cnum + fnum, vnum); 

	fc2 = fv_adjacent(fc); 
	if (!contains(fc2, done)) {
	  todo[fc2] = fc2; 
	}

	done[fc] = fc; 
	todo.erase(fc); 

	fc = fv_next(fc); 
      } while (fc!=fc0 && ++vnum < 15);

      if (vnum==15) {
	cout << "There were too many vertices on a face (more than 15).\n"; 
	return false; 
      }

    }

    if (fnum==200) { 
      cout << "There were too many faces (more than 200).\n"; 
      return false; 
    }
  }

  if (cnum==100) {
    cout << "There were too many cells (more than 100).\n"; 
    return false; 
  }

  return true; 
}


const char* type_letter = "BFTI"; 

bool operator == (face_corner const& a, face_corner const& b)
{
  return a.face==b.face && a.vertex==b.vertex && a.type==b.type; 
}

bool operator != (face_corner const& a, face_corner const& b)
{
  return a.face!=b.face || a.vertex!=b.vertex || a.type!=b.type;
}

bool operator < (face_corner const& a, face_corner const& b)
{
  if (a.face < b.face) return true; 
  if (a.face > b.face) return false; 
  if (a.vertex < b.vertex) return true; // now a.face==b.face
  if (a.vertex > b.vertex) return false; 
  return a.type < b.type; // now a.vertex==b.vertex

  /* 
     return a.face <  b.face || 
     (a.face == b.face && ( a.vertex <  b.vertex || 
     (a.vertex == b.vertex && a.type < b.type)));
   */
}

ostream& operator << (ostream& out, face_corner const& fc)
{
  if (fc.type==fc_back)
    return out << '(' << fc.face << ',' << fc.vertex << ')';
  return out << '(' << fc.face << ',' << fc.vertex << ',' << type_letter[fc.type] << ')';
}

void triangulation_builder::glue_face(int n)
{
  int m = paired[face_corner(n, 0)].face;
  face_corner adj;
  triangulation_builder T;
  fc_map::const_iterator it; 
  for (it = paired.begin(); it != paired.end(); it++) {
    if ((*it).first.face==m || (*it).first.face==n) 
      continue; 
    adj = adjacent[(*it).first];
    if (adj.face==m || adj.face==n)
      adj = adjacent[paired[adj]];

    T.add_face_corner((*it).first, (*it).second, adj, next[(*it).first]); 
  }
  paired = T.paired; 
  adjacent = T.adjacent;
  rev_adjacent = T.rev_adjacent;
  next = T.next; 
}

void triangulation_builder::add_face_corner(face_corner const& f, face_corner const& p, face_corner const& a, face_corner const& n)
{
  map<face_corner,face_corner>::const_iterator i; 

  if ((i = paired.find(f))==paired.end()) paired[f] = p;
  else { if ((*i).second != p) error("bad pairing"); }

  if ((i = adjacent.find(f))==adjacent.end()) { adjacent[f] = a; rev_adjacent[a] = f; }
  else { if ((*i).second != a) error("bad adjacency"); }

  if ((i = next.find(f))==next.end()) next[f] = n;
  else { if ((*i).second != n) error("bad next"); }
}

void triangulation_builder::error(const char* msg)
{
  if (!ok) return; 

  // err_msg = msg; 
  ok = false; 
  cout << "triangulation_builder: " << msg << ".\n"; 
}

bool triangulation_builder::all_cell_faces_triangular() const
{
  int i, j, nfaces, ncells = CX.size(); 
  for (i=0; i<ncells; i++) {
    nfaces = CX[i].size(); 
    for (j=0; j<nfaces; j++) {
      if (CX[i][j].size()!=3) return false; 
    }
  }
  return true; 
}

bool triangulation_builder::get_triangulation_data(TriangulationData& TD, int report) 
{
  if (!ok) return false; 

  int ncells = CX.size(), nfaces; 
  int i, j, k; 

  // Check all faces are triangular or if not, make it so. 

  if (!all_cell_faces_triangular()) {
    split_all_cell_faces(); 
    compute_cell_complex(); 
    if (!all_cell_faces_triangular()) {
      cout << "Problem splitting all non-triangular cell faces.\n";
      return false; 
    }
  }

  // Count tetrahedra needed.

  face_vertex fv, fv0; 
  int nf; 
  fc_map base_faces;
  int ntet = 0; 
  face_vertex f;
  for (i=0; i<ncells; i++) {

    fv0 = fv = CX[i][0][0]; 
    nf = 0; 
    do {
      fv = fv_adjacent(fv); ++nf; 
      f = face_vertex(index[fv].face,0); 
      base_faces[f] = f; 
    } while (fv!=fv0 && nf < 10); 

    if (fv!=fv0) {
      cout << "Too many faces around one cell vertex (more than 10).\n"; 
      return false;
    }

    ntet += CX[i].size() - nf; 
  }

  TetrahedronData* tet = NEW_ARRAY(ntet, TetrahedronData);

  // Set up the TriangulationData. 

  TD.num_tetrahedra = ntet; 
  TD.tetrahedron_data = tet; 

  TD.num_or_cusps = 0; 
  TD.num_nonor_cusps = 0; 
  TD.cusp_data = 0;

  TD.orientability = oriented_manifold;
  TD.CS_value_is_known = FALSE; 

  // Complete gluing fields for top triangles. 

  int tn = 0; 
  fc_map glued_tet; 
  for (i=0; i<ncells; i++) {
    nfaces = CX[i].size(); 
    for (j=0; j<nfaces; j++) {

      // Check if this is a top face. 
      if (contains(face_vertex(256*i+j,0),base_faces)) continue; 

      // Complete gluing info for a top face. 
      // This says that face 0 of tet tn is glued to face j of cell i. 

      tet[tn].neighbor_index[0] = i + ntet; 
      tet[tn].gluing[0][0] = j; 
      for (k=0; k<3; k++) {
	tet[tn].gluing[0][k+1] = k;
	glued_tet[CX[i][j][k]] = face_vertex(4*tn, k+1); 
	tfc_dual_edge[face_vertex(4*tn, k+1)] = CX[i][j][k];
      }

      if (report) cout << "Tet:" << tn << " face:0 <-> Cell:" << i << " face:" << j << endl;

      ++tn;
    }
  }

  // Now do gluing info for bottom triangles and internal gluings. 

  face_corner s1, s2, s3, c1, c2, c3;
  face_vertex tc1, tc2, tc, ts; 
  int tf, cf; 

  for (tn=0; tn<ntet; ++tn) {
    
    // Do top face glued to tet tn. 
    i = tet[tn].neighbor_index[0] - ntet;
    j = tet[tn].gluing[0][0]; 

    // Do each edge on this top face.
    for (k=0; k<3; k++) {

      c1 = CX[i][j][k];
      c2 = fv_next(c1); 
      c3 = fv_next(c2); 
      s1 = fv_adjacent(c2); 
      s2 = fv_next(s1); 
      s3 = fv_next(s2); 

      // Have c1,c2,c3,s1,s2,s3 around an edge.
      // c1,c2 adjacent to s2,s1.

      tc1 = glued_tet[c1];
      tc2 = glued_tet[c2]; 
      tc = glued_tet[c3];

      tf = tc.vertex; // Face opposite vertex at c3. 
      cf = index[s3].face; 

      if (contains(face_vertex(cf,0),base_faces)) {
	// Side edge: do bottom triangle gluing.

	glued_tet[s1] = face_vertex(4*tn+tf, tc2.vertex);
	glued_tet[s2] = face_vertex(4*tn+tf, tc1.vertex);
	glued_tet[s3] = face_vertex(4*tn+tf, 0); 

	tfc_dual_edge[face_vertex(4*tn+tf, tc2.vertex)] = s1;
	tfc_dual_edge[face_vertex(4*tn+tf, tc1.vertex)] = s2;
	tfc_dual_edge[face_vertex(4*tn+tf, 0)] = s3;

	tet[tn].neighbor_index[tf] = i + ntet; 
	tet[tn].gluing[tf][tf] = cf % 256; 
	tet[tn].gluing[tf][0] = index[s3].vertex;
	tet[tn].gluing[tf][tc2.vertex] = index[s1].vertex;
	tet[tn].gluing[tf][tc1.vertex] = index[s2].vertex;

	if (report) 
	  cout << "Tet:" << tn << " face:" << tf << 
	    " <-> Cell:" << i << " face:" << (cf % 256) << endl;
	
      } else {
	// Top edge: do internal gluing.

	ts = glued_tet[s3];

	tet[tn].neighbor_index[tf] = ts.face/4; 
	tet[tn].gluing[tf][tf] = ts.vertex; 
	tet[tn].gluing[tf][0] = 0;
	tet[tn].gluing[tf][tc1.vertex] = glued_tet[s2].vertex;
	tet[tn].gluing[tf][tc2.vertex] = glued_tet[s1].vertex;

	if (report)
	  cout << "Tet:" << tn << " face:" << tf << 
	    " <-> Tet:" << (ts.face/4) << " face:" << (ts.vertex) << endl; 

      }
    }
  }

  // Glue tetrahedra to each other rather than to cell faces. 

  face_vertex pfv; 
  int tv;
  for (tn=0; tn<ntet; ++tn) {
    for (tf=0; tf<4; ++tf) {

      i = tet[tn].neighbor_index[tf] - ntet;
      if (i < 0) continue; 
      j = tet[tn].gluing[tf][tf]; 

      for (tv=0; tv<4; ++tv) {
	if (tv==tf) continue; 
	k = tet[tn].gluing[tf][tv]; 
	fv = CX[i][j][k]; 
	pfv = fv_paired(fv); 
	ts = glued_tet[pfv];
	tet[tn].gluing[tf][tv] = ts.vertex; 
      }
      tet[tn].gluing[tf][tf] = ts.face%4;
      tet[tn].neighbor_index[tf] = ts.face/4; 
    }
  }

  return true; 
}


face_vertex triangulation_builder::fv_next(face_vertex const& fv) const
{
  return ((fc_map&)paired)[((fc_map&)adjacent)[fv]];
}

face_vertex triangulation_builder::fv_adjacent(face_vertex const& fv) const
{
  return ((fc_map&)rev_adjacent)[fv]; 
}

face_vertex triangulation_builder::fv_paired(face_vertex const& fv) const
{
  return ((fc_map&)next)[((fc_map&)adjacent)[fv]];
}

/* The following function operates on the combinatorial tube rather
   than on the cell complex. After calling it, it is necessary to
   recreate the complex. In spite of this its operation is easier to
   understand in terms of what it does to the complex. Its basic
   function is to take a face of the complex having more than 3 edges
   and split off a triangle (doing the same on the paired face
   simultaneously). In terms of the combinatorial tube what it does is
   take an edge where more than 3 faces meet and split it into two
   edges separated by a bigon.

   We start at corner X on the cell face, take two steps around it to
   find corner Y. We shall split the face from X to Y. On the paired
   face, the partners of X and Y are Z and W respectively. Splitting
   creates four new corners, A, B, C and D (reverse) adjacent to X, Y,
   Z and W respectively.

   In terms of the combinatorial tube, A and C form one bigon, pairing
   with B and D on its partner. A will be inserted between X and its
   adjacent corner X1.

   Since X and A really represent the same edge of the original tube
   they are numbered in such a way that we can recover X from A. In
   fact what we do is make A be X + 256*offset. Actually, in case X is
   already the result of an earlier splitting we make A be Xmod256 +
   256*offset. Then the original edge associated with a split tube can
   be recovered by taking the face_corner mod 256.  
*/

void triangulation_builder::split_cell_face(face_vertex const& X, int offset)
{
  face_vertex Y, Z, W, A, B, C, D;

  Y = fv_next(fv_next(X)); 
  Z = fv_paired(X);
  W = fv_paired(Y); 

  A = X; A.face = (X.face % 256) + 256*offset;
  B = Y; B.face = (Y.face % 256) + 256*offset;
  C = Z; C.face = (Z.face % 256) + 256*offset;
  D = W; D.face = (W.face % 256) + 256*offset;

  face_vertex X1, Y1, Z1, W1;

  X1 = adjacent[X];
  Y1 = adjacent[Y];
  Z1 = adjacent[Z];
  W1 = adjacent[W];

  // Now we've got all the new and existing corners we need, update
  // the fc_maps, so that A,B,C and D are inserted appropriately. 

  // Make 2 bigons.
  next[A] = C; next[C] = A; 
  next[B] = D; next[D] = B; 

  // Pair them. 
  paired[A] = B; paired[B] = A;
  paired[C] = D; paired[D] = C;

  // Insert them between existing faces. 
  adjacent[X] = A; adjacent[A] = X1; rev_adjacent[X1] = A; rev_adjacent[A] = X; 
  adjacent[Y] = B; adjacent[B] = Y1; rev_adjacent[Y1] = B; rev_adjacent[B] = Y; 
  adjacent[Z] = C; adjacent[C] = Z1; rev_adjacent[Z1] = C; rev_adjacent[C] = Z; 
  adjacent[W] = D; adjacent[D] = W1; rev_adjacent[W1] = D; rev_adjacent[D] = W; 
}

void triangulation_builder::split_all_cell_faces()
{
  fc_map done; 
  face_vertex fv, fv0; 
  int order, splitnum=1;  

  while (splitnum < 30) {

    // Find corner we haven't done yet. 
    if (!find_fv_in_set_diff(paired, done, fv)) break; 

    // Find the order of this edge. 
    fv0 = fv;
    order = 0;
    do {
      fv = fv_next(fv);
      order++;
      done[fv] = fv;
    } while (fv != fv0 && order < 10);
    if (fv != fv0) {
      cout << "Too many tube faces around one edge (more than 10).\n"; 
      return; 
    }

    if (order < 4) continue; 

    cout << "Splitting a cell face.\n";
    split_cell_face(fv, splitnum);
    
    splitnum++;
    done.clear();
  }
}

void triangulation_builder::eliminate(face_corner const& e)
{
  paired.erase(e);
  adjacent.erase(e);
  next.erase(e);
  rev_adjacent.erase(e); 
}

void triangulation_builder::eliminate_edge(face_corner const& e)
{
  face_corner ea = adjacent[e];
  face_corner er = rev_adjacent[e];
  next[ea] = next[e];
  paired[ea] = paired[e];
  adjacent[er] = ea;
  rev_adjacent[ea] = er; 
  eliminate(e); 
}

int triangulation_builder::vertex_order(face_corner const& fc)
{
  int order = 0;
  face_corner fc2 = fc; 
  do {
    fc2 = adjacent[fc2];
    order++;
  } while (fc2!=fc && order < 1000);
  return order; 
}

void triangulation_builder::remove_redundant_edges(int report)
{
  fc_map todo = paired; 
  face_corner fc, fc2;
  face_corner e[4];
  int ord, i, count=0; 

  while (todo.size()) {

    fc = (*todo.begin()).first; 

    // Find the order of this edge & erase it from todo list. 
    fc2 = fc; 
    ord = 0;
    do {
      todo.erase(fc2); 
      todo.erase(next[adjacent[fc2]]);
      ord++; 
      fc2 = paired[adjacent[fc2]]; 
    } while (fc2!=fc);

    if (ord==1) {
      cout << "Warning, edge of order 1 encountered!\n";
      continue; 
    }

    if (ord > 2) continue;

    // Get the 4 face edges which are identified with this order 2 edge. 
    e[0] = fc; 
    e[2] = paired[adjacent[fc]];
    e[1] = next[adjacent[fc]];
    e[3] = paired[adjacent[e[1]]];
    
    if (report) {
      cout << "Removing edge: "; 
      for (i=0; i<4; i++) cout << e[i] << ' ';
      cout << endl;
    }

    // Check for bad or redundant vertices at ends of our edge. 
#if 0
    for (i=0; i<4; i++) {
      if (vertex_order(e[i]) < 3) break; 
    }
    if (i < 4) {
      cout << "Bad or redundant vertex on redundant edge.. giving up.\n";
      return; 
    }
#endif

    // Eliminate the redundant edge.
    for (i=0; i<4; i++) eliminate_edge(e[i]); 

    count++;
  }
  if (count > 0) {
    if (count==1) 
      cout << "Removed 1 redundant edge.\n";
    else 
      cout << "Removed " << count << " redundant edges.\n";
  }
}


void triangulation_builder::remove_redundant_vertices()
{
  fc_map todo = paired; 
  face_corner fc, fc2, next_fc, rev_next, ad;
  int ord, count=0, vo; 
  bool redundant, non_redundant; 

  while (todo.size()) {

    fc = (*todo.begin()).first; 

    // Erase edge from todo list. 
    fc2 = fc; 
    ord = 0;
    redundant = false; 
    non_redundant = false; 
    do {
      vo = vertex_order(fc2); 
      if (vo==1) {
	cout << "Vertex of order 1 encountered, giving up!\n";
	return;
      } 
      if (vo==2) redundant = true;
      else non_redundant = true; 

      todo.erase(fc2); 
      ord++; 
      fc2 = paired[adjacent[fc2]]; 
    } while (fc2!=fc && ord < 1000);

    // Check we know what we've got. 
    if (ord==1000) {
      cout << "Edge order 1000+, giving up.\n";
      return; 
    }
    if (redundant && non_redundant) {
      cout << "Bad vertex encountered, giving up.\n";
      return;
    }

    if (non_redundant) continue; 

    // Step around the edge again, bypassing the redundant vertices. 
    do {
      ad = adjacent[fc2];
      todo.erase(ad); 

      // Bypass fc2 and ad. 
      rev_next = paired[next[paired[fc2]]];
      next[rev_next] = next[fc2]; 
      rev_next = paired[next[paired[ad]]];
      next[rev_next] = next[ad]; 

      fc2 = paired[ad]; 
    } while (fc2!=fc);

    // Step around the edge again, removing the redundant vertices. 
    do {
      ad = adjacent[fc2];
      next_fc = paired[ad]; 

      // Remove fc2 and ad entirely.
      eliminate(fc2); 
      eliminate(ad);
 
      fc2 = next_fc;
    } while (fc2!=fc);

    count++;
  }
  if (count > 0) {
    if (count==1) 
      cout << "Removed 1 redundant vertex.\n";
    else 
      cout << "Removed " << count << " redundant vertices.\n";
  }
}
