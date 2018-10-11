#include "dirichlet.hh"

#include <set>
#include <algorithm>
#include <iomanip>

#include "p_VS.hh"
#include "polytope_T.hh"
#include "O31_pari.hh"
#include "to_pari.hh"
#include "printable.hh"

#include "snappea/fundamental_group.h"
#include "snappea/winged_edge.h"
#include "kernel_extras.hh"

using std::min;
using std::max;
using std::set;
using std::string; 
using std::setw; 

struct group_element {
  FGWord w;
  MoebiusTransformation m;
  pari pt[2]; // pt[0] = m*origin;

  const face<p_VS> *f[2];

  void set_pt(); 

  group_element() { f[0]=f[1]=0; }
  group_element(FGWord const& w0, MoebiusTransformation const& m0)
    : w(w0), m(m0) { set_pt(); f[0]=f[1]=0; }
  group_element(int) : m(One)
  { m.acc_matrix = p_lisexpr("[1,0;0,1]"); set_pt(); f[0]=f[1]=0; }

  double height() const { return pt[0][0].double_value(); }
  void kill_face(int i) const { ((group_element*)this)->f[i]=0; }
  double accuracy(GroupPresentation* G) const; 
};

struct D_domain {
  polytope<p_VS> D;
  set<group_element> Live;

  D_domain();
  bool compute(GroupPresentation* G); 
  bool group_contains(group_element& e) const; 
  bool angle_sums_2pi(bool report) const;
};

static pari hdist(pari const& a, pari const& b)
{
  pari mindot = -O1n_innerprod(a,b);
  pari d = (mindot < pONE) ? pZERO : acosh(mindot); 
  return d; 
}

static pari dhangle(pari const& a, pari const& b)
{
  pari mindot = -O1n_innerprod(a,b)/sqrt(O1n_normsq(a)*O1n_normsq(b));
  return acos(mindot); 
}

static double match_error=0.;

static bool small(double x)
{
  double s = fabs(x);
  bool match = s < p_VS::eq_EPS; 
  if (s < p_VS::ne_EPS) {
    if (s > match_error) match_error = s; 
    if (!match) cout << "indeterminate size comparison\n";
  }
  return match; 
}

ostream& operator << (ostream& out, group_element const& e)
{ 
  return out << e.pt[0] << ' ' << e.w;
}

group_element inverse(group_element const& e)
{
  group_element ei;
  ei.w = inverse(e.w);
  Moebius_invert(&e.m, &ei.m);
  ei.pt[0] = e.pt[1];
  ei.pt[1] = e.pt[0];
  return ei; 
}

void group_element::set_pt()
{
  pari o31m = SL2C_to_O31(to_pari(m));
  pt[0] =  col(o31m,0);
  pt[1] = -row(o31m,0);
  pt[1][0] = -pt[1][0];
}

group_element operator * (group_element const& a, group_element const& b)
{
  group_element p;
  p.w = a.w * b.w;
  p.m = a.m * b.m; 
  p.set_pt(); 
  return p; 
}

static bool near(pari const& a, pari const& b)
{
  return small(sqrt(size(gnorml2(b/b[0] - a/a[0]))));
}


// Code for face_pair which enables us to iterate 
// over the edges forming a single edge of the glued domain. 

typedef list<edge<p_VS> > edge_list; 
typedef list<face<p_VS> > face_list;

static bool incident(edge<p_VS> const& e, face<p_VS> const& f)
{
  return (find(f.I.begin(), f.I.end(), e.a)!=f.I.end() && 
	  find(f.I.begin(), f.I.end(), e.b)!=f.I.end());
}

#if 0
static edge_list::iterator find_an_edge(edge_list& el, face<p_VS> const& f) 
{
  edge_list::iterator it;
  for (it=el.begin(); it!=el.end(); ++it) 
    if (incident(*it, f)) break; 
  return it; 
}
#endif

static bool find_edge(edge_list const& el, pari const& a, pari const& b, edge_list::const_iterator& it)
{
  for (it=el.begin(); it!=el.end(); ++it) 
    if ((near(it->a->V, a) && near(it->b->V, b)) ||
	(near(it->a->V, b) && near(it->b->V, a))) break; 
  return it!=el.end(); 
}

static void find_faces(face_list const& fl, edge<p_VS> const& e, 
		       face_list::const_iterator& fa, 
		       face_list::const_iterator& fb)
{
  for (fa=fl.begin(); fa!=fl.end(); ++fa)
    if (incident(e, *fa)) break; 
  if (fa==fl.end()) { fb=fl.end(); return; }
  for (fb=fa,++fb; fb!=fl.end(); ++fb)
    if (incident(e, *fb)) break;
  return; 
}

class face_pair {
  const D_domain *D;
  face_list::const_iterator fa, fb;
  edge_list::const_iterator e;
public:
  face_pair() {}
  face_pair(const D_domain *D0, edge_list::const_iterator const& it) 
    : D(D0), e(it)
  { find_faces(D->D.FL, *e, fa, fb); }

  edge_list::const_iterator edge() const { return e; }
  pari angle() const { return dhangle(fa->F, fb->F); }

  face_pair& operator ++();
  friend bool operator == (face_pair const& a, face_pair const& b);
};

bool operator == (face_pair const& a, face_pair const& b)
{
  return a.e==b.e && a.fa == b.fa;
}

face_pair& face_pair::operator ++()
{
  // these are the things we want to find. 
  face_list::const_iterator nfa, nfb;
  edge_list::const_iterator ne;

  // find partner for fb
  set<group_element>::const_iterator ei; 
  int k; 
  const face<p_VS> *fp = &*fb;
  for (ei=D->Live.begin(); ei!=D->Live.end(); ++ei) {
    for (k=0; k<2; ++k) 
      if (ei->f[k]==fp) break;
    if (k<2) break; 
  }
  if (k==2) {
    cout << "unable to find pairing for a face\n";
    return *this; 
  }
  for (nfa=D->D.FL.begin(); nfa!=D->D.FL.end(); ++nfa)
    if (&*nfa==ei->f[1-k]) break; 
  if (nfa==D->D.FL.end()) {
    cout << "unable to find face of a D_domain\n";
    return *this;
  }

  // find partner edge on nfa
  pari pairing = SL2C_to_O31(to_pari((k==1) ? ei->m : inverse(ei->m)));
  if (!find_edge(D->D.EL, pairing * e->a->V, pairing * e->b->V, ne)) {
    cout << "unable to find partner for a D_domain edge\n"; 
    return *this;
  }

  // find nfb
  face_list::const_iterator f1,f2; 
  find_faces(D->D.FL, *ne, f1, f2); 
  if (nfa==f1) nfb = f2; 
  else if (nfa==f2) nfb = f1; 
  else {
    cout << "paired edge and face do not appear to be incident\n"; 
    return *this; 
  }

  fa = nfa;
  fb = nfb;
  e  = ne; 
  return *this; 
}

// End of face_pair code. 


bool D_domain::angle_sums_2pi(bool report) const
{
  // Make a list of the edges. 
  edge_list::const_iterator e; 
  typedef list<edge_list::const_iterator> leli; 
  leli le, eclass;
  for (e=D.EL.begin(); e!=D.EL.end(); ++e) 
    le.push_back(e);
  
  leli::iterator lei; 
  face_pair fp, fp0; 
  pari asum; 
  while (le.size()) {
    asum = pZERO;
    e = le.front(); 
    eclass.push_back(e); 
    // le.pop_front(); 
    fp0 = fp = face_pair(this, e); 
    do {
      ++fp; 
      lei = find(le.begin(), le.end(), fp.edge());
      if (lei!=le.end()) 
	le.erase(lei);
      else { cout << "edge appears twice in cycle\n"; return false; }
      // if (report) cout << *fp.edge() << ' ';
      asum += fp.angle(); 
    } while (!(fp==fp0));
    if (report) asum.print(); 
  }
  return true; 
}


bool D_domain::group_contains(group_element& e) const
{
  pari origin = cvector(4); origin[0]=pONE; 

  if (!Live.size()) 
    return small(hdist(e.pt[0],origin).double_value()); 
  set<group_element>::const_iterator i, ci; 
  int k, ck, count=25; 
  pari min_dist, d; 
  while (true) {

    // find nearest to e^-1.
    ck = -1; 
    min_dist=hdist(e.pt[1], origin);
    if (small(min_dist.double_value())) break; 

    for (i=Live.begin(); i!=Live.end(); ++i) {
      for (k=0; k<2; ++k) {
	d = hdist(e.pt[0], i->pt[k]);
	if (small((min_dist - d).double_value())) continue; // equal closest
	if (d > min_dist) continue; 
	min_dist = d; 
	ci = i; 
	ck = k; 
      }
    }
    if (ck == -1) break; // origin is closest
    if (--count==0) break; 

    e = (ck==1) ? *ci * e : inverse(*ci) * e;
  }
  if (!count) 
    cout << "D_domain::group_contains failed to converge\n";
  return small(min_dist.double_value());
}

#if 0
static inline int real_prec(pari const& r)
{
  return int(double(length(abs(r))-2)* pariK);
}

inline static bool small(pari const& x)
{
  return small(size(x)); 
}
#endif

double group_element::accuracy(GroupPresentation* G) const
{
#if 0
  pari om = SL2C_to_O31(to_pari(m));
  pari mt = O31_to_SL2C(om)[0];
#endif
  pari mt = word_to_Moebius(G, w).acc_matrix; 

  if (gabs(mt[0][0]/m.acc_matrix[0][0]+pONE).double_value() < .1)
    mt = -mt; 

  double err = 0., e;
  int i, j;
  for (i=0; i<2; ++i) {
    for (j=0; j<2; ++j) {
      e = gabs(m.acc_matrix[i][j] - mt[i][j]).double_value(); 
      if (e > err) err = e; 
    }
  }
  return err; 
}

bool operator < (group_element const& a, group_element const& b)
{
  if (!small(log(a.height())-log(b.height())))
    return a.height() < b.height(); 
  if (near(a.pt[0],b.pt[0]) || near(a.pt[0],b.pt[1]))
    return false; 

  return a.w < b.w; 
}

#if 0
static void check_accuracy(pari a[2], pari const& b, double& err)
{
  double e = min(size(gsqrt(gnorml2(b/b[0] - a[0]/a[0][0]))),
		 size(gsqrt(gnorml2(b/b[0] - a[1]/a[1][0]))));
  if (e > err) err = e; 
}
#endif

inline static void add_element(group_element const& e, 
			set<group_element>& Old, 
			set<group_element>& New)
{
  if (Old.insert(e).second) New.insert(e);
}

void compute_products(set<group_element> const& GE, 
		      set<group_element> const& NGE, 
		      set<group_element>& Old,
		      set<group_element>& New)
{
  if (!NGE.size()) return; 
  New.erase(New.begin(), New.end()); 
  set<group_element>::iterator gi, ei, fi;
  group_element p, ginv; 

  match_error = 0.; 

  cout << "Finding products of " << GE.size() << " * " << NGE.size() 
       << " elements\n";

  for (ei=GE.begin(); ei!=GE.end(); ++ei) {
    cout << '.'; cout.flush(); 

    for (gi=NGE.begin(); gi!=NGE.end(); ++gi) {
      ginv = inverse(*gi); 

      add_element(*ei *  *gi, Old, New); 
      add_element(*ei * ginv, Old, New); 
      add_element(*gi  * *ei, Old, New); 
      add_element(ginv * *ei, Old, New); 
    }
  }

  for (ei=NGE.begin(); ei!=NGE.end(); ++ei) {
    cout << '.'; cout.flush(); 

    for (gi=ei; gi!=NGE.end(); ++gi) {
      add_element(*ei * *gi, Old, New); 

      if (ei==gi) continue; // if it was a square

      ginv = inverse(*gi); 

      add_element(*ei * ginv, Old, New); 
      add_element(*gi  * *ei, Old, New); 
      add_element(ginv * *ei, Old, New); 
    }
  }
  cout << endl;
  cout << "Found " << New.size() << " new elements, accuracy " << 
    match_error << endl;
}     

D_domain::D_domain()
{
  pari itet = p_lisexpr("[3,-1,-1,-1;3,-1,1,1;3,1,-1,1;3,1,1,-1]");
  D.initialize(itet);
}

#if 0
static void change_num_gens(GroupPresentation* G, int new_num)
{
  int i, n = fg_get_num_generators(G); 
  if (new_num==n) return; 

  MoebiusTransformation *mtl = NEW_ARRAY(new_num, MoebiusTransformation); 
  O31_matrix *mats = NEW_ARRAY(new_num, O31_matrix); 

  if (new_num < n) n = new_num; 
  for (i=0; i<n; ++i) {
    mtl[i] = fg_Moebius_generators(G)[i];
    mats[i] = G->itsMatrices[i]; 
  }
  my_free_array(fg_Moebius_generators(G));
  my_free_array(G->itsMatrices);
  fg_Moebius_generators(G) = mtl;
  G->itsMatrices = mats;
  G->itsNumGenerators = new_num; 
}
#endif

static int vertex_type(pari const& w)
{
  pari v = gcopy(w); 
  v /= v[0]; 
  v[0] = pZERO; 
  pari n = gnorml2(v);
  if (size(n - pONE) < p_VS::eq_EPS) 
    return 1; // ideal
  if (n < pONE) return 0; // finite
  return 2; 
}

static void count_vertices(polytope<p_VS> const& D, int& nf, int& ni, int& nh)
{
  nf=0, ni=0, nh=0; 
  polytope<p_VS>::vertex_cp vp;
  for (vp = D.VL.begin(); vp!= D.VL.end(); ++vp) {
    switch (vertex_type(vp->V)) {
    case 0:
      ++nf;
      break; 
    case 1:
      ++ni;
      break;
    case 2:
      ++nh;
      break;
    default:
      break;
    }
  }
}

static void to_O31(pari const& m, O31_matrix& mx)
{
  int i, j; 
  for (i=0; i<4; ++i) {
    for (j=0; j<4; ++j) {
      mx(i,j) = m[j][i].double_value(); 
    }
  }
}


static pari reflection(pari const& g)
{
  pari f = g/sqrt(abs(O1n_normsq(g)));
  pari mx = matrix(4,4);
  int i, j; 
  for (i=0; i<4; ++i) {
    for (j=0; j<4; ++j) {
      mx[j][i] = (i==j ? pONE : pZERO) - 
	(j==0 ? -pTWO:pTWO) * f[i] * f[j]; 
    }
  }
  return mx;
}

pari hyperbolic_volume(polytope<p_VS> const& P);

D_domain* compute_domain(GroupPresentation* G, WEPolyhedron** P)
{
  p_VS::eq_EPS = 1e-40;
  p_VS::ne_EPS = 1e-30; 

  cout << "Equality epsilon: " << p_VS::eq_EPS << endl; 

  D_domain* d = new D_domain;
  bool match = d->compute(G);
  cout << PSeq(d->Live, "\n") << endl; 
  if (!match) {
    cout << "\nDomain incomplete\n";
    delete d;
    return 0;
  }

  cout << "Domain has " << d->D.FL.size() << " faces\n";

  d->angle_sums_2pi(false); 

  cout << "Volume is: "; 
  hyperbolic_volume(d->D).print();

  int i, n_gens = fg_get_num_generators(G);
  int nf,ni,nh;
  count_vertices(d->D,nf,ni,nh);

  if (P && nh > 0) {
    pari rfl; 
    polytope<p_VS>::vertex_cp vp;

    cout << "Computing Dirichlet domain with " << nh << " mirrors\n";

    if (*P) free_Dirichlet_domain(*P); 
    O31_matrix* gens = NEW_ARRAY(n_gens+nh, O31_matrix); 

    for (i=0; i<n_gens; ++i) 
      gens[i] = fg_get_generators(G)[i]; 
    for (vp = d->D.VL.begin(); vp!= d->D.VL.end(); ++vp) {
      if (vertex_type(vp->V)!=2) continue;
      to_O31(reflection(vp->V), gens[i]);
      ++i; 
    }
    *P = Dirichlet_from_generators(gens, n_gens+nh, 1e-7, TRUE, FALSE); 
    if (!*P) *P = Dirichlet_from_generators(gens, n_gens+nh, 1e-7,TRUE,FALSE);
  }

  delete d; 
#if 0
  // Compare with Jeff's. 
  WEPolyhedron* dd = 
    Dirichlet_from_generators(fg_get_generators(G), n_gens, 1e-7,FALSE,FALSE);
  if (dd) {
    printf("\n"); 
    print_polyhedron(dd, G, 0xF0, stdout); 
    cout << "Domain has " << dd->num_faces << " faces\n"; 
    my_free(dd); 
  }
#endif
  return 0;
}

pari dir_dist(pari const& F)
{
  return F/gsqrt(gnorml2(extract(F,integer(14))));
}

static void remove_dead(set<group_element>& Live, const face<p_VS>* fp)
{
  set<group_element>::iterator i; 
  for (i=Live.begin(); i!=Live.end(); ++i)
    if (fp == i->f[0] || fp == i->f[1]) break; 
  if (i==Live.end()) return;

  if (fp == i->f[0]) i->kill_face(0); 
  if (fp == i->f[1]) i->kill_face(1); 
  if (i->f[0] || i->f[1]) return; 
  // cout << "Removing " << *i << endl; 
  Live.erase(i); 
}

static bool has_vertex(const face<p_VS>* f, pari const& pt)
{
  list<face<p_VS>::vertex_p>::const_iterator cp;
  for (cp=f->I.begin(); cp!=f->I.end(); ++cp)
    if (near(pt, (*cp)->V)) return true; 
  return false; 
}

static bool faces_match(group_element const& E)
{
  if (!E.f[0] || !E.f[1]) {
    // cout << "null face pointer in live group_element\n";
    return false; 
  }
  if (E.f[0]->I.size() != E.f[1]->I.size()) return false; 
  list<face<p_VS>::vertex_p>::const_iterator cp;
  pari M = SL2C_to_O31(to_pari(E.m));
  for (cp=E.f[1]->I.begin(); cp!=E.f[1]->I.end(); ++cp)
    if (!has_vertex(E.f[0], M * (*cp)->V)) return false; 
  return true; 
}

static double accuracy(set<group_element> const& GE, GroupPresentation* G)
{
  double err=0., e;
  set<group_element>::const_iterator it; 
  for (it=GE.begin(); it!=GE.end(); ++it) {
    e = it->accuracy(G);
    if (e > err) err = e; 
  }
  return err; 
}


bool D_domain::compute(GroupPresentation* G)
{
  int n_gens = fg_get_num_generators(G);
  FGWord w;
  MoebiusTransformation m;
  pari origin = cvector(4), F;
  origin[0] = pONE; 
#if 0
  origin[0] = cosh(pONE/real(100)); 
  origin[1] = sinh(pONE/real(100)); 
#endif
  // Make the generators
  set<group_element> New, Old, NLive;
  int i, cv;
  char gen[2];
  gen[1] = '\0';
  group_element E;
  for (i=0; i<n_gens; ++i) {
    gen[0] = 'a'+i;
    E = group_element(gen, fg_Moebius_generators(G)[i]);
    New.insert(E);
  }


  set<group_element>::iterator ei, li; 
  list<face<p_VS> >::iterator di; 
  bool check=true; 
  bool match=false; 
  int nv, ni, nh=0; 
  polytope<p_VS>::face_cp fp0, fp1;
  list<face<p_VS> > dead; 
  int nmatch; 
  bool full_group; 

  while (true) {

    for (ei=New.begin(); ei!=New.end(); ++ei) {
      if (ei->height() < 1.0001) continue; 
      
      E = *ei;
      E.f[0] = E.f[1] = 0; 
      F = origin - E.pt[0];
      F[0] = -F[0]; 
      if (D.cut(F, &fp0, &dead)!=0) continue;
      E.f[0] = &*fp0; 

      F = origin - E.pt[1];
      F[0] = -F[0]; 
      if (D.cut(F, &fp1, &dead)==0)
	E.f[1] = &*fp1;

      if (D.bad_cut) break; 

      Live.insert(E);
      NLive.insert(E); 

      for (di=dead.begin(); di!=dead.end(); ++di)
	remove_dead(Live, &*di);
      dead.erase(dead.begin(), dead.end()); 
    }

    if (D.bad_cut) break; 
    if (!NLive.size()) break; // no new faces

    // check if we've accidentally dropped to a cover
    full_group = true; 
    for (i=0; i<n_gens; ++i) {
      gen[0] = 'a'+i;
      E = group_element(gen, fg_Moebius_generators(G)[i]);
      if (!group_contains(E)) {
	New.insert(E);
	full_group = false; 
	cout << "Inserting : " << E << endl; 
      }
      if (!full_group) continue; 
    }

    // count hyperideal vertices. 
    if (nh) {
      count_vertices(D,nv,ni,nh); 
      cout << "Polytope has " << nh << " hyperideal vertices\n";
    }

    // check if faces match
    if (nh==0) {
      nmatch = 0; 
      for (li=Live.begin(); li!=Live.end(); ++li) 
	if (faces_match(*li)) ++nmatch; 
      cout << nmatch << '/' << Live.size() << " faces match\n"; 
      if (nmatch == Live.size()) {
	match = true; 
	break; 
      }
    }

    // compute more products
    compute_products(Live, NLive, Old, New);
    NLive.erase(NLive.begin(), NLive.end()); 
  }

  return match; 
}

// The following comes from the formula given in
// Ushijima A., "A Volume Formula for Generalized Hyperbolic Tetrahedra."
// arXiv:math.GT/0309216 v2

pari hyperbolic_tet_volume(pari const& cosangles)
{
  pari cA = cosangles[0];
  pari cB = cosangles[1];
  pari cC = cosangles[2];
  pari cD = cosangles[3];
  pari cE = cosangles[4];
  pari cF = cosangles[5];
  pari sA = sqrt(1-cA*cA);
  pari sB = sqrt(1-cB*cB);
  pari sC = sqrt(1-cC*cC);
  pari sD = sqrt(1-cD*cD);
  pari sE = sqrt(1-cE*cE);
  pari sF = sqrt(1-cF*cF);
  pari a = cA + sA*pI; // change to complex(cA,sA);
  pari b = cB + sB*pI;
  pari c = cC + sC*pI;
  pari d = cD + sD*pI;
  pari e = cE + sE*pI;
  pari f = cF + sF*pI;
  pari G = matrix(4,4);
  G[0][0] =pONE; G[1][0] = -cA; G[2][0] = -cB; G[3][0] = -cF;
  G[0][1] = -cA; G[1][1] =pONE; G[2][1] = -cC; G[3][1] = -cE;
  G[0][2] = -cB; G[1][2] = -cC; G[2][2] =pONE; G[3][2] = -cD;
  G[0][3] = -cF; G[1][3] = -cE; G[2][3] = -cD; G[3][3] =pONE;
  pari sd = sqrt(det(G));
  pari sp = sA*sD + sB*sE + sC*sF;
  pari den = a*d+b*e+c*f+a*b*f+a*c*e+b*c*d+d*e*f+a*b*c*d*e*f;
  pari z1 = pTWO*(sp - sd)/den;
  pari z2 = pTWO*(sp + sd)/den;
  pari U1 = pHALF*(dilog(z1)+dilog(a*b*d*e*z1)+
		   dilog(a*c*d*f*z1)+dilog(b*c*e*f*z1)
		   -dilog(-a*b*c*z1)-dilog(-a*e*f*z1)
		   -dilog(-b*d*f*z1)-dilog(-c*d*e*z1));
  pari U2 = pHALF*(dilog(z2)+dilog(a*b*d*e*z2)+
		   dilog(a*c*d*f*z2)+dilog(b*c*e*f*z2)
		   -dilog(-a*b*c*z2)-dilog(-a*e*f*z2)
		   -dilog(-b*d*f*z2)-dilog(-c*d*e*z2));
  return pHALF*gimag(U1-U2);
}

pari Lobachevsky(pari const& theta)
{
  return pHALF*gimag(dilog(exp(pTWO*pI*theta)));
}

// A tetrahedron [a,b,c,d] is birectangular if the
// dihedral angles are right angles at [a,c], [b,c] and [b,d]. 
// Writing aa, bb, cc, dd for the (opposite) face normals
// it follows that <bb,dd> == <aa,dd> == <aa,cc> == 0,
// and the tetrahedron is in fact parametrized by 
// <aa,bb>, <bb,cc> and <cc,dd>.

pari birectangular_tetrahedron_volume(pari const& m0)
{
  /*
   *  Compute the volume of the birectangular tetrahedron with vertices
   *  a, b, c and d, the columns of m0, using the method found in 
   *  section 4.3 of
   *
   *      E. B. Vinberg, Ob'emy neevklidovykh mnogogrannikov,
   *          Uspekhi Matematicheskix Nauk, May(?) 1993, 17-46.
   *
   *  Our a, b, c and d correspond to Vinberg's A, B, C and D, as shown
   *  in his Figure 9.  We need to compute the dual basis {aa, bb, cc, dd}
   *  defined by <aa, a> = 1, <aa, b> = <aa, c> = <aa, d> = 0, etc.
   *  Let m be the matrix whose rows are the vectors {a, b, c, d}, but
   *  with the entries in the first column negated to account for the
   *  indefinite inner product, and let mm be the matrix whose columns
   *  are the vectors {aa, bb, cc, dd}.  Then (m)(mm) = identity, by the
   *  definition of the dual basis.
   */

  pari m = gtrans(m0); 

  // If det==0 volume is zero. 
  pari ds = psign(det(m)); 
  if (ds == pZERO) return pZERO; 

  m[0] = -m[0];

  pari mm = ginv(m);

  pari aa=mm[0], bb=mm[1], cc=mm[2], dd=mm[3];

  /*
   *  Any pair of dual vectors lies in a positive definite 2-plane
   *  in E^(3,1).  Normalize them to have length one, so we can use
   *  their dot products to compute the dihedral angles.
   */
  aa /= sqrt(O1n_normsq(aa));
  bb /= sqrt(O1n_normsq(bb));
  cc /= sqrt(O1n_normsq(cc));
  dd /= sqrt(O1n_normsq(dd));

  /*
   *  Compute the angles alpha, beta and gamma, as shown in
   *  Vinberg's Figure 9.
   */
  pari alpha   = pPI - acos(O1n_innerprod(aa, bb));
  pari beta    = pPI - acos(O1n_innerprod(bb, cc));
  pari gamma   = pPI - acos(O1n_innerprod(cc, dd));

  /*
   *  Compute big_delta and delta, as in
   *  Vinberg's Sections 4.2 and 4.3.
   */
  pari big_delta = gsqr(sin(alpha)*sin(gamma)) - gsqr(cos(beta));

  if (big_delta > pZERO) {
    cout << "big_delta > 0 in birectangular_tetrahedron_volume\n";
    big_delta.print(); 
    return pZERO;
  }

  pari delta = atan(sqrt(-big_delta)/(cos(alpha) * cos(gamma)));

  pari tetrahedron_volume = pQUARTER * (
    Lobachevsky(alpha + delta)
    - Lobachevsky(alpha - delta)
    + Lobachevsky(gamma + delta)
    - Lobachevsky(gamma - delta)
    - Lobachevsky(pPI_OVER_2 - beta + delta)
    + Lobachevsky(pPI_OVER_2 - beta - delta)
    + pTWO * Lobachevsky(pPI_OVER_2 - delta)
    );

  return ds*tetrahedron_volume; // ds makes volume negative for negatively oriented tetrahedra. 
}

pari orientation(polytope<p_VS>::edge_p const& ep, 
		 polytope<p_VS>::face_cp const& fp)
{
  // orientation is sign of determinant of
  // a matrix whose columns are [origin, ep->a, ep->b, c]
  // where c is a point on the face other than a or b. 

  pari m = matrix(4,4);
  m[0][0]=pONE;
  m[1] = ep->a->V;
  m[2] = ep->b->V;

  // now find a point c on face, not a or b. 
  polytope<p_VS>::IV_p vp;
  for (vp=fp->I.begin(); vp!=fp->I.end(); ++vp)
    if (*vp != ep->a && *vp != ep->b) break;
  if (vp==fp->I.end()) {
    cout << "couldn't find a vertex not on given edge in orientation().\n";
    return 0;
  }

  m[3] = (*vp)->V; // this is c. 

  return psign(det(m));
}

pari hyperbolic_volume(polytope<p_VS> const& P)
{
  polytope<p_VS>::edge_p ep;
  face_list::const_iterator fa, fb;
  pari FA, FB, FBN, CFA, CFB, CE, Or, Vol, T, Eor;

  Or = cvector(4); Or[0] = pONE;

  T = matrix(4,4);
  T[3] = Or;

  // for each edge there are 4 birectangular tetrahedra
  for (ep = P.EL.begin(); ep != P.EL.end(); ++ep) {
    find_faces(P.FL, *ep, fa, fb);

    FA = fa->F/sqrt(O1n_normsq(fa->F)); FA[0] = -FA[0]; 
    FB = fb->F/sqrt(O1n_normsq(fb->F)); FB[0] = -FB[0]; 
    
    CFA = Or + FA[0]*FA;
    CFB = Or + FB[0]*FB;

    FBN = FB - O1n_innerprod(FB,FA)*FA;
    CE  = CFA - (O1n_innerprod(CFA, FBN)/O1n_normsq(FBN))*FBN;

    // At this point CFA, CFB, CE should be 
    // the (non-normalized) projections of the origin
    // onto the faces (normal to) FA, FB and their
    // intersection. 

#if 0
    // In particular the inner products <CE,FA> and <CE,FB>
    // should both be zero. 
    cout << "CE on FA & FB\n";
    O1n_innerprod(CE,FA).print();
    O1n_innerprod(CE,FB).print();

    cout << "Va on FA & FB\n"; 
    O1n_innerprod(ep->a->V,FA).print();
    O1n_innerprod(ep->a->V,FB).print();

    cout << "innerprod Or,CF,V\n";
    O1n_innerprod(CFA-ep->a->V, Or-CFA).print(); 

    cout << endl; 
#endif

    Eor = orientation(ep, fa);

    // Now add up the volume contributions for the tetrahedra
    // [ep->a->V, CE, CFA, Or], [ep->b->V, CE, CFA, Or],
    // [ep->a->V, CE, CFB, Or], [ep->b->V, CE, CFB, Or].

    T[1] = CE;
    T[0] = ep->a->V; 
    T[2] = CFA; 
    Vol -= Eor * birectangular_tetrahedron_volume(T);

    T[0] = ep->b->V; 
    Vol += Eor * birectangular_tetrahedron_volume(T);

    T[2] = CFB;
    Vol -= Eor * birectangular_tetrahedron_volume(T);

    T[0] = ep->a->V;
    Vol += Eor * birectangular_tetrahedron_volume(T);
  }

  return Vol;
}

